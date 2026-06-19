library(yaml)
library(jsonlite)
library(COCONUT)
library(DESeq2)
library(data.table)
library(dplyr)
library(rmeta)
library(forestplot)
source("run-meta/functions/meta_analysis_functions.R")
source("run-meta/functions/Sepsis_MC_analysis_functions.R")
source("meta/functions/expr_mtx.R")
source("meta/functions/filter_data.R")


# Loading config file
message("Loading config file...")
config <- yaml::read_yaml("run-meta/config/config.yaml")
message("Config file loaded!")

# Loading named list of studies
message("Loading list of studies...")
dna_studies <- jsonlite::fromJSON(config$analysis$input$dna_gse_file)
rna_studies <- jsonlite::fromJSON(config$analysis$input$rna_gse_file)
message("Studies loaded!")

# Loading expression matrices
message("Loading expression matrices...")
dna_matrices <-  readRDS("meta/matrices/dna_matrices.rds")
dna_matrices <- dna_matrices[!sapply(dna_matrices, is.null)]

# Function to get rid of methylation experiements
is_meth <- function(genes) {
  mean(grepl("^cg", genes, ignore.case = TRUE)) > 0.05
}

dna_matrices <- dna_matrices[!sapply(dna_matrices, function(m) is_meth(rownames(m)))]

#rna_matrices <- readRDS("meta/matrices/rna_matrices.rds") # the str is a list of 2: x$expr and x$pData
#rna_matrices <- rna_matrices[!sapply(rna_matrices, is.null)] # Removing NULL single cell datasets
rna_matrices <- readRDS("meta/matrices/norm_rna_mtxs.rds")
message("Matrices loaded!")

                                     
message("Loading pData...")
dna_pData <- readRDS("meta/pdata/dna_pData_with_condition.rds")
rna_pData <- readRDS("meta/pdata/rna_pData_with_condition.rds")
message("pData loaded!")


# ----------------------------
# NORMALIZATION (GLOBAL SAFE)
# ----------------------------
normalize <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}


# Need to add code that filters matrices to match dimmensions of pData         
# Filtering list
dna_matrices <- dna_matrices[!is.na(names(dna_matrices))]
dna_pData <- dna_pData[names(dna_matrices)]

dna_matrices <- mapply(function(mat, pd) {
  if (is.null(mat) || is.null(pd)) return(NULL)

  common_gsm <- intersect(colnames(mat), pd$gsm)
  mat <- mat[, common_gsm, drop = FALSE]

  return(mat)
},
                       dna_matrices, dna_pData, SIMPLIFY = FALSE)

                                     

dna_matrices <- lapply(dna_matrices, function(x) {
  if (is.null(x)) return(NULL)
  x[complete.cases(x), , drop = FALSE]
})

valid_studies <- sapply(dna_pData, function(pd) {
  tab <- table(pd$condition)

  length(tab) >= 2 %% all(tab >= 2)
})

valid_studies[is.na(valid_studies)] <- FALSE


dna_matrices <- dna_matrices[valid_studies]
dna_pData <- dna_pData[valid_studies]
                                                                                        
# Find common genes across all studies being used for meta-analysis
message("Searching for common genes...")
common_genes <- find_common_genes(DNA = config$analysis$modalities$DNA,
                                  RNA = config$analysis$modalities$RNA,
                                  list_of_dna_mtx = dna_matrices,
                                  list_of_rna_mtx = rna_matrices,
                                  use_DEG = config$analysis$use_DEG
                                 )
message(sprintf("%d common genes detected!", length(common_genes)))


message("Starting meta-analysis...")
                                     
combined_pData <- c(dna_pData, rna_pData)

study <- names(c(dna_matrices, rna_matrices))
                                     
meta_res <- generate_list_for_meta_analysis(
  DNA = config$analysis$modalities$DNA,
  RNA = config$analysis$modalities$RNA,
  list_of_dna_mtx = dna_matrices,
  list_of_rna_mtx = rna_matrices,
  list_of_pData = combined_pData,
  study = study,
  common_genes = common_genes
)

message("Meta-analysis completed!")
message("Saving meta results in run-meta/output folder...")
consistent_genes <- meta_res$meta$consistent_genes
meta_res_df <- meta_res$meta$pooled.estimates
meta_res_df <- meta_res_df[consistent_genes, ] 

which(meta_res_df$FDR <= 0.05)
                                     
write.csv(meta_res_df, "run-meta/output/meta_results.csv")
saveRDS(consistent_genes, "run-meta/output/consistent_genes.rds")

                                
message("Plotting foresplots...")
fp_data <- meta_res_df %>%
  mutate(
    logOR = summary,
    SE    = se.summary,
    OR    = exp(logOR),
    lower = exp(logOR - 1.96 * SE),
    upper = exp(logOR + 1.96 * SE),
    FDR = round(FDR, 3)
  )

fp_data <- fp_data %>%
  mutate(
    Gene = rownames(meta_res_df),
    CI = sprintf("%.2f - %.2f", lower, upper),
    OR_txt = sprintf("%.2f", OR)
  ) %>%
arrange(desc(OR), Gene)

label_mat <- cbind(
  as.character(fp_data$Gene),
  fp_data$OR_txt,
  fp_data$CI,
  fp_data$FDR
)
colnames(label_mat) <- c("Gene", "OR", "CI", "FDR")

png("run-meta/output/forestplot.png",
    width = 10,
    height = 14,
    units = "in",
    res = 800)

forestplot(
  labeltext = label_mat,
  mean  = fp_data$OR,
  lower = fp_data$lower,
  upper = fp_data$upper,
  
  boxsize = config$forestplot$boxsize,
  clip = as.numeric(config$forestplot$clip),
  xticks = as.numeric(config$forestplot$xticks),
  
  xlog = config$forestplot$xlog,
  zero = config$forestplot$zero,
  
  graphwidth = unit(config$forestplot$graphwidth, "npc"),
  
  colgap = unit(config$forestplot$colgap, "mm"),
  lineheight = unit(config$forestplot$lineheight, "mm"),
  
  col = fpColors(
    box = config$forestplot$colors$box,
    line = config$forestplot$colors$line,
    summary = config$forestplot$colors$summary
  ),
  txt_gp = fpTxtGp(
    ticks = gpar(cex = config$forestplot$text_gp$ticks$cex),
    xlab = gpar(cex = config$forestplot$text_gp$xlab$cex)
  )
) |>
  fp_add_header(
    Gene   = c("", "Gene"),
    OR     = c("", "OR"),
    CI     = c("", "95% CI"),
    FDR = c("", "FDR")
  ) |>
  fp_set_zebra_style(config$forestplot$zebrastyle) |>
  fp_set_favors(low = config$forestplot$set_favors$low,
                high = config$forestplot$set_favors$high,
                txt_gp = gpar(cex = config$forestplot$set_favors$text_gp),
                arrows = config$forestplot$set_favors$arrows)

dev.off()
message("Will find saved forestplots in run-meta/output folder")
