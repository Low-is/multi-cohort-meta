#### Meta-analysis function ####
meta_results <- function(list_of_studies) {

list_of_es <- lapply(list_of_studies, function( study ) {
  effect.sizes(study)
  })
  
  summary <- combine.effect.sizes(list_of_es)
  
  summary$g <- summary$g[!apply(is.na(summary$g) | is.infinite(summary$g), 1, any), ]
  summary$se.g <- summary$se.g[!apply(is.na(summary$se.g) | is.infinite(summary$se.g), 1, any), ]
  g_genes <- rownames(summary$g)
  
  summary$pooled.estimates$genes <- rownames(summary$pooled.estimates)
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(genes %in% g_genes)
  
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(!apply(is.na(.) | is.infinite(as.matrix(.)), 1, any))
  
  summary$pooled.estimates <- as.data.frame(lapply(summary$pooled.estimates, function(x) ifelse(x == 0, 1e-200, x)))
  rownames(summary$pooled.estimates) <- g_genes
  
  g <- summary$g
  se.g <- summary$se.g
  
  pool    <- summary$pooled.estimates[, "summary" ]
  names(pool) <- rownames(g)
  
  se.pool <- summary$pooled.estimates[, "se.summary" ]
  names(se.pool) <- rownames(g)
  
  x.label <- "Standardized Mean Difference (log2 scale)"
  
  consistent_genes <- c()
  
  for (gene in rownames(g)) {
    g_gene <- g[gene, ]
    if ((all(g_gene > 0) & pool[gene] > 0.5) | (all(g_gene < 0) & pool[gene] < -0.5)) {
      consistent_genes <- c(consistent_genes, gene)
    }
  }
  
  for (gene in consistent_genes) {
    g_gene <- g[gene, ]
    se.g_gene <- se.g[gene, ]
    study_names <- gsub("_g", "", names(g_gene))
    
    pdf(file = paste0(gene, ".pdf"), height = 3.5, width = 3.5)
    par(cex = 0.65)
    metaplot(g_gene, se.g_gene,
             labels = study_names,
             summn = pool[gene],
             sumse = se.pool[gene],
             sumnn = 1/se.pool[gene]^2,
             summlabel = "Summary effect",
             xlab = x.label,
             ylab = "",
             main = bquote(italic(.(gene))),
             colors = meta.colors(box = "violetred", lines = "plum", summary = "mediumpurple", 
                                  text = "black", axes = "black", zero = "black"),  
             boxsize = 1,
             lty.random = 1,  
             lwd.random = 2,   
             zero = 0,  
             col.zero = "black", 
             lty.zero = 3)
    dev.off()
  }
  
  return(list(summary = summary,
              genes = consistent_genes))
}

####
####

#### Function to store meta-analysis results ####
generate_list_for_meta_analysis <- function(DNA = FALSE, RNA = FALSE, dna_mtx = list(), rna_mtx = list(), pData = list(), study = list()) {
  
  #### Creating empty list to store meta-analysis results ####
  list_of_studies <- list()
  
  #### Preparing data for meta-analysis ####
  if(DNA && !RNA) {
    common_genes <- find_common_genes(DNA = TRUE, DNA_mtx = dna_mtx)
    for(i in 1:length(study)) {
      
      name <- names(study)[i]
      mtx <- dna_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
      
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  #### Preparing data for meta-analysis ####
  if(RNA && !DNA) {
    common_genes <- find_common_genes(RNA = TRUE, RNA_mtx = RNA_matrices$norm_counts)
    for(i in 1:length(study)) {
      name <- names(study)[i]
      mtx <- rna_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
     
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  
  if(DNA && RNA) {
    common_genes <- find_common_genes(DNA = TRUE, DNA_mtx = DNA_matrices, RNA = TRUE, RNA_mtx = RNA_matrices$norm_counts)
    for(i in 1:length(study)) {
      
      name <- names(study)[i]
      list_of_mtx <- c(DNA_mtx, RNA_mtx)
      mtx <- list_of_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
      
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  return(list(studies = list_of_studies,
              results = meta_results(list_of_studies)))
}

####
####

#### Return pvalues and qvalues ####
get_p_q_values <- function(meta_results) {
  
  cleaned_studies <- list(
    studies = lapply(meta_results$studies, function(study) {
      
      valid_rows <- complete.cases(study$class, study$pheno)
      expr_clean <- study$expr[, valid_rows]
      
      valid_genes <- apply(expr_clean, 1, function(gene_values) {
        all(!is.na(gene_values) & !is.nan(gene_values) & !is.infinite(gene_values))
      })
      expr_clean <- expr_clean[valid_genes, ]
      
      study_cleaned <- list(
        expr = expr_clean,
        pheno = study$pheno[valid_rows],
        keys = study$keys[valid_genes],
        class = study$class[valid_rows]
      )
      return(study_cleaned)
    }),
    
    results = meta_results$results
  )
  
  list.of.sigs <- lapply(cleaned_studies$studies, function( study ) {
    fisher_pval <- ttest.Pvalues(study)
    genes <- cleaned_studies$results$genes
    fisher_pval <- fisher_pval %>% dplyr::filter(keys%in%genes)
  })
  
  output.Fisher <- as.data.frame(sum.of.logs(list.of.sigs)) # For fisher p-value
  adj_pval <- adjust.fisher(output.Fisher=output.Fisher)[, c("F.Qval.up", "F.Qval.down")]
  
  genes <- cleaned_studies$results$genes
  pool <- cleaned_studies$results$summary$pooled.estimates[genes, "summary"]
  se.pool <- cleaned_studies$results$summary$pooled.estimates[genes, "se.summary" ]
  
  es_pval <- numeric()
  es_qval <- numeric()
  
  for(i in 1:nrow(output.Fisher)) {
    if(output.Fisher$F.pval.up[i] <= 0.05) {
      es_pval <- c(es_pval, output.Fisher$F.pval.up[i])
    }
    if(output.Fisher$F.pval.down[i] <= 0.05) {
      es_pval <- c(es_pval, output.Fisher$F.pval.down[i])
    }
  }
  
  for(i in 1:nrow(adj_pval)) {
    if(adj_pval$F.Qval.up[i] <= 0.05) {
      es_qval <- c(es_qval, adj_pval$F.Qval.up[i])
    }
    if(adj_pval$F.Qval.down[i] <= 0.05) {
      es_qval <- c(es_qval, adj_pval$F.Qval.down[i])
    }
  }
  
  sig <- data.frame(Genes = genes,
                    Pooled_es = pool,
                    Pooled_es_se = se.pool,
                    es_pval = es_pval,
                    es_qval = es_qval)
  return(sig)
}
