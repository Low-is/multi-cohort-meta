#library(jsonlite)
#library(dplyr)
#library(GEOquery)
#library(Biobase)
#library(limma)
#library(data.table)
#library(AnnotationDbi)
#library(org.Hs.eg.db)
#library(hta20transcriptcluster.db)


# -------------------------
# LOADING PDATA
# -------------------------
message("🧠 Loading pData...")

dna_pData <- get_pData(studies = dna_studies)
rna_pData <- get_pData(studies = rna_studies)

message("\n pData successfully loaded!")


# -------------------------
# SAVING RESULTS
# -------------------------
message("Saving pData...")
saveRDS(dna_pData, "meta/pdata/dna_pData.rds")
saveRDS(rna_pData, "meta/pdata/rna_pData.rds")
