library(jsonlite)
library(dplyr)
library(GEOquery)
library(Biobase)
library(limma)
library(data.table)
source("meta/functions/get_pData.R")


# -------------------------
# LOAD STUDY LISTS
# -------------------------
message("📥 Loading study ID JSON files...")

dna_studies <- fromJSON("miner/outputs/dna_gse_ids.json")
rna_studies <- fromJSON("miner/outputs/rna_gse_ids.json")

message(sprintf("✔ Loaded %d DNA studies", length(dna_studies)))
message(sprintf("✔ Loaded %d RNA studies", length(rna_studies)))


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
