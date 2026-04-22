library(jsonlite)
library(GEOquery)
library(Biobase)
library(limma)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

# -------------------------
# LOAD FUNCTIONS
# -------------------------
message("🚀 Loading expression matrix functions...")
source("meta/functions/expr_mtx.R")

# -------------------------
# LOAD STUDY LISTS
# -------------------------
message("📥 Loading study ID JSON files...")

dna_studies <- fromJSON("miner/outputs/dna_gse_ids.json")
rna_studies <- fromJSON("miner/outputs/rna_gse_ids.json")

message(sprintf("✔ Loaded %d DNA studies", length(dna_studies)))
message(sprintf("✔ Loaded %d RNA studies", length(rna_studies)))


# -------------------------
# RUN PIPELINE
# -------------------------
message("\n🧬 Starting DNA (microarray) processing...")

dna_exprs <- generate_exprs_mtx(
  DNA = TRUE,
  dna_studies = dna_studies
)

message("\n🧬 Starting RNA-seq processing...")

rna_exprs <- generate_exprs_mtx(
  RNA = TRUE,
  rna_studies = rna_studies
)


# -------------------------
# FINAL SUMMARY
# -------------------------
message("\n🎉 Pipeline finished successfully")
message(sprintf("Total DNA datasets processed: %d", length(dna_exprs)))
message(sprintf("Total RNA datasets processed: %d", length(rna_exprs)))
