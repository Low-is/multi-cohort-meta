library(jsonlite)
source("meta/functions/expr_mtx.R")

# -------------------------
# LOAD STUDY LISTS
# -------------------------
dna_studies <- fromJSON("miner/outputs/dna_gse_ids.json")
rna_studies <- fromJSON("miner/outputs/rna_gse_ids.json")

# -------------------------
# RUN PIPELINE
# -------------------------
dna_exprs <- generate_exprs_mtx(
  DNA = TRUE,
  dna_studies = dna_studies
)

rna_exprs <- generate_exprs_mtx(
  RNA = TRUE,
  rna_studies = rna_studies
)
