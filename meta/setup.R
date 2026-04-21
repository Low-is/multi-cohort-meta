# =========================================================
# SETUP SCRIPT FOR MULTI-COHORT META ANALYSIS PIPELINE
# (Run ONCE after cloning repo)
# =========================================================

# -----------------------
# 1. Install renv if missing
# -----------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# -----------------------
# 2. Initialize renv project
# -----------------------
# Use bare = TRUE so it doesn't try to snapshot immediately
renv::init(bare = TRUE)

# -----------------------
# 3. Install CRAN dependencies
# -----------------------
install.packages(c(
  "jsonlite",
  "data.table"
))

# -----------------------
# 4. Install Bioconductor manager (if needed)
# -----------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# -----------------------
# 5. Install Bioconductor dependencies
# -----------------------
BiocManager::install(c(
  "GEOquery",
  "Biobase",
  "limma",
  "AnnotationDbi",
  "org.Hs.eg.db"
), ask = FALSE)

# -----------------------
# 6. Save locked environment
# -----------------------
renv::snapshot()

# -----------------------
# DONE
# -----------------------
message("✔ Setup complete. Environment is now reproducible.")
