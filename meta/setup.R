# =========================================================
# SETUP SCRIPT FOR MULTI-COHORT META ANALYSIS PIPELINE
# (Run ONCE after cloning repo)
# =========================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

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
  "data.table",
  "dplyr",
  "yaml"
))

untar("meta/COCONUT_1.0.2.tar.gz", exdir = "meta/COCONUT")
#list.files("meta/COCONUT", recursive = TRUE)
install.packages("meta/COCONUT/COCONUT",
                repos = NULL,
                type = "source")

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
  "org.Hs.eg.db",
  "hta20transcriptcluster.db",
  "DESeq2"
), ask = FALSE)


# -----------------------
# 6. Install GEOmetadb + download SQLite
# -----------------------
if (!requireNamespace("GEOmetadb", quietly = TRUE)) {
  BiocManager::install("GEOmetadb", ask = FALSE)
}

message("📦 Downloading GEOmetadb SQLite...")

# Create a consistent directory inside project
dir.create("meta/db", recursive = TRUE, showWarnings = FALSE)

options(timeout = 600)
options(download.file.method = "libcurl")
# Download database
sqlite_path <- GEOmetadb::getSQLiteFile(destdir = "meta/db")

message("✔ GEOmetadb SQLite saved at: ", sqlite_path)


# -----------------------
# 7. Save locked environment
# -----------------------
renv::snapshot()

# -----------------------
# DONE
# -----------------------
message("✔ Setup complete. Environment is now reproducible.")
