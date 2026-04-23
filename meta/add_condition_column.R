# Loading functions
source("../functions/get_pData.R")

# Load config
library(yaml)
config <- yaml::read_yaml("../config/config.yaml")

# Define search patterns
case_patterns <- config$case_patterns
control_patterns <- config$control_patterns

# Loading pData
dna_pData <- readRDS("../pdata/dna_pData.rds")
rna_pData <- readRDS("../pdata/rna_pData.rds")
