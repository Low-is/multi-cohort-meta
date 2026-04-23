# Loading functions
source("meta/functions/get_pData.R")

# Load config
library(yaml)
config <- yaml::read_yaml("meta/config/config.yaml")

# Define search patterns
case_patterns <- config$case_patterns
control_patterns <- config$control_patterns

# Loading pData
dna_pData <- readRDS("meta/pdata/dna_pData.rds")
rna_pData <- readRDS("meta/pdata/rna_pData.rds")

# Apply condition detection + labeling
dna_pData_cond <- apply_codition_to_list(
  dna_pData,
  case_patterns,
  control_patterns
)

rna_pData_cond <- apply_condition_to_list(
  rna_pData,
  case_patterns,
  control_patterns
)

# Save results
saveRDS(dna_pData_cond, "meta/pdata/dna_pData_with_condition.rds")
saveRDS(rna_pData_cond, "meta/pdata/rna_pData_with_condition.rds")
# Saving as csv
write.csv(dna_pData_cond, "meta/pdata/dna_pData_with_condition.csv")
write.csv(rna_pData_cond, "meta/pdata/rna_pData_with_condition.csv")

message("✅ Condition column successfully added to all studies")
