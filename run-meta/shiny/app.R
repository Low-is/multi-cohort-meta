library(shiny)
library(shinythemes)
library(shinydashboard)
library(DT)
library(dplyr)
library(GEOquery)
library(limma)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(data.table)
library(patchwork)
library(stringr)
library(Biobase)
library(purrr)
library(GEOmetadb)
library(DBI)  # For database connectivity
library(RSQLite)
library(multtest)
library(meta)
library(rmeta)
library(caret)
library(org.Hs.eg.db)
library(AnnotationDbi)


# Source your custom functions
source("\\\\ifs.win.uthscsa.edu/M1509-AhujaS/MainShare/Lois/Lois_Local/Dr_M/Projects/Meta-Multicohort-Pipeline/functions/meta_analysis_functions.R")
source("\\\\ifs.win.uthscsa.edu/M1509-AhujaS/MainShare/Lois/Lois_Local/Dr_M/Projects/Meta-Multicohort-Pipeline/functions/Data_Processing.R")



source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)


