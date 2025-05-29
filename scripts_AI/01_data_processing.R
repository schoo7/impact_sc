# IMPACT-sc Script: 01_data_processing.R
# Purpose: Load raw data, create Seurat object, QC, and initial processing.
####
# --- Libraries ---
# Ensure all these libraries are installed as per README.md

library(AnnotationDbi)
library(CARD) 
library(celldex) 
library(decoupleR) 
library(dplyr)
library(ensembldb) 
library(ggplot2)
library(ggpubr)
library(homologene) 
library(Matrix)
library(msigdbr) 
library(patchwork)
library(reticulate) 
library(scDblFinder) 
library(scater)
library(scRNAtoolVis) 
library(SCpubr) 
library(Seurat)
library(SeuratDisk) 
library(SeuratExtend) 
library(SingleCellExperiment)
library(SingleR) 
library(SpatialExperiment) 
library(scran)
library(stringr)
library(tibble)
library(tidyr)
library(UCell) 
library(viridis)
library(reshape2)
library(pheatmap)
# --- IMPACT-sc Script Parameters (READ FROM ENVIRONMENT VARIABLES) ---
# Species: Read from environment variable, default to "human"
species <- Sys.getenv("IMPACT_SC_SPECIES", "human")
set.seed(1234)
# Base Output Path: Read from environment variable
base_output_path <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", "../output/") # Default to "../output/" if not set
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}
message(paste("Output directory set to:", base_output_path))

# Load species-specific annotation database
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    library(org.Hs.eg.db)
    org_db_object <- org.Hs.eg.db
    message("Loaded org.Hs.eg.db for human.")
  } else {
    stop("org.Hs.eg.db not installed. Please install it via BiocManager.")
  }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    org_db_object <- org.Mm.eg.db
    message("Loaded org.Mm.eg.db for mouse.")
  } else {
    stop("org.Mm.eg.db not installed. Please install it via BiocManager.")
  }
} else {
  stop("Species not supported by IMPACT-sc. Please choose 'human' or 'mouse'.")
}
# --- End Script Parameters ---

#### Module 1: Data Process ####
message("Starting Module 1: Data Processing")

# Input: Expects 'ori.RDS' in the parent directory or adjust path.
# Output: Saves 'module1_processed.RDS' in the 'output/' directory.

# Read input RDS path from environment variable, default to ../ori.RDS if not set.
obj_rds_path <- Sys.getenv("IMPACT_SC_INPUT_DATA_PATH", "../ori.RDS") 
obj <- tryCatch({
  readRDS(obj_rds_path)
}, error = function(e) {
  stop(paste("Error loading ori.RDS from", obj_rds_path, ":", e$message))
  return(NULL)
})
if (is.null(obj)) stop("Failed to load input data for IMPACT-sc.")

if (length(unique(obj$orig.ident)) > 1) {
    message("Multiple samples detected by 'orig.ident'.")
}
if (inherits(obj[["RNA"]], "list") || (is.list(obj@assays$RNA) && inherits(obj@assays$RNA[[1]], "Assay"))) {
    message("RNA assay is split. Joining layers for initial QC and processing.")
    obj <- JoinLayers(obj, assay="RNA")
}

# Perform basic quality control
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 )
message(paste("Cells after QC filtering:", ncol(obj)))

# Standard single-cell processing
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

module1_output_path <- file.path(base_output_path, "01_module1_processed.RDS") 
tryCatch({
  saveRDS(obj, module1_output_path)
  message("Module 1 output saved to: ", module1_output_path)
}, error = function(e) {
  warning(paste("Error saving Module 1 output:", e$message))
})

message("Finished Module 1: Data Processing")