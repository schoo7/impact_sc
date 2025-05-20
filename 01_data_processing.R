# IMPACT-sc Script: 01_data_processing.R
# Purpose: Load raw data, create Seurat object, QC, and initial processing.
####
# --- Libraries ---
# Ensure all these libraries are installed as per README.md
library(AnnotationDbi)
library(CARD) # Though not directly used in this script, loaded for consistency if running all
library(celldex) # For later modules
library(decoupleR) # For later modules
library(dplyr)
library(ensembldb) # For later modules
library(ggplot2)
library(ggpubr)
library(homologene) # For later modules
library(infercnv) # For later modules
library(Matrix)
library(msigdbr) # For later modules
library(patchwork)
library(presto) 
library(reticulate) 
library(rjags) # For CARD
library(scDblFinder) # For later modules
library(scater)
library(scRNAtoolVis) 
library(SCpubr) 
library(Seurat)
library(SeuratDisk) 
# library(SeuratExtend) # Ensure this is available if using its specific functions
library(SingleCellExperiment)
library(SingleR) # For later modules
library(SpatialExperiment) # For later modules
library(scran)
library(stringr)
library(tibble)
library(tidyr)
library(UCell) # For later modules
library(viridis)
library(reshape2)
library(pheatmap)

# --- IMPACT-sc Script Parameters ---
species <- "human" # MODIFY THIS PARAMETER AS NEEDED ("human" or "mouse")
set.seed(1234)
base_output_path <- "../output/" # Assuming scripts are in a 'scripts' subdirectory
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}

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

obj_rds_path <- "../ori.RDS" # Path relative to the script's location in 'scripts/'
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
mt_pattern <- if (species == "human") "^MT-" else if (species == "mouse") "^mt-" else "^MT-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
message(paste("Calculated mitochondrial percentage using pattern:", mt_pattern))

qc_plot_path <- file.path(base_output_path, "01_qc_violin_plot.png") # Numbered output
png(qc_plot_path, width=10, height=6, units="in", res=300)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
dev.off()
message("QC violin plot saved to: ", qc_plot_path)

# Filter low-quality cells
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
message(paste("Cells after QC filtering:", ncol(obj)))

# Standard single-cell processing
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

module1_output_path <- file.path(base_output_path, "01_module1_processed.RDS") # Numbered output
tryCatch({
  saveRDS(obj, module1_output_path)
  message("Module 1 output saved to: ", module1_output_path)
}, error = function(e) {
  warning(paste("Error saving Module 1 output:", e$message))
})

message("Finished Module 1: Data Processing")
