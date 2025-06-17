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

# --- NEW: Optional Processing Parameters ---
# Read user choices from environment variables set by the python orchestrator
remove_doublets <- tolower(Sys.getenv("IMPACT_SC_REMOVE_DOUBLETS", "false")) == "true"
regress_cell_cycle <- tolower(Sys.getenv("IMPACT_SC_REGRESS_CELL_CYCLE", "false")) == "true"
qc_min_nfeature_rna <- as.numeric(Sys.getenv("IMPACT_SC_QC_MIN_NFEATURE_RNA", 200))
qc_max_nfeature_rna <- as.numeric(Sys.getenv("IMPACT_SC_QC_MAX_NFEATURE_RNA", 6000))
qc_max_percent_mt <- as.numeric(Sys.getenv("IMPACT_SC_QC_MAX_PERCENT_MT", 10))
# --- NEW: User-configurable PCA dimensions ---
# Read from environment variable, default to 50 if not set
pca_dims <- as.numeric(Sys.getenv("IMPACT_SC_PCA_DIMS", 50))


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

# Input: Expects either 'ori.RDS' file or directory with 10x format files
# Output: Saves 'module1_processed.RDS' in the 'output/' directory.

# Read input data path from environment variable, default to ../ori.RDS if not set.
obj_rds_path <- Sys.getenv("IMPACT_SC_INPUT_DATA_PATH", "../ori.RDS") 

# --- UNCHANGED FILE LOADING BLOCK (FROM YOUR ORIGINAL FILE) ---
obj <- tryCatch({
  # Check if the path is a directory (10x format) or a file (RDS format)
  if (dir.exists(obj_rds_path)) {
    # Check for 10x format files
    matrix_file <- file.path(obj_rds_path, "matrix.mtx")
    genes_file <- file.path(obj_rds_path, "genes.tsv")
    barcodes_file <- file.path(obj_rds_path, "barcodes.tsv")
    
    if (file.exists(matrix_file) && file.exists(genes_file) && file.exists(barcodes_file)) {
      message(paste("Loading 10x format data from:", obj_rds_path))
      
      # Read 10x data
      data_10x <- Read10X(data.dir = obj_rds_path)
      
      # Create Seurat object
      obj <- CreateSeuratObject(counts = data_10x, project = "IMPACT_sc_demo")
      message(paste("Created Seurat object with", ncol(obj), "cells and", nrow(obj), "genes"))
      
      obj
    } else {
      stop("Directory does not contain required 10x files (matrix.mtx, genes.tsv, barcodes.tsv)")
    }
  } else if (file.exists(obj_rds_path)) {
    # Load RDS file
    message(paste("Loading RDS data from:", obj_rds_path))
    readRDS(obj_rds_path)
  } else {
    stop(paste("Input data path does not exist:", obj_rds_path))
  }
}, error = function(e) {
  stop(paste("Error loading data from", obj_rds_path, ":", e$message))
  return(NULL)
})
# --- END OF UNCHANGED FILE LOADING BLOCK ---

if (is.null(obj)) stop("Failed to load input data for IMPACT-sc.")

# --- UNCHANGED ---
if (length(unique(obj$orig.ident)) > 1) {
    message("Multiple samples detected by 'orig.ident'.")
}
if (inherits(obj[["RNA"]], "list") || (is.list(obj@assays$RNA) && inherits(obj@assays$RNA[[1]], "Assay"))) {
    message("RNA assay is split. Joining layers for initial QC and processing.")
    obj <- JoinLayers(obj, assay="RNA")
}

# --- REPLACEMENT FOR ORIGINAL QC ---
# Calculate mitochondrial percentage first
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
# Generate and save QC violin plot BEFORE filtering
message("Generating pre-filtering QC violin plot...")
qc_plot_path <- file.path(base_output_path, "01_qc_violin_plot_before_filtering.png")
tryCatch({
    png(qc_plot_path, width=12, height=6, units="in", res=300)
    print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
    dev.off()
    message(paste("QC plot saved to:", qc_plot_path))
}, error=function(e) { message(paste("Could not save QC plot:", e)) })
# Apply QC filtering based on user-defined parameters
message(paste("Applying QC filters: nFeature_RNA >", qc_min_nfeature_rna, "& nFeature_RNA <", qc_max_nfeature_rna, "& percent.mt <", qc_max_percent_mt))
obj <- subset(obj, subset = nFeature_RNA > qc_min_nfeature_rna & nFeature_RNA < qc_max_nfeature_rna & percent.mt < qc_max_percent_mt)
message(paste("Cells after QC filtering:", ncol(obj)))
# --- END OF QC REPLACEMENT ---


# --- NEW: OPTIONAL DOUBLET REMOVAL ---
if (remove_doublets) {
    message("Attempting to remove doublets with scDblFinder...")
    tryCatch({
        # To avoid issues with object states, convert to SCE, find doublets, then filter the original Seurat object
        sce <- as.SingleCellExperiment(obj)
        sce <- scDblFinder(sce)
        obj$scDblFinder.class <- sce$scDblFinder.class
        obj <- subset(obj, scDblFinder.class == "singlet")
        message(paste("Cells after doublet removal:", ncol(obj)))
    }, error = function(e) {
        warning(paste("Doublet removal step failed:", e$message, ". Continuing without it."))
    })
}


# --- UNCHANGED: Standard single-cell processing ---
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)


# --- MODIFICATION FOR OPTIONAL CELL CYCLE REGRESSION ---
vars_to_regress <- NULL
if (regress_cell_cycle) {
    message("Calculating cell cycle scores...")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    if(species == "mouse"){
        tryCatch({
            s.genes <- homologene(s.genes, inTax = 9606, outTax = 10090)$`10090`
            g2m.genes <- homologene(g2m.genes, inTax = 9606, outTax = 10090)$`10090`
        }, error = function(e){
            warning("Could not convert cell cycle genes to mouse homologs. Using human genes.")
        })
    }
    obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    vars_to_regress <- c("S.Score", "G2M.Score")
    message("Cell cycle scores will be regressed out during scaling.")
}
# Pass the variables to regress to ScaleData
obj <- ScaleData(obj, vars.to.regress = vars_to_regress, verbose = FALSE)
# --- END OF MODIFICATION ---


# --- MODIFIED PCA STEP ---
message(paste("Running PCA and computing", pca_dims, "principal components."))
obj <- RunPCA(obj, npcs = pca_dims, verbose = FALSE)
# --- END MODIFIED PCA STEP ---


module1_output_path <- file.path(base_output_path, "01_module1_processed.RDS") 
tryCatch({
  saveRDS(obj, module1_output_path)
  message("Module 1 output saved to: ", module1_output_path)
}, error = function(e) {
  warning(paste("Error saving Module 1 output:", e$message))
})

message("Finished Module 1: Data Processing")
