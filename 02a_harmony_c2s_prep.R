# IMPACT-sc Script: 02a_harmony_and_c2s_prep.R
# Purpose: Harmony batch correction and preparation for Cell2Sentence.

# --- Libraries ---
# (Same library loading block as in 01_data_processing.R)
library(AnnotationDbi)
library(CARD) 
library(celldex) 
library(decoupleR) 
library(dplyr)
library(ensembldb) 
library(ggplot2)
library(ggpubr)
library(homologene) 
library(infercnv) 
library(Matrix)
library(msigdbr) 
library(patchwork)
library(presto) 
library(reticulate) 
library(rjags) 
library(scDblFinder) 
library(scater)
library(scRNAtoolVis) 
library(SCpubr) 
library(Seurat)
library(SeuratDisk) 
# library(SeuratExtend) 
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

# --- IMPACT-sc Script Parameters ---
species <- "human" # MODIFY THIS PARAMETER AS NEEDED ("human" or "mouse")
set.seed(1234)
base_output_path <- "../output/" # Assuming scripts are in a 'scripts' subdirectory
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}
# Load species-specific annotation database (as in 01_data_processing.R)
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { library(org.Hs.eg.db); org_db_object <- org.Hs.eg.db } 
  else { stop("org.Hs.eg.db not installed.") }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { library(org.Mm.eg.db); org_db_object <- org.Mm.eg.db }
  else { stop("org.Mm.eg.db not installed.") }
} else { stop("Species not supported.") }
# --- End Script Parameters ---

#### Module 2.1: Harmony Batch Correction ####
message("Starting Module 2.1: Harmony Batch Correction")

# Input: Expects '01_module1_processed.RDS' from the 'output/' directory.
# Output: Saves '02_module2_harmony.RDS'.
obj_module1_path <- file.path(base_output_path, "01_module1_processed.RDS")
if (file.exists(obj_module1_path)) {
  obj <- tryCatch(readRDS(obj_module1_path), error = function(e) { 
    message(paste("Error loading 01_module1_processed.RDS:",e$message)); NULL
  })
  if(is.null(obj)) stop("Failed to load data for Module 2.1")
} else {
  stop("Module 1 output (01_module1_processed.RDS) not found at: ", obj_module1_path)
}

if (length(unique(obj$orig.ident)) > 1) {
    message("Multiple samples found. Proceeding with Harmony integration.")
    is_v5_layered <- !is.null(obj@assays$RNA$counts) && inherits(obj@assays$RNA$counts, "dgCMatrix") &&
                     length(obj@assays$RNA@layers) > 0 && inherits(obj@assays$RNA@layers[[1]], "dgCMatrix")
    if (!is_v5_layered && !is.list(obj@assays$RNA)) { # Check if already in list format (older Seurat)
         obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    }
    obj <- IntegrateLayers(
      object = obj, method = HarmonyIntegration, orig.reduction = "pca", 
      new.reduction = "harmony", group.by.vars = "orig.ident", verbose = FALSE
    )
} else {
    message("Only one 'orig.ident' found. Harmony integration skipped or run on single batch.")
    if (!("harmony" %in% Reductions(obj))) {
      obj[["harmony"]] <- obj[["pca"]]
      message("'harmony' reduction set to be same as 'pca'.")
    }
}

module2_harmony_path <- file.path(base_output_path, "02_module2_harmony.RDS")
tryCatch({
  saveRDS(obj, module2_harmony_path)
  message("Module 2.1 (Harmony) output saved to: ", module2_harmony_path)
}, error = function(e) {
  warning(paste("Error saving Module 2.1 output:", e$message))
})

#### Module 2.2: Cell2Sentence Data Preparation (R part) ####
message("Starting Module 2.2: Cell2Sentence Data Preparation")
# Input: Uses 'obj' from Harmony step (or '02_module2_harmony.RDS' if run separately).
# Output: Saves '02_module2_for_c2s.h5seurat' and '02_module2_for_c2s.h5ad'.
#         Also saves the Seurat object before H5AD conversion as '02_module2_c2s_prepped_object.RDS'.

reduction_for_clustering <- "harmony" 
message(paste("Using reduction '", reduction_for_clustering, "' for clustering before Cell2Sentence.", sep=""))

is_rna_split_c2s <- (is.list(obj@assays$RNA) && inherits(obj@assays$RNA[[1]], "Assay")) || 
                (length(obj@assays$RNA@layers) > 0 && !is.null(names(obj@assays$RNA@layers)))
if (is_rna_split_c2s) {
    message("Joining RNA assay layers before clustering for Cell2Sentence.")
    obj <- JoinLayers(obj, assay = "RNA")
}
DefaultAssay(obj) <- "RNA"

obj <- FindNeighbors(obj, dims = 1:min(50, ncol(obj[[reduction_for_clustering]])), reduction = reduction_for_clustering, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.1, cluster.name = "rna_snn_res.0.1", verbose = FALSE)
obj <- FindClusters(obj, resolution = 2, cluster.name = "rna_snn_res.2", verbose = FALSE)
obj$cell_type_low_res <- obj$rna_snn_res.0.1 # This will be used as input label for C2S

# Save the Seurat object at this stage, before H5AD conversion
prepped_obj_for_c2s_path <- file.path(base_output_path, "02_module2_c2s_prepped_object.RDS")
tryCatch({
  saveRDS(obj, prepped_obj_for_c2s_path)
  message("Seurat object prepped for C2S saved to: ", prepped_obj_for_c2s_path)
}, error = function(e) {
  warning(paste("Error saving C2S prepped Seurat object:", e$message))
})


h5seurat_path <- file.path(base_output_path, "02_module2_for_c2s.h5seurat")
h5ad_path <- file.path(base_output_path, "02_module2_for_c2s.h5ad")

tryCatch({
  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE, assay = "RNA")
  Convert(h5seurat_path, dest = "h5ad", assay = "RNA", overwrite = TRUE) 
  message(paste("H5Seurat saved to:", h5seurat_path))
  message(paste("H5AD (for Python Cell2Sentence) saved to:", h5ad_path))
}, error = function(e) {
  stop(paste("Error during H5Seurat/H5AD conversion:", e$message))
})

# Set environment variables for the Python script
Sys.setenv(H5AD_FILE_PATH = h5ad_path) # Python script will read this
c2s_output_dir_env <- file.path(base_output_path, "C2S_output") # Python script will read this
Sys.setenv(C2S_OUTPUT_DIR = c2s_output_dir_env)
Sys.setenv(C2S_EMBEDDINGS_CSV = file.path(c2s_output_dir_env, "c2s_cell_embeddings.csv"))
Sys.setenv(C2S_PREDICTED_CSV = file.path(c2s_output_dir_env, "predicted_cell_types.csv"))
Sys.setenv(SPECIES_FOR_C2S = species)

if (!dir.exists(c2s_output_dir_env)) {
  tryCatch(dir.create(c2s_output_dir_env, recursive = TRUE),
           error = function(e) warning(paste("Error creating C2S output directory:", e$message)))
}

message("Finished Module 2a: Harmony and Cell2Sentence Preparation.")
message(paste("Next step: Run '02b_run_cell2sentence.py' with H5AD_FILE_PATH set to:", h5ad_path))
