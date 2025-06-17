# IMPACT-sc Script: 02a_harmony_and_c2s_prep.R
# Purpose: Harmony batch correction and preparation for Cell2Sentence.

# --- Libraries ---
library(AnnotationDbi)
library(CARD)
library(celldex)
library(decoupleR)
library(dplyr)
library(ensembldb)
library(ggplot2)
library(ggpubr)
library(homologene)
library(harmony)
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
species <- Sys.getenv("IMPACT_SC_SPECIES", "human")
set.seed(1234)
base_output_path <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", "../output/")
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}
message(paste("Output directory set to:", base_output_path))

# --- NEW: User-configurable clustering parameters ---
# Read from environment variables, with defaults if not set
cluster_resolution <- as.numeric(Sys.getenv("IMPACT_SC_CLUSTER_RESOLUTION", 0.1))
dims_for_clustering_user <- as.numeric(Sys.getenv("IMPACT_SC_DIMS_FOR_CLUSTERING", 50))
# --- End new section ---

# Load species-specific annotation database
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

# Input path: Use base_output_path
obj_module1_path <- file.path(base_output_path, "01_module1_processed.RDS")
if (file.exists(obj_module1_path)) {
  obj <- tryCatch(readRDS(obj_module1_path), error = function(e) {
    message(paste("Error loading 01_module1_processed.RDS:",e$message)); NULL
  })
  if(is.null(obj)) stop("Failed to load data for Module 2.1")

  # Update Seurat object to ensure compatibility with Seurat v5 (especially Assay5 structure)
  message("Updating Seurat object to the latest structure...")
  obj <- UpdateSeuratObject(obj)
  message("Seurat object update complete.")

} else {
  stop("Module 1 output (01_module1_processed.RDS) not found at: ", obj_module1_path)
}

# Set default assay to RNA if it exists and is not already the default
if ("RNA" %in% Assays(obj) && DefaultAssay(obj) != "RNA") {
    DefaultAssay(obj) <- "RNA"
    message(paste("Default assay set to:", DefaultAssay(obj)))
}


if (length(unique(obj$orig.ident)) > 1) {
    message("Multiple samples found. Proceeding with Harmony integration.")

    # Ensure the PCA reduction exists, as specified by orig.reduction = "pca"
    if (!("pca" %in% Reductions(obj))) {
        stop("PCA reduction not found in the object. Run PCA before Harmony integration.")
    }
    if (ncol(obj[["pca"]]) == 0) {
        stop("PCA reduction has 0 dimensions. Check PCA computation.")
    }

    # --- Start Debug and Fix for orig.ident ---
    message("Debugging 'orig.ident' before HarmonyIntegration:")
    if (is.null(obj[["orig.ident", drop = TRUE]])) { # More robust check for metadata column
        stop("'orig.ident' metadata field is NULL or does not exist. This is required for group.by.vars.")
    }
    
    # Ensure 'orig.ident' is a factor
    if (!is.factor(obj$orig.ident)) {
        message("'orig.ident' is not a factor. Converting to factor.")
        obj$orig.ident <- as.factor(obj$orig.ident)
    }
    
    message("Levels of 'orig.ident':")
    print(levels(obj$orig.ident))
    message("Table of 'orig.ident':")
    print(table(obj$orig.ident, useNA = "ifany")) # Show NAs if they exist

    if (any(is.na(obj$orig.ident))) {
        warning("There are NA values in 'orig.ident'. This might cause issues. Consider removing or imputing them if Harmony fails.")
    }
    if (length(levels(obj$orig.ident)) < 1) {
        stop("'orig.ident' has no levels after factor conversion, or all values are NA. Check data.")
    }
    if (length(unique(na.omit(obj$orig.ident))) < 2) { # Harmony needs at least 2 groups to integrate
        warning("Harmony integration requires at least two distinct groups in 'orig.ident' (excluding NAs). Skipping Harmony if not met, or it might error.")
        # Potentially skip Harmony or handle as single batch if this condition isn't met
    }
    # --- End Debug and Fix ---

    message(paste("Using assay '", DefaultAssay(obj), "' for Harmony integration with 'group.by.vars = orig.ident'.", sep=""))

    obj <- IntegrateLayers(
      object = obj,
      method = HarmonyIntegration,
      assay = DefaultAssay(obj), # Explicitly use the default assay
      orig.reduction = "pca",
      new.reduction = "harmony",
      group.by.vars = "orig.ident", # IntegrateLayers will use this to find cells from different batches
      verbose = TRUE # Set to TRUE for more detailed output from Harmony
    )
    message("Harmony integration complete.")
} else {
    message("Only one unique 'orig.ident' found (or fewer than 2 after NA removal). Harmony integration skipped or run on single batch.")
    if (!("harmony" %in% Reductions(obj))) {
      # Ensure PCA exists before trying to copy it
      if ("pca" %in% Reductions(obj)) {
        obj[["harmony"]] <- obj[["pca"]]
        message("'harmony' reduction set to be same as 'pca'.")
      } else {
        warning("PCA reduction not found. Cannot set 'harmony' to 'pca' for single sample.")
      }
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

# Ensure DefaultAssay is "RNA" before clustering and downstream processing
if ("RNA" %in% Assays(obj) && DefaultAssay(obj) != "RNA") {
    DefaultAssay(obj) <- "RNA"
    message(paste("Default assay set to:", DefaultAssay(obj), "for C2S preparation."))
} else if (!("RNA" %in% Assays(obj))) {
    stop("RNA assay not found in the object for C2S preparation.")
}


if (inherits(obj@assays[[DefaultAssay(obj)]], "Assay5")) {
    assay_obj_check <- obj@assays[[DefaultAssay(obj)]]
    if (length(Layers(obj, assay = DefaultAssay(obj))) > 0 && !"data" %in% Layers(obj, assay = DefaultAssay(obj)) && !"counts" %in% Layers(obj, assay = DefaultAssay(obj))) {
        message("The default assay has multiple layers but no single 'data' or 'counts' layer. Consider JoinLayers().")
    }
}


# Check if the reduction exists and has enough dimensions
if (!(reduction_for_clustering %in% Reductions(obj))) {
    stop(paste("Reduction '", reduction_for_clustering, "' not found in the object. Check previous steps.", sep=""))
}
if (ncol(obj[[reduction_for_clustering]]) == 0) {
    stop(paste("Reduction '", reduction_for_clustering, "' has 0 dimensions. Check previous steps.", sep=""))
}

# --- MODIFIED CLUSTERING PARAMETERS ---
# Use the user-defined number of dimensions, ensuring it doesn't exceed available dimensions
dims_for_clustering <- 1:min(dims_for_clustering_user, ncol(obj[[reduction_for_clustering]]))
message(paste("Using dimensions", min(dims_for_clustering), "to", max(dims_for_clustering), "from reduction '", reduction_for_clustering, "' for FindNeighbors."))
message(paste("User-selected clustering resolution is", cluster_resolution, "."))
message(paste("Note for C2S users: recommended dims_for_clustering is 1024. Normal usage: 50."))
# --- END MODIFIED CLUSTERING PARAMETERS ---


obj <- FindNeighbors(obj, dims = dims_for_clustering, reduction = reduction_for_clustering, verbose = FALSE)

# --- MODIFIED CLUSTERING ---
# Dynamically create cluster name based on user-selected resolution
low_res_cluster_name <- paste0("rna_snn_res.", cluster_resolution)
message(paste("Finding low-resolution clusters with resolution =", cluster_resolution, "and storing in '", low_res_cluster_name, "'"))

# Run clustering with user-defined resolution
obj <- FindClusters(obj, resolution = cluster_resolution, cluster.name = low_res_cluster_name, verbose = FALSE)

# Also run a high-resolution clustering for other potential uses (as in original script)
message("Finding high-resolution clusters with resolution = 2 for additional analysis.")
obj <- FindClusters(obj, resolution = 2, cluster.name = "rna_snn_res.2", verbose = FALSE)

# Set cell_type_low_res from the user-defined resolution cluster results
obj$cell_type_low_res <- obj[[low_res_cluster_name, drop=TRUE]]
message(paste("'cell_type_low_res' column created from the clustering results of resolution", cluster_resolution))
# --- END MODIFIED CLUSTERING ---


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
  # Ensure the correct assay and its data are used for saving.
  current_default_assay <- DefaultAssay(obj)
  message(paste("Preparing assay '", current_default_assay, "' for H5AD conversion.", sep=""))
  
  if (length(Layers(obj, assay = current_default_assay)) > 0 && 
      any(grepl("\\.", Layers(obj, assay = current_default_assay)))) { 
      message(paste("Joining layers for assay '", current_default_assay, "' before H5AD conversion.", sep=""))
      obj <- JoinLayers(obj, assay = current_default_assay)
      message(paste("Layers for assay '", current_default_assay, "' after JoinLayers:", sep=""))
      print(Layers(obj, assay = current_default_assay))
  }
  
  message(paste("Saving H5Seurat using assay:", current_default_assay))
  message("Available layers in this assay for H5Seurat:")
  print(Layers(obj, assay = current_default_assay))

  # Convert to a standard Assay from Assay5 if needed for compatibility
  obj[["RNA3"]] <- as(obj = obj[["RNA"]], Class = "Assay")  
  DefaultAssay(obj) <- "RNA3"
  obj[["RNA"]] <- NULL
  obj <- RenameAssays(obj = obj, RNA3 = 'RNA') 

  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE, assay = 'RNA')
  Convert(h5seurat_path, dest = "h5ad", assay = 'RNA', overwrite = TRUE)
  message(paste("H5Seurat saved to:", h5seurat_path))
  message(paste("H5AD (for Python Cell2Sentence) saved to:", h5ad_path))
}, error = function(e) {
  stop(paste("Error during H5Seurat/H5AD conversion:", e$message))
})

# Set environment variables for the Python script
Sys.setenv(H5AD_FILE_PATH = h5ad_path)
c2s_output_dir_full_path <- file.path(base_output_path, "C2S_output")
Sys.setenv(C2S_OUTPUT_DIR = c2s_output_dir_full_path)
Sys.setenv(C2S_EMBEDDINGS_CSV = file.path(c2s_output_dir_full_path, "c2s_cell_embeddings.csv"))
Sys.setenv(C2S_PREDICTED_CSV = file.path(c2s_output_dir_full_path, "predicted_cell_types.csv"))
Sys.setenv(SPECIES_FOR_C2S = species)

if (!dir.exists(c2s_output_dir_full_path)) {
  tryCatch(dir.create(c2s_output_dir_full_path, recursive = TRUE),
           error = function(e) warning(paste("Error creating C2S output directory:", e$message)))
}

message("Finished Module 2a: Harmony and Cell2Sentence Preparation.")
message(paste("Next step: Run '02b_run_cell2sentence.py' with H5AD_FILE_PATH set to:", h5ad_path))
