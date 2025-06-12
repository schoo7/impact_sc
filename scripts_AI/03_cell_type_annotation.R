# IMPACT-sc Script: 03_cell_type_annotation.R
# Purpose: Annotate cell types using Seurat clustering, C2S predictions, and SingleR.
# This script now expects a local reference file to be provided via environment variable.
# The `download_data.sh` script is responsible for downloading any default reference files.

# --- Libraries ---
library(AnnotationDbi)
library(CARD)
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

# Get local SingleR reference path from environment variable. This is now the ONLY way to specify the reference.
local_singler_ref_rds_path <- Sys.getenv("IMPACT_SC_LOCAL_SINGLER_REF_PATH", "") 
# Get the name of the column containing cell type labels in the reference
ref_label_col_from_env <- Sys.getenv("IMPACT_SC_SINGLER_REF_LABEL_COL", "label.main")
# Get user's choice for the final cell_type column source
final_cell_type_source_choice <- tolower(Sys.getenv("IMPACT_SC_FINAL_CELL_TYPE_SOURCE", "auto"))

# Load species-specific annotation database (org.db might still be useful for other gene mapping if needed)
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { library(org.Hs.eg.db); org_db_object <- org.Hs.eg.db }
  else { stop("org.Hs.eg.db not installed.") }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { library(org.Mm.eg.db); org_db_object <- org.Mm.eg.db }
  else { stop("org.Mm.eg.db not installed.") }
} else { 
  message("Warning: Species not 'human' or 'mouse'. org.db object not loaded.")
  org_db_object <- NULL 
}
options(future.globals.maxSize = 4 * 1024^3) # Increased memory limit for certain operations
# --- End Script Parameters ---

#### Module 3: Cell Type Annotation ####
message("Starting Module 3: Cell Type Annotation")

# Input: '02_module2_c2s_processed.RDS' (preferred) or '02_module2_harmony.RDS'.
# Output: Saves '03_module3_final_annotated.RDS'.

obj_module2_path <- file.path(base_output_path, "02_module2_c2s_processed.RDS")
if (!file.exists(obj_module2_path)) {
    obj_module2_path <- file.path(base_output_path, "02_module2_harmony.RDS")
    message("C2S processed object not found, loading Harmony processed object for Module 3.")
}
if (file.exists(obj_module2_path)) {
  obj <- tryCatch(readRDS(obj_module2_path), error=function(e){
    message(paste("Error loading Module 2 output from", obj_module2_path, ":", e$message)); NULL
  })
  if(is.null(obj)) stop("Failed to load data for Module 3.")
} else {
  stop(paste("Processed Seurat object from Module 2 not found at:", obj_module2_path))
}

# This step ensures the Seurat object's RNA assay is compatible with various functions.
message("Processing Seurat object assay for compatibility...")
obj[["RNA3"]] <- as(object = obj[["RNA"]], Class = "Assay")
DefaultAssay(obj) <- "RNA3"
obj[["RNA"]] <- NULL
obj <- RenameAssays(object = obj, RNA3 = 'RNA')
message("Assay processing complete.")


# Ensure layers are joined if they were split by sample in previous steps
if (inherits(obj@assays$RNA, "Assay5")) { # Check if it's an Assay5 object
    rna_layers <- Layers(obj, assay = "RNA")
    # Check if layers are split (e.g., "counts.sample1", "data.sample1") or if standard layers are missing
    has_split_layers <- any(grepl("\\.", rna_layers)) 
    has_standard_layers <- "counts" %in% rna_layers && "data" %in% rna_layers
    
    if (length(rna_layers) > 0 && (has_split_layers || !has_standard_layers)) {
         message("RNA assay layers seem to be split or non-standard. Attempting JoinLayers for Module 3 operations.")
         obj <- JoinLayers(obj, assay = "RNA") # Join all layers in RNA assay
         message("RNA assay layers after JoinLayers:")
         print(Layers(obj, assay="RNA"))
    }
}
DefaultAssay(obj) <- "RNA" # Ensure RNA assay is default

#### Module 3.1: Annotation with Seurat Clustering (Refined) ####
message("Module 3.1: Seurat Clustering based annotation")
reduction_for_annot_clustering <- "harmony"
if (!("harmony" %in% Reductions(obj))) {
    if ("pca" %in% Reductions(obj)) {
        reduction_for_annot_clustering <- "pca"
        message("Harmony reduction not found, using PCA for clustering.")
    } else {
        stop("Neither harmony nor pca reduction found for clustering in Module 3.1.")
    }
}
message(paste("Using reduction '", reduction_for_annot_clustering, "' for clustering.", sep=""))

dims_available <- ncol(obj[[reduction_for_annot_clustering]])
dims_for_annot_clustering <- 1:min(30, dims_available) # Use up to 30 dims or available
message(paste("Using dimensions", min(dims_for_annot_clustering), "to", max(dims_for_annot_clustering), "for FindNeighbors."))


obj <- FindNeighbors(obj, dims = dims_for_annot_clustering, reduction = reduction_for_annot_clustering,
                     graph.name = "snn_annot", verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, graph.name = "snn_annot",
                    cluster.name = "seurat_clusters_for_annot", verbose = FALSE)
obj$cell_type_seurat_clustering <- obj$seurat_clusters_for_annot # Store in a dedicated column
message("Seurat clustering performed, results in 'cell_type_seurat_clustering'.")

# Determine UMAP reduction to use for plotting
umap_to_plot_m3 <- if ("umap_c2s" %in% Reductions(obj)) "umap_c2s" else if ("umap" %in% Reductions(obj)) "umap" else reduction_for_annot_clustering
if (umap_to_plot_m3 %in% Reductions(obj)) { # Check if chosen reduction exists
    p_clusters_m3 <- DimPlot(obj, reduction = umap_to_plot_m3, group.by = "cell_type_seurat_clustering", label = TRUE) + NoLegend()
    ggsave(file.path(base_output_path, "03_umap_seurat_clusters.png"), plot = p_clusters_m3, width = 8, height = 7)
    message(paste("UMAP of Seurat clusters (Module 3.1) using reduction '", umap_to_plot_m3, "' saved.", sep=""))
} else {
    message(paste("Chosen UMAP reduction '", umap_to_plot_m3, "' not found. Skipping UMAP plot for Seurat clusters.", sep=""))
}

#### Module 3.2: Annotation with Cell2Sentence Prediction (if available) ####
message("Module 3.2: Adding Cell2Sentence direct predictions (if available from Module 2c)")
if ("c2s_direct_predicted_label" %in% colnames(obj@meta.data)) {
    message("Found 'c2s_direct_predicted_label' in metadata from Cell2Sentence.")
    # Values are already in obj$c2s_direct_predicted_label
} else {
    message("'c2s_direct_predicted_label' not found. This column is added if C2S predictions were successfully loaded in script 02c.")
    # Ensure column exists with NAs if not found, for consistent downstream processing
    obj$c2s_direct_predicted_label <- NA 
}

#### Module 3.3: Annotation with SingleR Prediction ####
message("Module 3.3: Annotation with SingleR")
ref.data <- NULL
ref_label_col <- NULL

message(paste("Attempting to load local SingleR reference from path:", local_singler_ref_rds_path))
if (!nzchar(local_singler_ref_rds_path) || !file.exists(local_singler_ref_rds_path)) {
    stop(paste("Local SingleR reference file not found or path is empty. Path provided via IMPACT_SC_LOCAL_SINGLER_REF_PATH was:", local_singler_ref_rds_path))
}

# Load the reference data (now always from a local path)
ref.data.seurat <- tryCatch(readRDS(local_singler_ref_rds_path), error = function(e) {
    stop(paste("Error loading local SingleR reference RDS from", local_singler_ref_rds_path, ":", e$message))
})

# Special processing for the default Azimuth demo file
if (basename(local_singler_ref_rds_path) == "bmcite_demo.rds") {
    message("Default reference 'bmcite_demo.rds' detected. Subsetting to 50 cells per cell type for demonstration.")
    meta <- ref.data.seurat@meta.data
    barcodes <- meta %>%
      mutate(cell = rownames(.)) %>%
      group_by(celltype.l1) %>%
      slice_sample(n = 50) %>%
      pull(cell)
    ref.data.seurat <- subset(ref.data.seurat, cells = barcodes)
    message("Default reference subsetting complete.")
}

# Convert to SingleCellExperiment and validate
ref.data <- as.SingleCellExperiment(ref.data.seurat)
ref_label_col <- ref_label_col_from_env
if (!ref_label_col %in% colnames(colData(ref.data))) {
    stop(paste0("The specified label column '", ref_label_col, "' was not found in the metadata of the local reference file."))
}
message("Successfully loaded and validated local SingleR reference.")
message(paste("Using label column:", ref_label_col))
print(table(ref.data[[ref_label_col]]))

# Proceed with SingleR annotation
# Convert query Seurat object to SingleCellExperiment for SingleR
obj_sce <- as.SingleCellExperiment(obj, assay = "RNA")

# Ensure logcounts are present, normalizing if necessary
if (!("logcounts" %in% SummarizedExperiment::assayNames(obj_sce))) {
    message("Logcounts slot not found. Normalizing data to create logcounts for SingleR.")
    obj_sce <- logNormCounts(obj_sce)
}

message("Running SingleR prediction...")
predictions <- tryCatch({
    SingleR(test = obj_sce, assay.type.test = 1, # Use the first assay (logcounts)
            ref = ref.data, labels = ref.data[[ref_label_col]])
}, error = function(e) { 
    message(paste("SingleR annotation step failed:", e$message))
    NULL
})

if (!is.null(predictions) && "labels" %in% names(predictions)) {
    message("SingleR predicted cell type distribution:")
    print(table(predictions$labels))
    obj$singler_cell_type_main <- predictions$labels
} else {
    warning("SingleR annotation was not successful. 'singler_cell_type_main' will be NA.")
    obj$singler_cell_type_main <- NA
}

#### Module 3.4: Final Cell Type Assignment based on User Choice or Priority ####
message("Module 3.4: Assigning final 'cell_type' column based on user choice or priority.")
obj$cell_type <- NA # Initialize the final cell_type column

# Check for existence of annotation columns and assign NAs if not present
if (!"cell_type_seurat_clustering" %in% colnames(obj@meta.data)) {
    obj$cell_type_seurat_clustering <- NA
    message("Column 'cell_type_seurat_clustering' not found, initialized with NAs.")
}
if (!"c2s_direct_predicted_label" %in% colnames(obj@meta.data)) {
    obj$c2s_direct_predicted_label <- NA
    message("Column 'c2s_direct_predicted_label' not found, initialized with NAs.")
}
if (!"singler_cell_type_main" %in% colnames(obj@meta.data)) {
    obj$singler_cell_type_main <- NA
    message("Column 'singler_cell_type_main' not found, initialized with NAs.")
}

# Assign final cell_type based on choice
if (final_cell_type_source_choice == "singler") {
    if (sum(!is.na(obj$singler_cell_type_main)) > 0) {
        obj$cell_type <- obj$singler_cell_type_main
        message("Using 'singler_cell_type_main' as the final 'cell_type'.")
    } else {
        warning("'singler' was chosen for final_cell_type_source, but 'singler_cell_type_main' has no valid annotations. Falling back to auto-priority.")
        final_cell_type_source_choice <- "auto" # Force auto fallback
    }
} else if (final_cell_type_source_choice == "c2s") {
    if (sum(!is.na(obj$c2s_direct_predicted_label)) > 0) {
        obj$cell_type <- obj$c2s_direct_predicted_label
        message("Using 'c2s_direct_predicted_label' as the final 'cell_type'.")
    } else {
        warning("'c2s' was chosen for final_cell_type_source, but 'c2s_direct_predicted_label' has no valid annotations. Falling back to auto-priority.")
        final_cell_type_source_choice <- "auto" # Force auto fallback
    }
} else if (final_cell_type_source_choice == "seurat") {
    if (sum(!is.na(obj$cell_type_seurat_clustering)) > 0) {
        obj$cell_type <- obj$cell_type_seurat_clustering
        message("Using 'cell_type_seurat_clustering' as the final 'cell_type'.")
    } else {
        warning("'seurat' was chosen for final_cell_type_source, but 'cell_type_seurat_clustering' has no valid annotations. Falling back to auto-priority.")
        final_cell_type_source_choice <- "auto" # Force auto fallback
    }
}

# Auto-priority (applies if choice is "auto" or if chosen method had no valid annotations)
if (final_cell_type_source_choice == "auto") {
    message("Using automatic priority for final 'cell_type': SingleR > C2S > Seurat Clustering.")
    if (sum(!is.na(obj$singler_cell_type_main)) > 0) {
        obj$cell_type <- obj$singler_cell_type_main
        message("Auto-selected 'singler_cell_type_main' as 'cell_type'.")
    } else if (sum(!is.na(obj$c2s_direct_predicted_label)) > 0) {
        obj$cell_type <- obj$c2s_direct_predicted_label
        message("Auto-selected 'c2s_direct_predicted_label' as 'cell_type' (SingleR was NA).")
    } else if (sum(!is.na(obj$cell_type_seurat_clustering)) > 0) {
        obj$cell_type <- obj$cell_type_seurat_clustering
        message("Auto-selected 'cell_type_seurat_clustering' as 'cell_type' (SingleR and C2S were NA).")
    } else {
        warning("All potential annotation sources ('singler_cell_type_main', 'c2s_direct_predicted_label', 'cell_type_seurat_clustering') are NA. 'cell_type' column will be all NA.")
    }
}

# Final check and summary for the 'cell_type' column
if (sum(is.na(obj$cell_type)) == ncol(obj)) {
    warning("The final 'cell_type' column is all NA. Please check annotation steps and input data.")
} else {
    message("Final 'cell_type' column has been assigned. Distribution:")
    print(table(obj$cell_type, useNA = "ifany"))
}


# Save the final annotated Seurat object
module3_final_annotated_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
tryCatch({
  saveRDS(obj, module3_final_annotated_path)
  message("Module 3 (All annotations and final cell_type choice) output saved to: ", module3_final_annotated_path)
}, error = function(e) {
  warning(paste("Error saving Module 3 output:", e$message))
})

message("Finished Module 3: Cell Type Annotation")
