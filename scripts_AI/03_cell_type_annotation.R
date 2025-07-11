# IMPACT-sc Script: 03_cell_type_annotation.R
# Purpose: Annotate cell types using Seurat clustering, C2S predictions, SingleR, or the ceLLama method.
# This version reads the annotation method, ceLLama parameters, and Ollama model from environment variables.

# --- Libraries ---
# Ensure all necessary packages are installed before running.
# install.packages(c("ceLLama", "Seurat", "tidyverse", "httr", "thinkR", "AnnotationDbi", "CARD", "decoupleR", "ensembldb", "ggpubr", "homologene", "harmony", "msigdbr", "patchwork", "reticulate", "scDblFinder", "scater", "scRNAtoolVis", "SCpubr", "SeuratDisk", "SeuratExtend", "SingleCellExperiment", "SingleR", "scran", "UCell", "viridis", "pheatmap"))

# Suppress startup messages for cleaner output
suppressPackageStartupMessages({
    library(ceLLama)
    library(Seurat)
    library(tidyverse)
    library(httr)
    library(thinkR)
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
})


# --- IMPACT-sc Script Parameters (READ FROM ENVIRONMENT VARIABLES) ---
species <- Sys.getenv("IMPACT_SC_SPECIES", "human")
set.seed(1234)
base_output_path <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", "../output/")
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}
message(paste("Output directory set to:", base_output_path))

# --- Annotation Method Parameters (NEW) ---
annotation_method <- tolower(Sys.getenv("IMPACT_SC_ANNOTATION_METHOD", "auto"))
cellama_temperature <- as.numeric(Sys.getenv("IMPACT_SC_CELLAMA_TEMPERATURE", "0.0"))
ollama_model <- Sys.getenv("IMPACT_SC_OLLAMA_MODEL", "gemma3:12b-it-qat") # Default to your desired model

message(paste("Annotation method selected:", annotation_method))
if (annotation_method == "cellama") {
    message(paste("  - ceLLama Temperature:", cellama_temperature))
    message(paste("  - Ollama Model to be used (must be running):", ollama_model))
}

# --- SingleR Specific Parameters ---
local_singler_ref_rds_path <- Sys.getenv("IMPACT_SC_LOCAL_SINGLER_REF_PATH", "")
ref_label_col_from_env <- Sys.getenv("IMPACT_SC_SINGLER_REF_LABEL_COL", "label.main")
final_cell_type_source_choice <- tolower(Sys.getenv("IMPACT_SC_FINAL_CELL_TYPE_SOURCE", "auto"))

# Load species-specific annotation database
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
options(future.globals.maxSize = 4 * 1024^3)
# --- End Script Parameters ---


#### Module 3: Cell Type Annotation ####
message("Starting Module 3: Cell Type Annotation")

# Load the Seurat object from the previous step (harmony output)
obj_module2_path <- file.path(base_output_path, "02_module2_harmony.RDS")
if (file.exists(obj_module2_path)) {
  obj <- tryCatch(readRDS(obj_module2_path), error=function(e){
    message(paste("Error loading Module 2 output from", obj_module2_path, ":", e$message)); NULL
  })
  if(is.null(obj)) stop("Failed to load data for Module 3.")
} else {
  stop(paste("Processed Seurat object from Module 2 not found at:", obj_module2_path))
}

# Ensure the RNA assay is correctly formatted
message("Processing Seurat object assay for compatibility...")
DefaultAssay(obj) <- "RNA"
if (inherits(obj@assays$RNA, "Assay")) {
    obj <- JoinLayers(obj, assay = "RNA")
}
message("Assay processing complete.")


# --- Clustering (Common for all annotation methods) ---
message("Performing Seurat clustering...")
# Read clustering parameters from environment variables
dims_for_clustering <- 1:as.numeric(Sys.getenv("IMPACT_SC_DIMS_FOR_CLUSTERING", 50))
cluster_resolution <- as.numeric(Sys.getenv("IMPACT_SC_CLUSTER_RESOLUTION", 0.1))

obj <- FindNeighbors(obj, dims = dims_for_clustering, reduction = "harmony")
obj <- FindClusters(obj, resolution = cluster_resolution)

# Initialize the final celltype column
obj$celltype_final <- NA


# --- Run Annotation method if chosen ---
if (annotation_method == "cellama") {
    #### ceLLama Annotation ####
    message("--- Starting ceLLama Annotation ---")
    message("IMPORTANT: This script assumes the Ollama application is running and serving a model in a separate terminal (e.g., via 'ollama run gemma3:12b-it-qat').")
    
    # Find markers for the seurat clusters
    message("Finding markers for ceLLama...")
    obj.markers <- FindAllMarkers(obj, group.by = "seurat_clusters", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    obj.markers.list <- split(obj.markers, obj.markers$cluster)

    # Run ceLLama, which will communicate with the active Ollama service
    message(paste("Running ceLLama with temperature =", cellama_temperature))
    res <- ceLLama(
        obj.markers.list,
        temperature = cellama_temperature,
        seed = 101,
        get_reason = TRUE
    )

    # Extract annotations
    annotations <- map_chr(res, 1)

    # Apply annotations to the Seurat object
    Idents(obj) <- "seurat_clusters"
    names(annotations) <- levels(obj)
    obj <- RenameIdents(obj, annotations)

    # Store the annotations in a new metadata column
    obj$celltype_cellama <- Idents(obj)
    message("ceLLama annotation complete. Results stored in 'celltype_cellama' metadata.")

} else if (annotation_method == "singler") {
    #### SingleR Annotation ####
    message("--- Starting SingleR Annotation ---")
    
    message(paste("Attempting to load local SingleR reference from path:", local_singler_ref_rds_path))
    if (!nzchar(local_singler_ref_rds_path) || !file.exists(local_singler_ref_rds_path)) {
        stop(paste("Local SingleR reference file not found. Path provided was:", local_singler_ref_rds_path))
    }
    
    ref.data.seurat <- tryCatch(readRDS(local_singler_ref_rds_path), error = function(e) {
        stop(paste("Error loading local SingleR reference RDS:", e$message))
    })
    
    ref.data <- as.SingleCellExperiment(ref.data.seurat)
    if (!ref_label_col_from_env %in% colnames(colData(ref.data))) {
        stop(paste0("Specified label column '", ref_label_col_from_env, "' not found in reference metadata."))
    }

    obj_sce <- as.SingleCellExperiment(obj, assay = "RNA")
    if (!("logcounts" %in% SummarizedExperiment::assayNames(obj_sce))) {
        obj_sce <- logNormCounts(obj_sce)
    }

    message("Running SingleR prediction...")
    predictions <- tryCatch({
        SingleR(test = obj_sce, assay.type.test = 1, ref = ref.data, labels = ref.data[[ref_label_col_from_env]])
    }, error = function(e) {
        message(paste("SingleR annotation step failed:", e$message)); NULL
    })

    if (!is.null(predictions) && "labels" %in% names(predictions)) {
        obj$celltype_singler <- predictions$labels
        message("SingleR annotation complete. Results stored in 'celltype_singler' metadata.")
    } else {
        warning("SingleR annotation was not successful. 'celltype_singler' will be NA.")
        obj$celltype_singler <- NA
    }
}

# --- Set Final Cell Type Annotation (MODIFIED LOGIC) ---
message(paste("Setting final cell type column based on choice:", final_cell_type_source_choice))

if (final_cell_type_source_choice == "seurat") {
    obj$celltype_final <- obj$seurat_clusters
    message("Final cell type set to 'seurat_clusters'.")
} else if (final_cell_type_source_choice == "c2s") {
    if ("c2s_direct_predicted_label" %in% colnames(obj@meta.data)) {
        obj$celltype_final <- obj$c2s_direct_predicted_label
        message("Final cell type set to 'c2s_direct_predicted_label'.")
    } else {
        warning("'c2s' was chosen as final source, but 'c2s_direct_predicted_label' column not found. 'celltype_final' will be NA.")
    }
} else if (final_cell_type_source_choice == "singler") {
    if ("celltype_singler" %in% colnames(obj@meta.data) && !all(is.na(obj$celltype_singler))) {
        obj$celltype_final <- obj$celltype_singler
        message("Final cell type set to 'celltype_singler'.")
    } else {
        warning("'singler' was chosen as final source, but 'celltype_singler' column is missing or empty. This may be because SingleR was not run. 'celltype_final' will be NA.")
    }
} else if (final_cell_type_source_choice == "cellama") {
    if ("celltype_cellama" %in% colnames(obj@meta.data) && !all(is.na(obj$celltype_cellama))) {
        obj$celltype_final <- obj$celltype_cellama
        message("Final cell type set to 'celltype_cellama'.")
    } else {
        warning("'cellama' was chosen as final source, but 'celltype_cellama' column is missing or empty. This may be because ceLLama was not run. 'celltype_final' will be NA.")
    }
} else { # This handles "auto" or any other case
    message("Using 'auto' logic to determine final cell type.")
    # Auto-prioritization: ceLLama > SingleR > C2S > Seurat Clusters
    if ("celltype_cellama" %in% colnames(obj@meta.data) && !all(is.na(obj$celltype_cellama))) {
        obj$celltype_final <- obj$celltype_cellama
        message("Auto-selected 'celltype_cellama' as final cell type.")
    } else if ("celltype_singler" %in% colnames(obj@meta.data) && !all(is.na(obj$celltype_singler))) {
        obj$celltype_final <- obj$celltype_singler
        message("Auto-selected 'celltype_singler' as final cell type.")
    } else if ("c2s_direct_predicted_label" %in% colnames(obj@meta.data) && !all(is.na(obj$c2s_direct_predicted_label))) {
        obj$celltype_final <- obj$c2s_direct_predicted_label
        message("Auto-selected 'c2s_direct_predicted_label' as final cell type.")
    } else {
        obj$celltype_final <- obj$seurat_clusters
        message("Auto-selected 'seurat_clusters' as fallback final cell type.")
    }
}


# Final check and summary for the 'celltype_final' column
if (sum(is.na(obj$celltype_final), na.rm=TRUE) == ncol(obj)) {
    warning("The final 'celltype_final' column is all NA. Please check annotation steps and input data.")
} else {
    message("Final cell type distribution in 'celltype_final':")
    print(table(obj$celltype_final, useNA = "ifany"))
}

obj$cell_type <- obj$celltype_final
message("Unified cell type variable created: 'cell_type'")

# Save the final annotated Seurat object
module3_final_annotated_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
tryCatch({
  saveRDS(obj, module3_final_annotated_path)
  message(paste("Module 3 (Final Annotated Object) saved to:", module3_final_annotated_path))
}, error = function(e) {
  warning(paste("Error saving Module 3 output:", e$message))
})

# Generate and save a UMAP plot of the final cell types
if (!all(is.na(obj$celltype_final))) {
    umap_to_plot_m3 <- if ("umap" %in% Reductions(obj)) "umap" else "harmony"
    if (umap_to_plot_m3 %in% Reductions(obj)) {
        p_final_celltypes <- DimPlot(obj, reduction = umap_to_plot_m3, group.by = "celltype_final", label = TRUE, repel = TRUE) + NoLegend() + labs(title = "Final Cell Type Annotations")
        ggsave(file.path(base_output_path, "03_umap_final_celltypes.png"), plot = p_final_celltypes, width = 10, height = 8)
        message(paste("UMAP of final cell types saved to 03_umap_final_celltypes.png using reduction '", umap_to_plot_m3, "'.", sep=""))
    }
}


message("Finished Module 3: Cell Type Annotation")
