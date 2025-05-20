# IMPACT-sc Script: 03_cell_type_annotation.R
# Purpose: Annotate cell types using Seurat clustering, C2S predictions, and SingleR.

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
options(future.globals.maxSize = 4 * 1024^3) # 4GB, adjust as needed
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

# Ensure layers are joined if they were split
is_rna_split_m3 <- (is.list(obj@assays$RNA) && inherits(obj@assays$RNA[[1]], "Assay")) || 
                   (length(obj@assays$RNA@layers) > 0 && !is.null(names(obj@assays$RNA@layers)))
if (is_rna_split_m3) {
    message("Joining RNA assay layers for Module 3 operations.")
    obj <- JoinLayers(obj, assay = "RNA")
}
DefaultAssay(obj) <- "RNA"

#### Module 3.1: Annotation with Seurat Clustering (Refined) ####
message("Module 3.1: Seurat Clustering based annotation")
reduction_for_annot_clustering <- "harmony" 
if (!("harmony" %in% Reductions(obj)) && ("pca" %in% Reductions(obj))) {
    reduction_for_annot_clustering <- "pca"
    message("Harmony reduction not found, using PCA for clustering.")
} else if (!("harmony" %in% Reductions(obj)) && !("pca" %in% Reductions(obj))) {
    stop("Neither harmony nor pca reduction found for clustering in Module 3.1.")
}
message(paste("Using reduction '", reduction_for_annot_clustering, "' for clustering.", sep=""))

dims_for_annot_clustering <- 1:min(30, ncol(obj[[reduction_for_annot_clustering]]))

obj <- FindNeighbors(obj, dims = dims_for_annot_clustering, reduction = reduction_for_annot_clustering, 
                     graph.name = "snn_annot", verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, graph.name = "snn_annot", 
                    cluster.name = "seurat_clusters_for_annot", verbose = FALSE)
obj$cell_type_seurat_clustering <- obj$seurat_clusters_for_annot
message("Seurat clustering performed, results in 'cell_type_seurat_clustering'.")

# Visualization of clusters
umap_to_plot_m3 <- if ("umap_c2s" %in% Reductions(obj)) "umap_c2s" else if ("umap" %in% Reductions(obj)) "umap" else reduction_for_annot_clustering
if (!is.null(umap_to_plot_m3)) {
    p_clusters_m3 <- DimPlot(obj, reduction = umap_to_plot_m3, group.by = "cell_type_seurat_clustering", label = TRUE) + NoLegend()
    ggsave(file.path(base_output_path, "03_umap_seurat_clusters.png"), plot = p_clusters_m3, width = 8, height = 7)
    message("UMAP of Seurat clusters (Module 3.1) saved.")
}

#### Module 3.2: Annotation with Cell2Sentence Prediction (if available) ####
message("Module 3.2: Adding Cell2Sentence direct predictions (if available from Module 2c)")
# This step assumes 'c2s_direct_predicted_label' might have been added in script 02c.
if ("c2s_direct_predicted_label" %in% colnames(obj@meta.data)) {
    message("Found 'c2s_direct_predicted_label' in metadata from Cell2Sentence.")
} else {
    message("'c2s_direct_predicted_label' not found. This column is added if C2S predictions were successfully loaded in script 02c.")
}

#### Module 3.3: Annotation with SingleR Prediction ####
message("Module 3.3: Annotation with SingleR")
ref_data_celldex <- NULL
celldex_ref_name <- ""

if (species == "human") {
  celldex_ref_name <- "HumanPrimaryCellAtlasData"
  ref_data_celldex <- tryCatch({ celldex::HumanPrimaryCellAtlasData(ensembl = FALSE) }, 
    error = function(e) { message(paste("Failed to download", celldex_ref_name, "with symbols:", e$message, "\nTrying with ENSEMBL IDs."));
    hpca_ensembl <- tryCatch(celldex::HumanPrimaryCellAtlasData(ensembl = TRUE), error=function(e2){message(e2);NULL});
    if(is.null(hpca_ensembl)) return(NULL); ensembl_ids <- rownames(hpca_ensembl);
    gene_symbols <- mapIds(org_db_object, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first");
    valid_map_indices <- !is.na(gene_symbols); hpca_filtered <- hpca_ensembl[valid_map_indices, ];
    unique_gene_symbols <- make.unique(gene_symbols[valid_map_indices]); rownames(hpca_filtered) <- unique_gene_symbols;
    return(hpca_filtered[!duplicated(rownames(hpca_filtered)), ])
  })
} else if (species == "mouse") {
  celldex_ref_name <- "MouseRNAseqData"
  ref_data_celldex <- tryCatch({ celldex::MouseRNAseqData(ensembl = FALSE) }, 
    error = function(e) { message(paste("Failed to download", celldex_ref_name, "with symbols:", e$message, "\nTrying with ENSEMBL IDs."));
    mrna_ensembl <- tryCatch(celldex::MouseRNAseqData(ensembl = TRUE), error=function(e2){message(e2);NULL});
    if(is.null(mrna_ensembl)) return(NULL); ensembl_ids <- rownames(mrna_ensembl);
    gene_symbols <- mapIds(org_db_object, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first");
    valid_map_indices <- !is.na(gene_symbols); mrna_filtered <- mrna_ensembl[valid_map_indices, ];
    unique_gene_symbols <- make.unique(gene_symbols[valid_map_indices]); rownames(mrna_filtered) <- unique_gene_symbols;
    return(mrna_filtered[!duplicated(rownames(mrna_filtered)), ])
  })
} else { stop("SingleR reference data not configured for the specified species.") }

if(is.null(ref_data_celldex) || nrow(ref_data_celldex) == 0) {
    stop(paste("Failed to load or process", celldex_ref_name, "for SingleR."))
}
message(paste("Using", celldex_ref_name, "for SingleR."))
if(nrow(ref_data_celldex) > 0) print(table(ref_data_celldex$label.main))


# Ensure obj has logcounts for SingleR
if (!("logcounts" %in% SummarizedExperiment::assayNames(as.SingleCellExperiment(obj, assay="RNA")))) {
    if ("data" %in% slotNames(obj@assays$RNA)) {
        message("Logcounts not found, attempting to use 'data' slot from RNA assay for SingleR.")
    } else {
        message("Normalizing data to create logcounts for SingleR as they are missing.")
        obj <- NormalizeData(obj, verbose = FALSE) # This creates obj@assays$RNA$data
    }
}

obj_sce <- as.SingleCellExperiment(obj, assay = "RNA") 
predictions_singler <- tryCatch({
  SingleR(test = obj_sce, assay.type.test = "logcounts", ref = ref_data_celldex, labels = ref_data_celldex$label.main)
}, error = function(e) { message(paste("SingleR annotation failed:", e$message)); NULL })

if (!is.null(predictions_singler)) {
  message("SingleR predicted cell type distribution:")
  print(table(predictions_singler$labels))
  obj$singler_cell_type_main <- predictions_singler$labels
} else {
  obj$singler_cell_type_main <- NA 
  warning("SingleR annotation was not successful.")
}

# Save the fully annotated Seurat object
module3_final_annotated_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
tryCatch({
  saveRDS(obj, module3_final_annotated_path)
  message("Module 3 (All annotations) output saved to: ", module3_final_annotated_path)
}, error = function(e) {
  warning(paste("Error saving Module 3 output:", e$message))
})

message("Finished Module 3: Cell Type Annotation")
