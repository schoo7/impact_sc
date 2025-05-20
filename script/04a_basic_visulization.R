# IMPACT-sc Script: 04a_basic_visualization.R
# Purpose: Perform basic visualizations of the annotated Seurat object.

# --- Libraries ---
# (Ensure these are loaded, typically via a shared setup or at the start of each script)
library(AnnotationDbi); library(CARD); library(celldex); library(decoupleR); library(dplyr)
library(ensembldb); library(ggplot2); library(ggpubr); library(homologene); library(infercnv)
library(Matrix); library(msigdbr); library(patchwork); library(presto); library(reticulate)
library(rjags); library(scDblFinder); library(scater); library(scRNAtoolVis); library(SCpubr)
library(Seurat); library(SeuratDisk); # library(SeuratExtend) 
library(SingleCellExperiment); library(SingleR); library(SpatialExperiment); library(scran)
library(stringr); library(tibble); library(tidyr); library(UCell); library(viridis)
library(reshape2); library(pheatmap)

# --- IMPACT-sc Script Parameters ---
species <- "human" # MODIFY THIS PARAMETER AS NEEDED ("human" or "mouse")
set.seed(1234)
base_output_path <- "../output/" 
if (!dir.exists(base_output_path)) dir.create(base_output_path, recursive = TRUE)
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { library(org.Hs.eg.db); org_db_object <- org.Hs.eg.db } else { stop("org.Hs.eg.db not installed.") }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { library(org.Mm.eg.db); org_db_object <- org.Mm.eg.db } else { stop("org.Mm.eg.db not installed.") }
} else { stop("Species not supported.") }
options(future.globals.maxSize = 4 * 1024^3)
# --- End Script Parameters ---

#### Module 4.1: Basic Visualization ####
message("Starting Module 4.1: Basic Visualization")

# Input: '03_module3_final_annotated.RDS'
# Output: Various plots saved to 'output/' directory. The 'data' object is modified in place.
obj_module3_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
if (!file.exists(obj_module3_path)) {
  stop(paste("Annotated Seurat object from Module 3 not found at:", obj_module3_path, 
             "\nPlease run previous modules first."))
}
data <- tryCatch(readRDS(obj_module3_path), error=function(e){
  message(paste("Error loading 03_module3_final_annotated.RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.1.")
message("Loaded final annotated object for Module 4.1.")

# Set data$cell_type (prioritizing SingleR, then C2S, then Seurat clustering)
if ("singler_cell_type_main" %in% colnames(data@meta.data) && sum(!is.na(data$singler_cell_type_main)) > 0) {
  data$cell_type <- data$singler_cell_type_main
  message("Using 'singler_cell_type_main' as 'cell_type'.")
} else if ("c2s_direct_predicted_label" %in% colnames(data@meta.data) && sum(!is.na(data$c2s_direct_predicted_label)) > 0) {
  data$cell_type <- data$c2s_direct_predicted_label
  message("Using 'c2s_direct_predicted_label' as 'cell_type'.")
} else if ("cell_type_seurat_clustering" %in% colnames(data@meta.data)) {
  data$cell_type <- data$cell_type_seurat_clustering
  message("Using 'cell_type_seurat_clustering' as 'cell_type'.")
} else {
  warning("No suitable 'cell_type' column found. Please set one manually. Using 'seurat_clusters' if available.")
  if("seurat_clusters" %in% colnames(data@meta.data)) data$cell_type <- data$seurat_clusters else stop("No cell type column available for downstream analysis.")
}
data$cell_type <- as.factor(data$cell_type)

# Ensure layers are joined if they were split
is_rna_split_m4a <- (is.list(data@assays$RNA) && inherits(data@assays$RNA[[1]], "Assay")) || 
                   (length(data@assays$RNA@layers) > 0 && !is.null(names(data@assays$RNA@layers)))
if (is_rna_split_m4a) { 
    message("Joining RNA assay layers for Module 4.1 operations.")
    data <- JoinLayers(data, assay="RNA") 
}
DefaultAssay(data) <- "RNA"

# Determine UMAP reduction to use for plotting
umap_to_use_viz_m4 <- if ("umap_c2s" %in% Reductions(data)) "umap_c2s" else if ("umap" %in% Reductions(data)) "umap" else "pca"
message(paste("Using reduction '", umap_to_use_viz_m4, "' for DimPlots.", sep=""))

# Dimplot (UMAP) - Seurat's DimPlot
p_dimplot_seurat <- DimPlot(data, reduction = umap_to_use_viz_m4, group.by = "cell_type", label = TRUE) + NoLegend()
ggsave(file.path(base_output_path, "04a_basicviz_dimplot_seurat.png"), plot = p_dimplot_seurat, width = 8, height = 7)

# Proportion Plots - SCpubr's do_BarPlot
if ("orig.ident" %in% colnames(data@meta.data) && nlevels(as.factor(data$orig.ident)) > 0) {
  p_prop_scpubr <- SCpubr::do_BarPlot(sample = data, group.by = "cell_type", split.by = "orig.ident", legend.position = "top")
  ggsave(file.path(base_output_path, "04a_basicviz_proportion_barplot_scpubr.png"), plot = p_prop_scpubr, width = 10, height = 7)
} else {message("Skipping proportion plot by orig.ident as 'orig.ident' is missing or has only one level.")}

# Feature Plots - Seurat's FeaturePlot
genes_to_plot <- c("CD3D", "CD14", "CD79A", "MS4A1", "CD8A", "FCGR3A", "S100A8", "IL7R", "PTPRC", "EPCAM") # Example genes
genes_to_plot_present <- genes_to_plot[genes_to_plot %in% rownames(data[["RNA"]])]
if(length(genes_to_plot_present) >= 1) {
    p_featureplot_seurat <- FeaturePlot(data, features = head(genes_to_plot_present, min(4, length(genes_to_plot_present))), reduction = umap_to_use_viz_m4)
    ggsave(file.path(base_output_path, "04a_basicviz_featureplot_seurat.png"), plot = p_featureplot_seurat, width = 12, height = 4 * ceiling(length(head(genes_to_plot_present,4))/2) )
} else {message("None of the example genes for FeaturePlot are present in the data.")}

# Dot Plot - Seurat's DotPlot
grouped_features_list <- list(
  "B_cell_markers" = c("MS4A1", "CD79A", "CD19"),
  "T_cell_markers" = c("CD3D", "CD3E", "CD8A", "IL7R", "CD4"),
  "Myeloid_markers" = c("CD14", "FCGR3A", "S100A8", "CD68", "LYZ"),
  "NK_cells" = c("NKG7", "KLRD1", "NCR1"),
  "Epithelial_cells" = c("EPCAM", "KRT18", "KRT8")
)
grouped_features_present <- lapply(grouped_features_list, function(g) g[g %in% rownames(data[["RNA"]])])
grouped_features_present <- Filter(function(x) length(x) > 0, grouped_features_present) 

if (length(grouped_features_present) > 0) {
  p_dotplot_seurat <- DotPlot(data, features = grouped_features_present, group.by = "cell_type") + RotatedAxis()
  ggsave(file.path(base_output_path, "04a_basicviz_dotplot_seurat.png"), plot = p_dotplot_seurat, width = max(10, 2*length(unlist(grouped_features_present))), height = max(7, 0.5*nlevels(data$cell_type)))
} else {message("No example genes for DotPlot are present in the data or groups are empty.")}

# Heatmap of top variable genes - Seurat's DoHeatmap
if (length(VariableFeatures(data, assay = "RNA")) == 0) data <- FindVariableFeatures(data, assay = "RNA", verbose=FALSE)
if (length(VariableFeatures(data, assay = "RNA")) > 0) {
    top_hvgs <- head(VariableFeatures(data), 20)
    p_heatmap_seurat_hvgs <- DoHeatmap(data, features = top_hvgs, group.by = "cell_type")
    ggsave(file.path(base_output_path, "04a_basicviz_heatmap_seurat_hvgs.png"), plot = p_heatmap_seurat_hvgs, width = 12, height = 10)
} else {message("No variable features found for HVG heatmap.")}

# Save the Seurat object (data) after this module if subsequent modules need its state
# Or, only save at the very end of all 04x scripts. For now, let's save.
intermediate_rds_path_04a <- file.path(base_output_path, "04a_data_after_basicviz.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04a)
  message(paste("Data object after basic visualization saved to:", intermediate_rds_path_04a))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.1:", e$message))
})

message("Finished Module 4.1: Basic Visualization")
