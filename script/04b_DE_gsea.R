# IMPACT-sc Script: 04b_diffex_gsea.R
# Purpose: Perform Differential Expression and Gene Set Enrichment Analysis.

# --- Libraries ---
# (Standard library loading block)
library(AnnotationDbi); library(CARD); library(celldex); library(decoupleR); library(dplyr)
library(ensembldb); library(ggplot2); library(ggpubr); library(homologene); library(infercnv)
library(Matrix); library(msigdbr); library(patchwork); library(presto); library(reticulate)
library(rjags); library(scDblFinder); library(scater); library(scRNAtoolVis); library(SCpubr)
library(Seurat); library(SeuratDisk); # library(SeuratExtend) 
library(SingleCellExperiment); library(SingleR); library(SpatialExperiment); library(scran)
library(stringr); library(tibble); library(tidyr); library(UCell); library(viridis)
library(reshape2); library(pheatmap)

# --- IMPACT-sc Script Parameters ---
species <- "human" 
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

#### Module 4.2: Differential Expression & GSEA ####
message("Starting Module 4.2: Differential Expression & GSEA")

# Input: '04a_data_after_basicviz.RDS' (or '03_module3_final_annotated.RDS' if 04a is skipped)
# Output: DGE results, plots. Modifies 'data' object.
obj_prev_module_path <- file.path(base_output_path, "04a_data_after_basicviz.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
  message(paste("04a output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required input Seurat object not found at:", obj_prev_module_path))
}
data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading input RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.2.")
message("Loaded Seurat object for Module 4.2.")
DefaultAssay(data) <- "RNA" # Ensure RNA is default for DGE

# Gene expression across cell types (Violin/Box plots)
example_gene_dge <- "CD3D" 
if (example_gene_dge %in% rownames(data[["RNA"]])) {
  p_vln_seurat <- VlnPlot(data, features = example_gene_dge, group.by = "cell_type", pt.size = 0.1) 
  if(nlevels(data$cell_type) > 1) { 
      p_vln_seurat <- p_vln_seurat + stat_compare_means(aes(group = cell_type), label = "p.signif", method="wilcox.test", label.y.npc = 0.95) 
  }
  ggsave(file.path(base_output_path, "04b_dge_vlnplot_seurat.png"), plot = p_vln_seurat, width = max(7, 2*nlevels(data$cell_type)), height = 6)
}

# Identify marker genes
if (!("cell_type" %in% colnames(data@meta.data))) stop ("'cell_type' column not found in metadata for DGE.")
Idents(data) <- "cell_type"
all_markers <- tryCatch(FindAllMarkers(
  data, min.cells.feature = 10, min.cells.group = 10,
  min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE, verbose = FALSE
), error = function(e){ message(paste("FindAllMarkers failed:", e$message)); NULL})

if(!is.null(all_markers) && nrow(all_markers) > 0){
  write.csv(all_markers, file.path(base_output_path, "04b_all_celltype_markers.csv"), row.names = FALSE)
  message("Differential expression markers saved.")
  top10_markers_m4.2 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  print(head(top10_markers_m4.2))

  # Volcano plot using scRNAtoolVis::jjVolcano
  if (nlevels(data$cell_type) > 0) {
    first_cluster_name <- levels(data$cell_type)[1]
    first_cluster_markers <- dplyr::filter(all_markers, cluster == first_cluster_name)
    if (nrow(first_cluster_markers) > 0 && all(c("avg_log2FC", "p_val_adj") %in% colnames(first_cluster_markers))) {
      volcano_plot_m4.2 <- jjVolcano(
        diffData = first_cluster_markers, log2FC.cutoff = 0.5, FDR.cutoff = 0.05, topGeneN = 10, fontface = 'italic'
      )
      ggsave(file.path(base_output_path, "04b_volcano_plot_example.png"), plot = volcano_plot_m4.2, width = 8, height = 7)
      message("Example volcano plot saved.")
    } else { message("Not enough data or required columns for jjVolcano in Module 4.2.") }
  }
} else {message("No markers found or FindAllMarkers failed.")}

# GSEA with SeuratExtend (GeneSetAnalysisGO, RenameGO, GSEAplot) - Commented out due to dependency
message("GSEA with SeuratExtend is commented out due to specific dependencies. Standard GSEA tools like clusterProfiler or fgsea can be used as an alternative. These tools typically require a ranked list of genes from DGE results.")
# Example placeholder for where clusterProfiler might be used:
# if(!is.null(all_markers) && nrow(all_markers) > 0 && "avg_log2FC" %in% colnames(all_markers) && "gene" %in% colnames(all_markers)){
#   # ranked_genes <- all_markers %>% arrange(-avg_log2FC) %>% dplyr::select(gene, avg_log2FC) %>% tibble::deframe()
#   # if(requireNamespace("clusterProfiler", quietly=TRUE) && requireNamespace("enrichplot", quietly=TRUE)){
#   #   # gse_result <- clusterProfiler::gseGO(geneList=ranked_genes, OrgDb=org_db_object, keyType="SYMBOL", ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05)
#   #   # if(!is.null(gse_result) && nrow(gse_result) > 0) {
#   #   #   dotplot_gsea <- enrichplot::dotplot(gse_result)
#   #   #   ggsave(file.path(base_output_path, "04b_gsea_dotplot_example_clusterprofiler.png"), plot=dotplot_gsea)
#   #   # }
#   # }
# }

intermediate_rds_path_04b <- file.path(base_output_path, "04b_data_after_diffex.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04b)
  message(paste("Data object after DGE saved to:", intermediate_rds_path_04b))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.2:", e$message))
})

message("Finished Module 4.2: Differential Expression & GSEA")
