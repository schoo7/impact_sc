# IMPACT-sc Script: 04c_decoupler_analysis.R
# Purpose: Add Pathway and TF Activity Scores with DecoupleR.

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

#### Module 4.3: Add Pathway and TF Activity Scores with DecoupleR ####
message("Starting Module 4.3: DecoupleR Analysis")

# Input: '04b_data_after_diffex.RDS' (or '03_module3_final_annotated.RDS' if 04a/b skipped)
# Output: Modifies 'data' object by adding DecoupleR assays, saves plots.
obj_prev_module_path <- file.path(base_output_path, "04b_data_after_diffex.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS") # Fallback
  message(paste("04b output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required input Seurat object not found at:", obj_prev_module_path))
}
data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading input RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.3.")
message("Loaded Seurat object for Module 4.3.")
DefaultAssay(data) <- "RNA"

# Ensure log-normalized data is present for DecoupleR
if (!("data" %in% slotNames(data@assays$RNA))) { 
    message("Normalizing data as 'data' slot is missing for DecoupleR.")
    data <- NormalizeData(data, verbose=FALSE) 
} 
mat_norm_decoupler <- GetAssayData(data, assay = "RNA", layer = "data")
decoupler_organism <- if (species == "human") "human" else if (species == "mouse") "mouse" else "human"
message(paste("Using organism '", decoupler_organism, "' for DecoupleR.", sep=""))

# TF activity (CollecTRI)
net_collectri <- tryCatch(decoupleR::get_collectri(organism = decoupler_organism, split_complexes = FALSE), 
                          error = function(e){message(paste("Error getting CollecTRI network:",e$message)); NULL})
if(!is.null(net_collectri)){
  tf_acts_ulm <- decoupleR::run_ulm(as.matrix(mat_norm_decoupler), net = net_collectri, 
                                    .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
  tf_scores_wide <- tf_acts_ulm %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source')
  data[['tfs_ulm']] <- CreateAssayObject(counts = tf_scores_wide) # Store raw scores
  DefaultAssay(data) <- "tfs_ulm"
  data <- ScaleData(data, verbose = FALSE) # Scale for visualization
  
  n_top_tfs_decoupler <- 25 
  avg_tf_activity_decoupler <- AverageExpression(data, assays = "tfs_ulm", features = rownames(data[["tfs_ulm"]]), group.by = "cell_type", verbose = FALSE)$tfs_ulm
  tf_variances_decoupler <- apply(avg_tf_activity_decoupler, 1, var, na.rm = TRUE) 
  top_variable_tfs_decoupler <- names(sort(tf_variances_decoupler, decreasing = TRUE)[1:min(n_top_tfs_decoupler, length(tf_variances_decoupler), na.rm=TRUE)])
  avg_tf_activity_top_decoupler <- avg_tf_activity_decoupler[top_variable_tfs_decoupler, , drop=FALSE]

  if(nrow(avg_tf_activity_top_decoupler) > 0 && ncol(avg_tf_activity_top_decoupler) > 0) {
      tf_heatmap_path_decoupler <- file.path(base_output_path, "04c_tf_activity_heatmap_decoupler.png")
      png(tf_heatmap_path_decoupler, width=10, height=max(6, 0.3*nrow(avg_tf_activity_top_decoupler) + 2), units="in", res=300)
      pheatmap::pheatmap(avg_tf_activity_top_decoupler, scale = "row", 
                         main = "Mean TF Activity per Cell Type (DecoupleR Top Variable TFs)", 
                         fontsize_row = max(5, 10 - nrow(avg_tf_activity_top_decoupler) * 0.1), # Adjust font size
                         fontsize_col = 10) 
      dev.off()
      message("DecoupleR TF activity heatmap saved.")
  } else { message("Not enough data to plot TF activity heatmap.")}
} else { message("Skipping TF activity due to CollecTRI network loading failure.")}

# Pathway activity (PROGENy)
net_progeny <- tryCatch(decoupleR::get_progeny(organism = decoupler_organism, top = 500), 
                        error = function(e){message(paste("Error getting PROGENy network:",e$message)); NULL})
if(!is.null(net_progeny)){
  pathway_acts_mlm <- decoupleR::run_mlm(as.matrix(mat_norm_decoupler), net = net_progeny, 
                                         .source = 'source', .target = 'target', .mor = 'weight', minsize = 5)
  pathway_scores_wide <- pathway_acts_mlm %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source')
  data[['pathways_mlm']] <- CreateAssayObject(counts = pathway_scores_wide)
  DefaultAssay(data) <- "pathways_mlm"
  data <- ScaleData(data, verbose = FALSE)
  
  avg_pathway_activity_decoupler <- AverageExpression(data, assays = "pathways_mlm", features = rownames(data[["pathways_mlm"]]), group.by = "cell_type", verbose=FALSE)$pathways_mlm
  if(nrow(avg_pathway_activity_decoupler) > 0 && ncol(avg_pathway_activity_decoupler) > 0) {
      pathway_heatmap_path_decoupler <- file.path(base_output_path, "04c_pathway_activity_heatmap_decoupler.png")
      png(pathway_heatmap_path_decoupler, width=10, height=max(6, 0.3*nrow(avg_pathway_activity_decoupler) + 2), units="in", res=300)
      pheatmap::pheatmap(avg_pathway_activity_decoupler, scale = "row", 
                         main = "Mean Pathway Activity per Cell Type (PROGENy)",
                         fontsize_row = max(5, 10 - nrow(avg_pathway_activity_decoupler) * 0.1),
                         fontsize_col = 10)
      dev.off()
      message("DecoupleR Pathway (PROGENy) activity heatmap saved.")
  } else { message("Not enough data to plot PROGENy activity heatmap.")}
} else { message("Skipping PROGENy pathway activity due to network loading failure.")}
DefaultAssay(data) <- "RNA" # Reset default assay

intermediate_rds_path_04c <- file.path(base_output_path, "04c_data_after_decoupler.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04c)
  message(paste("Data object after DecoupleR analysis saved to:", intermediate_rds_path_04c))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.3:", e$message))
})

message("Finished Module 4.3: DecoupleR Analysis")
