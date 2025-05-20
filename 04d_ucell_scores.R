# IMPACT-sc Script: 04d_ucell_genescores.R
# Purpose: Calculate Gene Scores with UCell.

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

#### Module 4.4: Calculate Gene Scores with UCell ####
message("Starting Module 4.4: UCell Gene Score Calculation")

# Input: '04c_data_after_decoupler.RDS' (or previous if 04c skipped)
# Output: Modifies 'data' object by adding UCell scores, saves plot.
obj_prev_module_path <- file.path(base_output_path, "04c_data_after_decoupler.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS") # Fallback
  message(paste("04c output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required input Seurat object not found at:", obj_prev_module_path))
}
data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading input RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.4.")
message("Loaded Seurat object for Module 4.4.")
DefaultAssay(data) <- "RNA"

msigdbr_species_param <- if (species == "human") "Homo sapiens" else if (species == "mouse") "Mus musculus" else "Homo sapiens"
message(paste("Using species '", msigdbr_species_param, "' for MSigDB (UCell).", sep=""))

hallmark_sets_ucell <- tryCatch(msigdbr(species = msigdbr_species_param, category = "H"), 
                                error=function(e){message(paste("Error getting MSigDB Hallmark sets:", e$message));NULL})

if(!is.null(hallmark_sets_ucell) && nrow(hallmark_sets_ucell) > 0){
  hallmark_list_ucell <- hallmark_sets_ucell %>%
    dplyr::select(gs_name, gene_symbol) %>% 
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(gene_symbols = list(gene_symbol), .groups = "drop") %>% 
    tibble::deframe()
  
  # Ensure 'data' slot has log-normalized counts for UCell
  if (!("data" %in% slotNames(data@assays$RNA))) { 
    message("Normalizing data as 'data' slot is missing for UCell.")
    data <- NormalizeData(data, verbose=FALSE) 
  }
  data <- AddModuleScore_UCell(data, features = hallmark_list_ucell, assay = "RNA", slot = "data", name = "_UCell", BPPARAM = BiocParallel::SerialParam()) # Added BPPARAM for non-parallel execution
  
  ucell_score_cols_full_ucell <- paste0(names(hallmark_list_ucell), "_UCell")
  ucell_score_cols_present_ucell <- ucell_score_cols_full_ucell[ucell_score_cols_full_ucell %in% colnames(data@meta.data)]

  reduction_for_smoothing_ucell <- "harmony" 
  if (!(reduction_for_smoothing_ucell %in% Reductions(data))) {
      reduction_for_smoothing_ucell <- "pca" # Fallback
      if(!(reduction_for_smoothing_ucell %in% Reductions(data))) {
          warning(paste("Reduction", reduction_for_smoothing_ucell, "not found for UCell KNN smoothing. Skipping smoothing."))
          reduction_for_smoothing_ucell <- NULL # Explicitly nullify
      }
  }

  if (length(ucell_score_cols_present_ucell) > 0 && !is.null(reduction_for_smoothing_ucell)) {
      data <- SmoothKNN(obj = data, signature.names = ucell_score_cols_present_ucell, reduction = reduction_for_smoothing_ucell)
      message("UCell scores smoothed using KNN.")
      smoothed_ucell_names_ucell <- paste0(ucell_score_cols_present_ucell, "_kNN")
      
      example_hallmark_set_smoothed_ucell <- head(smoothed_ucell_names_ucell[smoothed_ucell_names_ucell %in% colnames(data@meta.data)], 1)
      umap_to_use_viz_m4_ucell <- if ("umap_c2s" %in% Reductions(data)) "umap_c2s" else if ("umap" %in% Reductions(data)) "umap" else "pca"
      if (length(example_hallmark_set_smoothed_ucell) > 0 ) {
        p_ucell_smooth <- FeaturePlot(data, features = example_hallmark_set_smoothed_ucell, reduction = umap_to_use_viz_m4_ucell)
        ggsave(file.path(base_output_path, "04d_ucell_example_hallmark_smooth_umap.png"), plot = p_ucell_smooth, width = 7, height = 6)
        message("Example smoothed UCell score UMAP saved.")
      }
  } else { message("No UCell scores to smooth or reduction for smoothing not found.") }
} else { message("Failed to get Hallmark gene sets from MSigDB. Skipping UCell.")}

intermediate_rds_path_04d <- file.path(base_output_path, "04d_data_after_ucell.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04d)
  message(paste("Data object after UCell analysis saved to:", intermediate_rds_path_04d))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.4:", e$message))
})

message("Finished Module 4.4: UCell Gene Score Calculation")
