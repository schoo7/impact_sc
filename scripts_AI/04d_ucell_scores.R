# IMPACT-sc Script: 04d_ucell_genescores.R
# Purpose: Calculate Gene Scores with UCell.

# --- Libraries ---
library(AnnotationDbi)
library(CARD)
# library(celldex) # Celldex is no longer used for downloading reference data
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
if (!dir.exists(base_output_path)) dir.create(base_output_path, recursive = TRUE)
message(paste("Output directory set to:", base_output_path))

# New environment variables for UCell customization
msigdb_category_to_use <- Sys.getenv("IMPACT_SC_MSIGDB_CATEGORY", "H") # Default to Hallmark if not set
ucell_plot_pathway_name_env <- Sys.getenv("IMPACT_SC_UCELL_PLOT_PATHWAY_NAME", "") # Default to empty (plot first available)

message(paste("MSigDB category for UCell set to:", msigdb_category_to_use))
if (ucell_plot_pathway_name_env != "") {
  message(paste("Requested UCell pathway for plotting:", ucell_plot_pathway_name_env))
} else {
  message("No specific UCell pathway requested for plotting, will use the first available.")
}

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

msigdbr_species_param <- if (species == "human") "Homo sapiens" else if (species == "mouse") "Mus musculus" else "Homo sapiens" # Default to human if species is weird
message(paste("Using species '", msigdbr_species_param, "' for MSigDB (UCell).", sep=""))

# Use the user-defined MSigDB category
gene_sets_ucell <- tryCatch(msigdbr(species = msigdbr_species_param, category = msigdb_category_to_use),
                                error=function(e){message(paste("Error getting MSigDB sets for category '", msigdb_category_to_use, "':", e$message));NULL})

if(!is.null(gene_sets_ucell) && nrow(gene_sets_ucell) > 0){
  gene_list_ucell <- gene_sets_ucell %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(gene_symbols = list(gene_symbol), .groups = "drop") %>%
    tibble::deframe()

  # Ensure 'data' slot has log-normalized counts for UCell
  if (!("data" %in% slotNames(data@assays$RNA))) {
    message("Normalizing data as 'data' slot is missing for UCell.")
    data <- NormalizeData(data, verbose=FALSE)
  }
  data <- AddModuleScore_UCell(data, features = gene_list_ucell, assay = "RNA", slot = "data", name = "_UCell", BPPARAM = BiocParallel::SerialParam())

  ucell_score_cols_full_ucell <- paste0(names(gene_list_ucell), "_UCell")
  ucell_score_cols_present_ucell <- ucell_score_cols_full_ucell[ucell_score_cols_full_ucell %in% colnames(data@meta.data)]

  reduction_for_smoothing_ucell <- "harmony"
  if (!(reduction_for_smoothing_ucell %in% Reductions(data))) {
      reduction_for_smoothing_ucell <- "pca"
      if(!(reduction_for_smoothing_ucell %in% Reductions(data))) {
          warning(paste("Reduction", reduction_for_smoothing_ucell, "not found for UCell KNN smoothing. Skipping smoothing."))
          reduction_for_smoothing_ucell <- NULL
      }
  }

  if (length(ucell_score_cols_present_ucell) > 0 && !is.null(reduction_for_smoothing_ucell)) {
      data <- SmoothKNN(obj = data, signature.names = ucell_score_cols_present_ucell, reduction = reduction_for_smoothing_ucell)
      message("UCell scores smoothed using KNN.")
      smoothed_ucell_names_ucell <- paste0(ucell_score_cols_present_ucell, "_kNN") # These are the names of the smoothed scores

      pathway_to_plot_ucell_smooth <- NULL
      # Check if user specified a pathway and if it exists (after smoothing suffix is added)
      if (ucell_plot_pathway_name_env != "") {
          potential_plot_name_smooth <- paste0(ucell_plot_pathway_name_env, "_UCell_kNN")
          if (potential_plot_name_smooth %in% smoothed_ucell_names_ucell) {
              pathway_to_plot_ucell_smooth <- potential_plot_name_smooth
              message(paste("Found requested pathway for plotting (smoothed):", pathway_to_plot_ucell_smooth))
          } else {
              # Check if the raw (unsmoothed) name was provided and if its smoothed version exists
              potential_raw_plot_name_smooth <- paste0(gsub("_UCell$", "", ucell_plot_pathway_name_env), "_UCell_kNN") # try to match if user gave raw name
               if (potential_raw_plot_name_smooth %in% smoothed_ucell_names_ucell) {
                  pathway_to_plot_ucell_smooth <- potential_raw_plot_name_smooth
                  message(paste("Found requested pathway by matching raw name to smoothed version:", pathway_to_plot_ucell_smooth))
              } else {
                message(paste("Requested pathway '", ucell_plot_pathway_name_env, "' (or its smoothed version '", potential_plot_name_smooth, "') not found among smoothed UCell scores. Will plot the first available.", sep=""))
              }
          }
      }

      # If no specific pathway is chosen or found, pick the first available smoothed score
      if (is.null(pathway_to_plot_ucell_smooth) || pathway_to_plot_ucell_smooth == "") {
          if (length(smoothed_ucell_names_ucell) > 0) {
            pathway_to_plot_ucell_smooth <- head(smoothed_ucell_names_ucell, 1)
            message(paste("Plotting the first available smoothed UCell score:", pathway_to_plot_ucell_smooth))
          } else {
            message("No smoothed UCell scores available to plot.")
          }
      }

      umap_to_use_viz_m4_ucell <- if ("umap_c2s" %in% Reductions(data)) "umap_c2s" else if ("umap" %in% Reductions(data)) "umap" else "pca"

      if (!is.null(pathway_to_plot_ucell_smooth) && pathway_to_plot_ucell_smooth != "" && (pathway_to_plot_ucell_smooth %in% colnames(data@meta.data))) {
        p_ucell_smooth <- FeaturePlot(data, features = pathway_to_plot_ucell_smooth, reduction = umap_to_use_viz_m4_ucell)
        plot_filename <- paste0("04d_ucell_featureplot_", gsub("[^A-Za-z0-9_]", "_", pathway_to_plot_ucell_smooth), "_", umap_to_use_viz_m4_ucell, ".png") # Sanitize filename
        ggsave(file.path(base_output_path, plot_filename), plot = p_ucell_smooth, width = 7, height = 6)
        message(paste("UCell score UMAP for '", pathway_to_plot_ucell_smooth, "' saved to ", plot_filename, sep=""))
      } else {
        message(paste("Could not generate UMAP plot. Pathway '", pathway_to_plot_ucell_smooth, "' not found in data metadata or no pathway selected.", sep=""))
      }

  } else { message("No UCell scores to smooth or reduction for smoothing not found.")}
} else { message(paste("Failed to get gene sets from MSigDB for category '", msigdb_category_to_use, "'. Skipping UCell.", sep=""))}

intermediate_rds_path_04d <- file.path(base_output_path, "04d_data_after_ucell.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04d)
  message(paste("Data object after UCell analysis saved to:", intermediate_rds_path_04d))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.4:", e$message))
})

message("Finished Module 04d: UCell Gene Score Calculation")
