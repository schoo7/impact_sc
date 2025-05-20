# IMPACT-sc Script: 04h_card_spatial.R
# Purpose: Annotate Spatial Data with CARD.
# THIS MODULE IS OPTIONAL AND REQUIRES A SPATIAL DATASET (e.g., Visium).

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

#### Module 4.8: Annotate Spatial Data with CARD ####
message("Starting Module 4.8: Annotate Spatial Data with CARD")

# Input: 
#   - '04g_data_after_queryproj.RDS' (or previous) as scRNA-seq REFERENCE ('sc_ref_data')
#   - SpatialExperiment RDS file (USER-PROVIDED, e.g., 'Acinar_Cell_Carcinoma.RDS')
# Output: CARD results, plots. Saves final processed Seurat object.

obj_prev_module_path <- file.path(base_output_path, "04g_data_after_queryproj.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS") # Fallback
  message(paste("04g output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required scRNA-seq reference Seurat object not found at:", obj_prev_module_path))
}
sc_ref_data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading scRNA-seq reference RDS:",e$message));NULL
})
if(is.null(sc_ref_data)) stop("Failed to load scRNA-seq reference data for Module 4.8.")
message("Loaded scRNA-seq reference object (as 'sc_ref_data') for Module 4.8.")
DefaultAssay(sc_ref_data) <- "RNA"


spatial_rds_path <- "../Acinar_Cell_Carcinoma.RDS" # Path to SpatialExperiment RDS, relative to 'scripts/'
if (!file.exists(spatial_rds_path)) {
  message(paste("SpatialExperiment RDS file not found at:", spatial_rds_path, ". Skipping Module 4.8."))
} else {
  spe <- tryCatch(readRDS(spatial_rds_path), error = function(e){
    message(paste("Error loading SpatialExperiment RDS:", e$message)); NULL
  })
  if(is.null(spe)){
    message("Failed to load SpatialExperiment object. Skipping CARD analysis.")
  } else {
    message("Loaded SpatialExperiment object.")
    
    # Prepare spatial data for CARD
    if (!("Spatial" %in% names(assays(spe))) || is.null(spatialCoords(spe))) {
        message("Spatial assay or coordinates not found in SpatialExperiment object. Skipping CARD.")
    } else {
        spatial_counts <- assays(spe)[["Spatial"]] # Or specific layer like counts(spe)
        if(is.null(spatial_counts)) spatial_counts <- counts(spe) # Fallback to counts if "Spatial" assay is complex

        spatial_coords_df <- as.data.frame(spatialCoords(spe))
        
        # CARD expects x, y. Check colnames and rename if necessary.
        if (!all(c("x", "y") %in% colnames(spatial_coords_df))) {
            if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(spatial_coords_df))) {
                spatial_coords_df <- dplyr::rename(spatial_coords_df, x = pxl_col_in_fullres, y = pxl_row_in_fullres)
            } else if (all(c("imagecol", "imagerow") %in% colnames(spatial_coords_df))) {
                spatial_coords_df <- dplyr::rename(spatial_coords_df, x = imagecol, y = imagerow)
            } else {
                stop("Spatial coordinates do not have 'x' and 'y' columns or common alternatives. Please check and rename.")
            }
        }
        spatial_location_df <- spatial_coords_df[, c("x", "y"), drop=FALSE]
        
        common_spots <- intersect(colnames(spatial_counts), rownames(spatial_location_df))
        if(length(common_spots) == 0) stop("No common spots/cells between spatial counts and coordinates.")
        spatial_counts <- spatial_counts[, common_spots, drop=FALSE]
        spatial_location_df <- spatial_location_df[common_spots, , drop=FALSE]

        # Prepare single-cell reference
        if (is.list(sc_ref_data@assays$RNA) || (length(sc_ref_data@assays$RNA@layers) > 0 && !is.null(names(sc_ref_data@assays$RNA@layers)))) { 
            sc_ref_data <- JoinLayers(sc_ref_data, assay="RNA") 
        }
        sc_counts <- GetAssayData(sc_ref_data, assay = "RNA", layer = "counts")
        
        if (!("cell_type" %in% colnames(sc_ref_data@meta.data))) stop("'cell_type' not found in scRNA-seq reference metadata.")
        if (!("orig.ident" %in% colnames(sc_ref_data@meta.data))) {
            message("orig.ident not found in scRNA-seq reference, using a dummy sample ID.")
            sc_ref_data$orig.ident <- "sc_sample1"
        }

        sc_meta_df <- data.frame(
          cellID = colnames(sc_ref_data),
          ident = sc_ref_data$cell_type,
          orig.ident = sc_ref_data$orig.ident 
        )
        rownames(sc_meta_df) <- sc_meta_df$cellID

        # Harmonize gene names (assuming SYMBOLs for both for simplicity here)
        common_genes_card <- intersect(rownames(sc_counts), rownames(spatial_counts))
        if(length(common_genes_card) < 50) { # Arbitrary threshold
            warning(paste("Very few common genes (", length(common_genes_card), ") between scRNA-seq and spatial data. CARD results may be unreliable. Check gene identifiers (e.g. SYMBOL vs ENSEMBL)."))
        }
        if(length(common_genes_card) == 0) stop("No common genes found between scRNA-seq and spatial data for CARD.")

        sc_counts_card <- sc_counts[common_genes_card, , drop=FALSE]
        spatial_counts_card <- spatial_counts[common_genes_card, , drop=FALSE]

        CARD_obj <- tryCatch(createCARDObject(
          sc_count = sc_counts_card, sc_meta = sc_meta_df,
          spatial_count = spatial_counts_card, spatial_location = spatial_location_df,
          ct.varname = "ident", ct.select = unique(sc_meta_df$ident),
          sample.varname = "orig.ident", minCountGene = 10, minCountSpot = 2 
        ), error = function(e){message(paste("Error creating CARD object:", e$message)); NULL})
        
        if(!is.null(CARD_obj)){
          CARD_obj <- tryCatch(CARD_deconvolution(CARD_object = CARD_obj),
                               error = function(e){message(paste("CARD deconvolution failed:", e$message)); NULL})
          if(!is.null(CARD_obj) && !is.null(CARD_obj@Proportion_CARD)){
            message("CARD deconvolution results (first 2 spots):")
            print(head(CARD_obj@Proportion_CARD, 2))
            
            # Visualization
            num_cell_types_card <- length(unique(sc_meta_df$ident))
            card_colors <- scales::hue_pal()(num_cell_types_card)
            
            p1_card_pie <- CARD.visualize.pie(
              proportion = CARD_obj@Proportion_CARD, spatial_location = CARD_obj@spatial_location,
              colors = card_colors[1:ncol(CARD_obj@Proportion_CARD)], radius = 30 # Adjust radius
            )
            ggsave(file.path(base_output_path, "04h_card_pie_plot.png"), plot = p1_card_pie, width = 8, height = 7)
            message("CARD pie plot saved.")
            
            saveRDS(CARD_obj, file.path(base_output_path, "04h_card_results.RDS"))
            message("CARD analysis complete and results saved.")
          } else { message("CARD deconvolution did not produce results.")}
        } else { message("Skipping CARD deconvolution due to object creation failure.")}
    }
  }
}

# This is the last script in the 04x sequence. Save the main 'sc_ref_data' object.
# It might have been modified by previous 04x scripts (e.g., UMAP model in 04g).
final_rds_path_module4 <- file.path(base_output_path, "04_module4_fully_processed_final.RDS")
tryCatch({
  saveRDS(sc_ref_data, final_rds_path_module4) # Save sc_ref_data as it's the main Seurat object
  message(paste("IMPACT-sc Module 4 processing complete. Final Seurat object saved to:", final_rds_path_module4))
}, error = function(e) {
  warning(paste("Error saving final Module 4 Seurat object:", e$message))
})

message("Finished Module 4.8: Annotate Spatial Data with CARD")
message("End of IMPACT-sc analysis workflow scripts.")
