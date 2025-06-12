#
# IMPACT-sc: 04g_card
#
# Description:
# This script performs spatial deconvolution using the CARD algorithm.
# It integrates single-cell RNA-seq (scRNA-seq) data with spatial transcriptomics
# data to infer the proportions of different cell types at each spatial location.
#
# IMPORTANT: This module requires R version 4.3.0 or higher to run correctly.
#
# Inputs:
# 1. IMPACT_SC_SPATIAL_RDS_PATH: Environment variable pointing to the spatial data
#    Seurat object saved as an .RDS file. This object should contain spatial
#    transcriptomics data.
# 2. 03_module3_final_annotated.RDS: The annotated single-cell Seurat object from
#    module 03. This is used as the reference. Located in IMPACT_SC_OUTPUT_DIR.
#
# Outputs:
# - All outputs are saved directly into the user-specified output directory:
#   - 04g_CARD_obj.rds: The CARD object after deconvolution and imputation.
#   - 04g_CARD_proportions.csv: The inferred cell type proportions.
#   - 04g_plot_pie.png: Pie chart visualization of cell type proportions.
#   - 04g_plot_proportions.png: Spatial map of individual cell type proportions.
#   - 04g_plot_proportions_2CT.png: Co-localization plot for two selected cell types.
#   - 04g_plot_correlation.png: Heatmap of cell type proportion correlations.
#

# --- 1. Load libraries and check environment ---
cat("--- Loading libraries for CARD ---\n")
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(CARD)
})

# R version check
if (as.numeric(R.version$major) < 4 || (as.numeric(R.version$major) == 4 && as.numeric(R.version$minor) < 3.0)) {
    stop("FATAL: CARD module (04g) requires R version 4.3.0 or higher.")
}

# Get environment variables
output_dir <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", unset = ".")
spatial_rds_path <- Sys.getenv("IMPACT_SC_SPATIAL_RDS_PATH")
# Corrected path for the single-cell reference from module 03
sc_ref_rds_path <- file.path(output_dir, "03_module3_final_annotated.RDS")

# Validate inputs
if (spatial_rds_path == "" || !file.exists(spatial_rds_path)) {
    stop("FATAL: Spatial data RDS path not provided or file does not exist. Check IMPACT_SC_SPATIAL_RDS_PATH.")
}
if (!file.exists(sc_ref_rds_path)) {
    stop(paste("FATAL: Single-cell reference from module 03 (", sc_ref_rds_path, ") not found in the output directory."))
}

# --- 2. Prepare Spatial Data ---
cat("--- Preparing Spatial Data ---\n")
cat("Loading spatial data from:", spatial_rds_path, "\n")
spe <- readRDS(spatial_rds_path)
print("Spatial object summary:")
print(spe)

# Extract spatial count matrix
spatial_count <- LayerData(spe, assay = "Spatial", layer = "counts")
# Get spatial coordinates
coords_df <- GetTissueCoordinates(spe)
# Extract x and y coordinates
spatial_location <- as.matrix(coords_df[, c("x", "y")])
# Convert to data frame
spatial_location <- as.data.frame(spatial_location)
# Get spot/cell names from spatial data
location_cells <- rownames(spatial_location)
# Subset count matrix to ensure it only contains cells present in the spatial locations
spatial_count <- spatial_count[, colnames(spatial_count) %in% location_cells]

cat("Spatial data prepared. Dimensions of count matrix:", dim(spatial_count)[1], "genes,", dim(spatial_count)[2], "spots.\n")

# --- 3. Prepare Single-Cell Reference Data ---
cat("--- Preparing Single-Cell Reference Data ---\n")
cat("Loading scRNA-seq reference from:", sc_ref_rds_path, "\n")
sce <- readRDS(sc_ref_rds_path)
print("Single-cell reference object summary:")
print(sce)

# Join layers if needed (e.g., counts, logcounts)
if ("counts" %in% Layers(sce, search = "RNA") && "scale.data" %in% Layers(sce, search = "RNA")) {
    sce <- JoinLayers(sce)
}

# Extract single-cell count matrix
sc_count <- LayerData(sce, assay = "RNA", layer = "counts")
# Set cell type as identity label
if (!"cell_type" %in% colnames(sce@meta.data)) {
    stop("FATAL: 'cell_type' column not found in the metadata of the single-cell reference object.")
}
sce$ident <- sce$cell_type
# Extract metadata
metadata <- sce@meta.data
# Select sample and identity info
result <- metadata %>% select(orig.ident, ident)
# Add cell IDs
result$cellID <- rownames(result)
# Prepare metadata for CARD
sc_meta <- result

cat("Single-cell reference prepared. Dimensions of count matrix:", dim(sc_count)[1], "genes,", dim(sc_count)[2], "cells.\n")

# --- 4. Run CARD Deconvolution ---
cat("--- Running CARD ---\n")
# Create CARD object
CARD_obj <- createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "ident",
    ct.select = unique(sc_meta$ident),
    sample.varname = "orig.ident",
    minCountGene = 100,
    minCountSpot = 5
)

# Run CARD deconvolution
CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
cat("Deconvolution finished. Proportions matrix summary:\n")
print(CARD_obj@Proportion_CARD[1:min(5, nrow(CARD_obj@Proportion_CARD)), 1:min(5, ncol(CARD_obj@Proportion_CARD))])

# --- 5. Run Spatial Imputation ---
cat("--- Running CARD Imputation ---\n")
CARD_obj <- CARD.imputation(CARD_obj, NumGrids = 2000, ineibor = 10, exclude = NULL)
cat("Imputation finished.\n")

# Save the final CARD object and proportions directly to the main output directory
cat("--- Saving CARD results ---\n")
saveRDS(CARD_obj, file.path(output_dir, "04g_CARD_obj.rds"))
write.csv(CARD_obj@Proportion_CARD, file.path(output_dir, "04g_CARD_proportions.csv"))

# --- 6. Visualization ---
cat("--- Generating Visualizations ---\n")

# Visualization 1: Pie charts
p1 <- CARD.visualize.pie(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    radius = 50 
)
ggsave(filename = file.path(output_dir, "04g_plot_pie.png"), plot = p1, width = 8, height = 8, dpi = 300)
cat("Saved pie chart plot.\n")

# Visualization 2: Individual Proportions
# Select all cell types to visualize
ct.visualize_prop <- colnames(CARD_obj@Proportion_CARD) 
p2 <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    ct.visualize = ct.visualize_prop,
    colors = c("lightblue", "lightyellow", "red"),
    NumCols = 4,
    pointSize = 0.2
)
ggsave(filename = file.path(output_dir, "04g_plot_proportions.png"), plot = p2, width = 12, height = 3 * ceiling(length(ct.visualize_prop)/4), dpi = 300, limitsize = FALSE)
cat("Saved individual proportion plots.\n")

# Visualization 3: Two Cell Type Co-localization (if more than one cell type)
if (ncol(CARD_obj@Proportion_CARD) >= 2) {
    ct2.visualize_prop <- colnames(CARD_obj@Proportion_CARD)[1:2]
    p3 <- CARD.visualize.prop.2CT(
        proportion = CARD_obj@Proportion_CARD,
        spatial_location = CARD_obj@spatial_location,
        ct2.visualize = ct2.visualize_prop,
        colors = list(c("lightblue", "lightyellow", "red"), c("lightblue", "lightyellow", "black"))
    )
    ggsave(filename = file.path(output_dir, "04g_plot_proportions_2CT.png"), plot = p3, width = 10, height = 5, dpi = 300)
    cat("Saved 2-celltype co-localization plot.\n")
} else {
    cat("Skipping 2-celltype plot; only one cell type present.\n")
}


# Visualization 4: Correlation Heatmap
p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD, colors = NULL)
ggsave(filename = file.path(output_dir, "04g_plot_correlation.png"), plot = p4, width = 7, height = 6, dpi = 300)
cat("Saved correlation heatmap.\n")

cat("--- Module 04g: CARD finished successfully ---\n")

