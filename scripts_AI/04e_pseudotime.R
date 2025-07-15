# IMPACT-sc Script: 04e_pseudotime_analysis.R
# Purpose: Perform Pseudotime Analysis (e.g., with Palantir via SeuratExtend).
# THIS MODULE IS OPTIONAL AND DEPENDS ON SeuratExtend OR ALTERNATIVE TOOLS.

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

if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { library(org.Hs.eg.db); org_db_object <- org.Hs.eg.db } else { stop("org.Hs.eg.db not installed.") }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { library(org.Mm.eg.db); org_db_object <- org.Mm.eg.db } else { stop("org.Mm.eg.db not installed.") }
} else { stop("Species not supported.") }
options(future.globals.maxSize = 4 * 1024^3)
# --- End Script Parameters ---

#### Module 4.5: Pseudotime Analysis with SeuratExtend (Palantir) ####
message("Starting Module 4.5: Pseudotime Analysis (Palantir/SeuratExtend)")

# Input: '04d_data_after_ucell.RDS' (or previous if 04d skipped)
# Output: Modifies 'data' object with pseudotime results, saves plots.
obj_prev_module_path <- file.path(base_output_path, "04d_data_after_ucell.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS") # Fallback
   message(paste("04d output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required input Seurat object not found at:", obj_prev_module_path))
}
data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading input RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.5.")
message("Loaded Seurat object for Module 4.5.")
DefaultAssay(data) <- "RNA"

# Check if SeuratExtend functions for Palantir are available
if (exists("Palantir.RunDM") && exists("Palantir.Pseudotime") && 
    exists("GeneTrendCurve.Palantir") && exists("GeneTrendHeatmap.Palantir")) {
  message("Palantir functions (from SeuratExtend) found. Proceeding with pseudotime analysis.")
  
  data_palantir <- tryCatch(Palantir.RunDM(data), error = function(e){
    message(paste("Palantir.RunDM failed:", e$message)); NULL
  })
  
  if(!is.null(data_palantir) && "ms" %in% Reductions(data_palantir)) {
    p_dm_plot_seurat <- DimPlot(data_palantir, reduction = "ms", group.by="cell_type", label=TRUE)
    ggsave(file.path(base_output_path, "04f_palantir_diffusionmap.png"), plot = p_dm_plot_seurat, width=8, height=7)
    
    # **USER ACTION REQUIRED: Define a valid start cell for Palantir.**
    # This should ideally be passed as an environment variable from params.json
    start_cell_palantir <- Sys.getenv("IMPACT_SC_PALANTIR_START_CELL", colnames(data_palantir)[1]) # Default to first cell
    if(start_cell_palantir == colnames(data_palantir)[1]) {
        warning("IMPACT_SC_PALANTIR_START_CELL environment variable not set. Using first cell as placeholder. Please define a biologically relevant start cell in your params.json or manually.")
    }
    if(!(start_cell_palantir %in% colnames(data_palantir))) {
        stop(paste("Defined Palantir start cell '", start_cell_palantir, "' not found in Seurat object. Please check your parameter.", sep=""))
    }
    message(paste("Using start cell for Palantir:", start_cell_palantir))
    
    data_palantir <- tryCatch(Palantir.Pseudotime(data_palantir, start_cell = start_cell_palantir),
                              error = function(e){message(paste("Palantir.Pseudotime failed:", e$message)); NULL})
    
    if(!is.null(data_palantir) && !is.null(data_palantir@misc$Palantir$Pseudotime)) {
      ps_data <- data_palantir@misc$Palantir$Pseudotime
      # Ensure ps_data rownames match cell names in data_palantir for safe cbind
      ps_data_ordered <- ps_data[match(colnames(data_palantir), rownames(ps_data)), , drop=FALSE]
      data_palantir@meta.data <- cbind(data_palantir@meta.data, ps_data_ordered)
      
      if ("Pseudotime" %in% colnames(data_palantir@meta.data)) {
        p_pseudo_dm <- FeaturePlot(data_palantir, features = "Pseudotime", reduction = "ms") + scale_color_viridis_c()
        ggsave(file.path(base_output_path, "04f_palantir_pseudotime_on_dm.png"), plot = p_pseudo_dm, width=8, height=7)
      }
      data <- data_palantir # Update main data object
    } else { message("Palantir pseudotime calculation did not yield expected results.") }
  } else { message("Palantir.RunDM failed or did not produce 'ms' reduction.") }
} else {
  message("Skipping Pseudotime Analysis with Palantir as SeuratExtend functions are not available. Consider alternatives like Monocle3 or Slingshot.")
}

intermediate_rds_path_04f <- file.path(base_output_path, "04f_data_after_pseudotime.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04f)
  message(paste("Data object after (attempted) Pseudotime analysis saved to:", intermediate_rds_path_04f))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.5:", e$message))
})

message("Finished Module 4.5: Pseudotime Analysis")
