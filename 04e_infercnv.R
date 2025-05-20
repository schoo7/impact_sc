# IMPACT-sc Script: 04e_infercnv_analysis.R
# Purpose: Predict CNV Activity with InferCNV.

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

#### Module 4.5: Predict CNV Activity with InferCNV ####
message("Starting Module 4.5: InferCNV Analysis")

# Input: '04d_data_after_ucell.RDS' (or previous if 04d skipped)
# Output: InferCNV results in a subdirectory, potentially modifies 'data' object.
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

ref_group_names_infercnv <- "Normal_Epithelial" # MODIFY THIS to actual normal cell type(s) in data$cell_type
message("For InferCNV, ensure 'ref_group_names_infercnv' is correctly set for your data: ", paste(ref_group_names_infercnv, collapse=", "))

# **USER ACTION REQUIRED FOR GENE ORDER FILE**
gene_order_file_path_m4.5 <- if (species == "human") {
  "../hg38_gencode_v27.txt" # Example for human, path relative to 'scripts/' dir
} else if (species == "mouse") {
  "../mm10_gencode_vM25.txt" # Example for mouse, path relative to 'scripts/' dir
} else {
  stop("Species not supported for InferCNV gene order file example.")
}
message(paste("IMPORTANT: The gene_order_file for InferCNV is set to '", gene_order_file_path_m4.5, 
              "'. You MUST ensure this file is correct for species '", species, 
              "' and is present in your project directory (one level above 'scripts/') or provide the full path.", sep=""))

if (!file.exists(gene_order_file_path_m4.5)) {
  warning(paste("Gene order file '", gene_order_file_path_m4.5, "' not found. InferCNV may not run correctly or gene ordering will be arbitrary."))
}

if (!(any(ref_group_names_infercnv %in% unique(data$cell_type)))) {
  warning(paste("Reference group(s) for InferCNV ('", paste(ref_group_names_infercnv, collapse=","), "') not found in data$cell_type. Skipping InferCNV.", sep=""))
} else if (!file.exists(gene_order_file_path_m4.5)) {
  warning(paste("Gene order file not found for InferCNV. Skipping InferCNV."))
} else {
  infercnv_out_dir_m4.5 <- file.path(base_output_path, "04e_InferCNV_output/")
  if (!dir.exists(infercnv_out_dir_m4.5)) dir.create(infercnv_out_dir_m4.5, recursive = TRUE)

  counts_matrix_infercnv_m4.5 <- GetAssayData(data, assay = "RNA", layer = "counts") 
  counts_file_path_m4.5 <- file.path(infercnv_out_dir_m4.5, "counts_matrix.txt")
  tryCatch(write.table(as.matrix(counts_matrix_infercnv_m4.5), file = counts_file_path_m4.5, quote = FALSE, sep = "\t", col.names = NA),
           error = function(e) warning(paste("Error writing InferCNV counts matrix:", e$message)))

  annotations_df_infercnv_m4.5 <- data.frame(CellID = colnames(data), CellType = data$cell_type)
  annotations_file_path_m4.5 <- file.path(infercnv_out_dir_m4.5, "cellAnnotations.txt")
  tryCatch(write.table(annotations_df_infercnv_m4.5, file = annotations_file_path_m4.5, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE),
           error = function(e) warning(paste("Error writing InferCNV annotations file:", e$message)))

  infercnv_obj_m4.5 <- tryCatch(CreateInfercnvObject(
    raw_counts_matrix = counts_file_path_m4.5, annotations_file = annotations_file_path_m4.5, delim = "\t",
    gene_order_file = gene_order_file_path_m4.5, ref_group_names = ref_group_names_infercnv
  ), error = function(e) {message(paste("Error creating InferCNV object:", e$message)); NULL})
  
  if(!is.null(infercnv_obj_m4.5)){
      num_threads_infercnv <- max(1, floor(parallel::detectCores() / 2))
      message(paste("Running InferCNV with", num_threads_infercnv, "threads. This may take a long time."))
      infercnv_run_obj_m4.5 <- tryCatch({
        infercnv::run(
          infercnv_obj_m4.5, cutoff = 0.1, out_dir = infercnv_out_dir_m4.5,
          cluster_by_groups = TRUE, denoise = TRUE, HMM = TRUE, 
          num_threads = num_threads_infercnv, analysis_mode = "cells"
        )
      }, error = function(e) { message(paste("InferCNV run failed:", e$message)); return(NULL) })

      if (!is.null(infercnv_run_obj_m4.5)) {
        message("InferCNV analysis completed. Output in: ", infercnv_out_dir_m4.5)
        # To add results back to Seurat, use infercnv::add_to_seurat or parse output files.
        # data <- infercnv::add_to_seurat( infercnv_output_path = infercnv_out_dir_m4.5, seurat_obj = data, top_n = 0 )
        # message("InferCNV results added to Seurat object (example, needs verification of output).")
      } else { message("InferCNV run was not successful.") }
  } else { message("Skipping InferCNV run due to object creation failure.")}
}

intermediate_rds_path_04e <- file.path(base_output_path, "04e_data_after_infercnv.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04e)
  message(paste("Data object after (attempted) InferCNV analysis saved to:", intermediate_rds_path_04e))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.5:", e$message))
})

message("Finished Module 4.5: InferCNV Analysis")
