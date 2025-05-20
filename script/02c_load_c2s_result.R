# IMPACT-sc Script: 02c_load_c2s_results.R
# Purpose: Load Cell2Sentence embeddings and predictions back into Seurat.

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
# --- End Script Parameters ---

#### Module 2.2: Cell2Sentence - Load Results (R part) ####
message("Starting Module 2c: Loading Cell2Sentence Results into R")

# Input: 
#   - '02_module2_c2s_prepped_object.RDS' (Seurat object before H5AD conversion)
#   - 'c2s_cell_embeddings.csv' and 'predicted_cell_types.csv' from C2S_output/
# Output: Saves '02_module2_c2s_processed.RDS' with C2S embeddings and UMAP.

# Load the Seurat object that was used to generate the H5AD file
obj_prepped_for_c2s_path <- file.path(base_output_path, "02_module2_c2s_prepped_object.RDS")
if (file.exists(obj_prepped_for_c2s_path)) {
  obj <- tryCatch(readRDS(obj_prepped_for_c2s_path), error = function(e) {
    message(paste("Error loading 02_module2_c2s_prepped_object.RDS:", e$message)); NULL
  })
  if(is.null(obj)) stop("Failed to load prepped Seurat object for C2S results.")
} else {
  stop("Prepped Seurat object (02_module2_c2s_prepped_object.RDS) not found at: ", obj_prepped_for_c2s_path)
}
DefaultAssay(obj) <- "RNA" # Ensure correct default assay

# Get paths from environment variables (should be set if 02a was run in same session, or set manually)
c2s_output_dir_env <- file.path(base_output_path, "C2S_output") # Reconstruct if not from env
c2s_embeddings_csv_path <- Sys.getenv("C2S_EMBEDDINGS_CSV", file.path(c2s_output_dir_env, "c2s_cell_embeddings.csv"))
c2s_predicted_csv_path <- Sys.getenv("C2S_PREDICTED_CSV", file.path(c2s_output_dir_env, "predicted_cell_types.csv"))


if (!file.exists(c2s_embeddings_csv_path)) {
  warning(paste("C2S cell embeddings CSV file not found at:", c2s_embeddings_csv_path, ". Skipping embedding loading."))
} else {
  message(paste("Loading C2S embeddings from:", c2s_embeddings_csv_path))
  embedding_df <- tryCatch(read.csv(c2s_embeddings_csv_path, row.names = 1), error = function(e){
    warning(paste("Error reading embeddings CSV:", e$message)); NULL
  })
  
  if (!is.null(embedding_df)) {
    embedding_matrix <- as.matrix(embedding_df)
    if (!all(rownames(embedding_matrix) %in% colnames(obj))) {
        warning("Mismatch or missing cell IDs between C2S embeddings and Seurat object. Check Python script output. Attempting to use common cells.")
        common_cells_emb <- intersect(rownames(embedding_matrix), colnames(obj))
        if (length(common_cells_emb) == 0) stop("No common cells found between embeddings and Seurat object.")
        obj <- subset(obj, cells = common_cells_emb) # Subset Seurat object to common cells
        embedding_matrix <- embedding_matrix[common_cells_emb, , drop = FALSE]
    } else {
       embedding_matrix <- embedding_matrix[colnames(obj), , drop = FALSE] # Reorder to match Seurat object
    }

    colnames(embedding_matrix) <- paste0("C2S_", seq_len(ncol(embedding_matrix)))
    c2s_reduction <- CreateDimReducObject(
      embeddings = embedding_matrix, 
      key = "C2S_", 
      assay = DefaultAssay(obj)
    )
    obj[["c2s_embeddings"]] <- c2s_reduction
    message("C2S embeddings added to Seurat object.")

    num_c2s_dims <- ncol(embedding_matrix)
    if (num_c2s_dims > 0) {
        obj <- RunUMAP(obj, reduction = "c2s_embeddings", dims = 1:num_c2s_dims, 
                       reduction.name = "umap_c2s", min.dist = 0.1, verbose = FALSE)
        message("UMAP run on C2S embeddings, stored as 'umap_c2s'.")
    } else {
        warning("No C2S embedding dimensions found to run UMAP.")
    }
  }
}

# Optionally, load predicted cell types from C2S (if prediction was successful)
if (!file.exists(c2s_predicted_csv_path)) {
  warning(paste("C2S predicted cell types CSV file not found at:", c2s_predicted_csv_path, ". Skipping prediction loading."))
} else {
  message(paste("Loading C2S predicted cell types from:", c2s_predicted_csv_path))
  predicted_df_c2s <- tryCatch(read.csv(c2s_predicted_csv_path, header = TRUE, stringsAsFactors = FALSE), error = function(e){
    warning(paste("Error reading predicted types CSV:", e$message)); NULL
  })

  if (!is.null(predicted_df_c2s) && nrow(predicted_df_c2s) > 0) {
    # Assuming predicted_df_c2s has 'text_ids' (cell barcodes) and 'predicted_label'
    if ("text_ids" %in% colnames(predicted_df_c2s) && "predicted_label" %in% colnames(predicted_df_c2s)) {
        if (all(predicted_df_c2s$text_ids %in% colnames(obj))) {
            obj_cell_order <- colnames(obj)
            c2s_preds_ordered <- character(length(obj_cell_order))
            names(c2s_preds_ordered) <- obj_cell_order
            
            c2s_preds_map <- setNames(predicted_df_c2s$predicted_label, predicted_df_c2s$text_ids)
            common_cells_preds <- intersect(obj_cell_order, names(c2s_preds_map))
            c2s_preds_ordered[common_cells_preds] <- c2s_preds_map[common_cells_preds]
            
            obj$c2s_direct_predicted_label <- c2s_preds_ordered # Add as a new metadata column
            message("C2S direct predicted labels added to Seurat object metadata as 'c2s_direct_predicted_label'.")
        } else {
            warning("Cell IDs in C2S predicted types do not fully match current Seurat object. Skipping direct prediction loading.")
        }
    } else {
        warning("C2S predicted types CSV does not have expected 'text_ids' or 'predicted_label' columns.")
    }
  } else {
    warning("C2S predicted cell types CSV is empty or could not be read.")
  }
}


module2_c2s_processed_path <- file.path(base_output_path, "02_module2_c2s_processed.RDS")
tryCatch({
  saveRDS(obj, module2_c2s_processed_path)
  message("Module 2c (C2S results loaded) output saved to: ", module2_c2s_processed_path)
}, error = function(e) {
  warning(paste("Error saving Module 2c output:", e$message))
})

message("Finished Module 2c: Loading Cell2Sentence Results.")
