# IMPACT-sc Script: 02c_load_c2s_result.R
# Purpose: Load Cell2Sentence embeddings and predictions back into a Seurat object.
# This script is designed to run after '02b_c2s.py' has generated the output files.
# Version: 2.0 (Modified for robust prediction file reading)

# --- Libraries ---
library(AnnotationDbi)
library(decoupleR) 
library(dplyr)
library(ensembldb) 
library(ggplot2)
library(ggpubr)
library(homologene) 
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
# These parameters are configured by 'interactive_setup.py' and passed by the main pipeline script.
species <- Sys.getenv("IMPACT_SC_SPECIES", "human")
set.seed(1234)

# The base output directory for all pipeline results.
base_output_path <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", "demo_output")
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}
message(paste("Base output directory set to:", base_output_path))

# Load species-specific annotation database
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { 
    library(org.Hs.eg.db)
    org_db_object <- org.Hs.eg.db 
  } else { 
    stop("The 'org.Hs.eg.db' annotation package is not installed. Please install it to proceed.") 
  }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { 
    library(org.Mm.eg.db)
    org_db_object <- org.Mm.eg.db 
  } else { 
    stop("The 'org.Mm.eg.db' annotation package is not installed. Please install it to proceed.") 
  }
} else { 
  stop("Species not supported. Please use 'human' or 'mouse'.") 
}
# --- End Script Parameters ---


#### Module 2c: Cell2Sentence - Load Results ####
message("Starting Module 2c: Loading Cell2Sentence Results into R")

# --- Path Definitions ---
# This section defines the expected locations for input and output files.

# Path to the Seurat object that was prepped for C2S in module 02a.
obj_prepped_for_c2s_path <- file.path(base_output_path, "02_module2_c2s_prepped_object.RDS")

# --- HARDCODED PATHS (QUICK FIX) ---
# NOTE: This is a temporary fix based on user feedback. The robust solution 
# is to define these paths in interactive_setup.py and read them from environment variables.
# We construct the path by appending the known subdirectory to the base_output_path.
c2s_embeddings_csv_path <- file.path(base_output_path, "cell2sentence_module_outputs", "c2s_embeddings.csv")
c2s_predicted_csv_path <- file.path(base_output_path, "cell2sentence_module_outputs", "c2s_predictions.csv")

message(paste("INFO: Using hardcoded path for embeddings:", c2s_embeddings_csv_path))
message(paste("INFO: Using hardcoded path for predictions:", c2s_predicted_csv_path))
# --- END OF HARDCODED PATHS ---

# Path for the final output RDS file of this script.
module2_c2s_processed_path <- file.path(base_output_path, "02_module2_c2s_processed.RDS")
# --- End Path Definitions ---


# Load the Seurat object that was used to generate the H5AD file for Python
if (file.exists(obj_prepped_for_c2s_path)) {
  obj <- tryCatch(
    readRDS(obj_prepped_for_c2s_path), 
    error = function(e) {
      message(paste("Error loading '02_module2_c2s_prepped_object.RDS':", e$message))
      return(NULL)
    }
  )
  if(is.null(obj)) stop("Failed to load the prepped Seurat object for C2S results.")
} else {
  stop("Prepped Seurat object ('02_module2_c2s_prepped_object.RDS') not found at: ", obj_prepped_for_c2s_path)
}
DefaultAssay(obj) <- "RNA" 

# Load Cell2Sentence Embeddings
if (!file.exists(c2s_embeddings_csv_path)) {
  stop(paste("CRITICAL ERROR: C2S cell embeddings CSV file not found at the hardcoded path:", c2s_embeddings_csv_path))
} else {
  message(paste("Loading C2S embeddings from:", c2s_embeddings_csv_path))
  embedding_df <- tryCatch(
    read.csv(c2s_embeddings_csv_path, row.names = 1), # row.names = 1 because cell IDs are the first column from Python
    error = function(e){
      stop(paste("Error reading embeddings CSV:", e$message))
    }
  )
  
  if (!is.null(embedding_df)) {
    embedding_matrix <- as.matrix(embedding_df)
    
    # Ensure cell IDs in the embedding matrix match the Seurat object
    if (!all(rownames(embedding_matrix) %in% colnames(obj))) {
        warning("Mismatch or missing cell IDs between C2S embeddings and Seurat object. Check Python script output. Attempting to use common cells.")
        common_cells_emb <- intersect(rownames(embedding_matrix), colnames(obj))
        if (length(common_cells_emb) == 0) stop("No common cells found between embeddings and Seurat object.")
        
        # Subset both the object and the matrix to only the common cells
        obj <- subset(obj, cells = common_cells_emb) 
        embedding_matrix <- embedding_matrix[common_cells_emb, , drop = FALSE]
    } else {
       # Ensure the matrix is in the same order as the object
       embedding_matrix <- embedding_matrix[colnames(obj), , drop = FALSE] 
    }

    # Add the embeddings as a new dimensional reduction object
    colnames(embedding_matrix) <- paste0("C2S_", seq_len(ncol(embedding_matrix)))
    c2s_reduction <- CreateDimReducObject(
      embeddings = embedding_matrix, 
      key = "C2S_", 
      assay = DefaultAssay(obj)
    )
    obj[["c2s_embeddings"]] <- c2s_reduction
    message("C2S embeddings successfully added to the Seurat object.")

    # Run UMAP on the new C2S embeddings
    num_c2s_dims <- ncol(embedding_matrix)
    if (num_c2s_dims > 0) {
        obj <- RunUMAP(obj, reduction = "c2s_embeddings", dims = 1:num_c2s_dims, 
                       reduction.name = "umap_c2s", min.dist = 0.1, verbose = FALSE)
        message("UMAP run on C2S embeddings and stored as 'umap_c2s'.")
    } else {
        warning("No C2S embedding dimensions were found to run UMAP.")
    }
  }
}

# --- MODIFIED SECTION: Load Predicted Cell Types from Cell2Sentence ---
if (file.exists(c2s_predicted_csv_path)) {
  message(paste("Loading C2S predicted cell types from:", c2s_predicted_csv_path))
  predicted_df_c2s <- tryCatch(
    read.csv(c2s_predicted_csv_path, header = TRUE, stringsAsFactors = FALSE), 
    error = function(e){
      warning(paste("Error reading predicted cell types CSV:", e$message))
      return(NULL)
    }
  )

  # The predictions CSV should have at least two columns:
  # 1. Cell IDs (barcodes)
  # 2. Predicted labels
  # We will assume the first column contains the cell IDs and the second contains the predictions,
  # regardless of the column headers, to make the script more robust.
  if (!is.null(predicted_df_c2s) && ncol(predicted_df_c2s) >= 2) {
      message("INFO: Assuming first column of predictions CSV contains cell IDs and second column contains predicted labels.")

      # Rename columns for clarity and consistency
      colnames(predicted_df_c2s)[1] <- "cell_id"
      colnames(predicted_df_c2s)[2] <- "predicted_label"

      cell_ids_from_csv <- predicted_df_c2s$cell_id
      
      # Check if the cell IDs from the CSV match the Seurat object
      if (all(cell_ids_from_csv %in% colnames(obj))) {
          # Create a named vector for easy mapping
          c2s_preds_map <- setNames(predicted_df_c2s$predicted_label, cell_ids_from_csv)
          
          # Add the predictions to the Seurat object metadata
          # The map is applied to the object's cell names to ensure correct order
          obj$c2s_predicted_label <- c2s_preds_map[colnames(obj)]
          
          message("C2S predicted labels successfully added to Seurat object metadata as 'c2s_predicted_label'.")
      } else {
          warning("Cell IDs in C2S predictions CSV do not fully match the current Seurat object. Prediction loading skipped.")
      }
  } else {
      warning("C2S predictions CSV does not have at least two columns (for cell IDs and predictions) or could not be read.")
  }
} else {
    warning(paste("WARNING: C2S predicted cell types CSV file not found at the hardcoded path:", c2s_predicted_csv_path))
}
# --- END OF MODIFIED SECTION ---


# Save the final Seurat object with the new C2S data
tryCatch({
  saveRDS(obj, module2_c2s_processed_path)
  message("Module 2c (C2S results loaded) output saved to: ", module2_c2s_processed_path)
}, error = function(e) {
  warning(paste("Error saving final Module 2c output:", e$message))
})

message("Finished Module 2c: Loading Cell2Sentence Results.")
