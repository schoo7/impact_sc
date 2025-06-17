#
# IMPACT-SC Single-Cell Analysis Pipeline
#
# Module 04h: Cell-Cell Communication Analysis with LIANA
#
# Version: 1.5.0
#
# Description:
# This module performs cell-cell communication analysis using the LIANA+ framework.
# It reads the annotated Seurat object from Module 3, runs the LIANA workflow
# with one or more user-specified methods. It conditionally aggregates results
# only when multiple methods are used, preventing errors with single-method runs.
# It uses liana's native plotting functions.
# All plots are saved to PNG files using a standardized png()-dev.off() approach.
#
# Environment Variables:
#   - IMPACT_SC_OUTPUT_DIR: The directory where output files will be saved.
#   - IMPACT_SC_CELLCHAT_SOURCE_GROUPS: Comma-separated string of source cell type clusters (e.g., "0,1,2").
#   - IMPACT_SC_CELLCHAT_TARGET_GROUPS: Comma-separated string of target cell type clusters (e.g., "3,4,5").
#   - IMPACT_SC_LIANA_METHOD: Comma-separated string of LIANA methods to use (e.g., "natmi,logfc").
#

# --- 1. Load Libraries and Set Up ---
# -------------------------------------
# Use suppressPackageStartupMessages to keep the console clean.
suppressPackageStartupMessages({
  library(liana)
  library(tidyverse)
  library(readr)
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(magrittr)
  library(circlize) # Dependency for chord_freq
})

# --- 2. Get Parameters from Environment Variables ---
# ---------------------------------------------------
# Retrieve the output directory from the environment variable set by the orchestrator script.
output_dir <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", unset = "results")
# Construct the path to the input Seurat object from Module 3.
seurat_object_path <- file.path(output_dir, "03_module3_final_annotated.RDS")

# Retrieve source and target group strings.
source_groups_str <- Sys.getenv("IMPACT_SC_CELLCHAT_SOURCE_GROUPS", unset = "")
target_groups_str <- Sys.getenv("IMPACT_SC_CELLCHAT_TARGET_GROUPS", unset = "")

# Retrieve the LIANA method(s) to be used.
liana_method_str <- Sys.getenv("IMPACT_SC_LIANA_METHOD", unset = "logfc")


# --- 3. Validate Inputs ---
# ---------------------------
# Check if the required Seurat object exists.
if (!file.exists(seurat_object_path)) {
  stop(paste("FATAL: Input Seurat object not found at:", seurat_object_path))
}

# Check if source and target groups are provided.
if (nchar(source_groups_str) == 0 || nchar(target_groups_str) == 0) {
    stop("FATAL: Source and/or target groups were not provided. Please specify them in the setup script.")
}

# --- 4. Load Data ---
# ---------------------
# Read the annotated Seurat object.
message("Loading annotated data from Module 3...")
testdata <- readRDS(seurat_object_path)
message("Data loaded successfully.")

# Set the default identity for LIANA to use.
# It is assumed that the cell clusters/types are in the 'seurat_clusters' column.
# Change this if your cluster information is stored elsewhere.
Idents(testdata) <- "seurat_clusters"


# --- 5. Run LIANA Analysis ---
# ------------------------------
message(paste("Running LIANA workflow with method(s):", liana_method_str))
# Split the comma-separated string of methods into a character vector
liana_methods <- unlist(strsplit(liana_method_str, ","))

# liana_wrap is a wrapper that runs the specified methods.
liana_results <- liana_wrap(testdata, method = liana_methods)

# Conditionally aggregate results ONLY if more than one method was used.
if (length(liana_methods) > 1) {
  message("Multiple methods detected. Aggregating results...")
  liana_results <- liana_results %>%
    liana_aggregate()
  message("LIANA aggregation complete.")
} else {
  message("Single method detected. Skipping aggregation.")
  # If only one method is used, liana_aggregate fails.
  # We manually create the 'aggregate_rank' column for compatibility with plotting functions.
  # The ranking is based on the score of the single method used (e.g., logfc_comb, sca_score).
  # We find the score column and rank by it.
  score_col <- names(liana_results)[grepl("_score$|logfc_comb", names(liana_results))]
  if(length(score_col) > 0){
      # Use the first detected score column for ranking
      score_col_symbol <- rlang::sym(score_col[1])
      liana_results <- liana_results %>%
        arrange(!!score_col_symbol) %>%
        mutate(aggregate_rank = dense_rank(!!score_col_symbol))
  } else {
      # Fallback if no standard score column is found
      warning("Could not find a standard score column to create aggregate_rank. Downstream plots might fail.")
      liana_results$aggregate_rank <- 1 # Assign a dummy rank
  }
}

# Save the results for potential future use, directly in the main output directory.
saveRDS(liana_results, file.path(output_dir, "04h_liana_results.RDS"))

# --- 6. Process and Visualize Results ---
# -----------------------------------------
message("Generating visualizations...")

# Convert comma-separated string inputs to character vectors
source_groups_vec <- unlist(strsplit(source_groups_str, ","))
target_groups_vec <- unlist(strsplit(target_groups_str, ","))

# Generate and save the Dot Plot
# This plot shows the top N interactions, ranked by aggregate score.
message("Generating LIANA dotplot...")
tryCatch({
    p_dotplot <- liana_dotplot(liana_results,
                               source_groups = source_groups_vec,
                               target_groups = target_groups_vec,
                               ntop = 10) # Show top 10 interactions

    # Save the dot plot using the png() function for consistency
    png(file.path(output_dir, "04h_liana_dotplot.png"), width = 10, height = 8, units = 'in', res = 300)
    print(p_dotplot) # ggplot objects must be explicitly printed when using png()
    dev.off()
    
    message("LIANA dotplot saved.")
}, error = function(e) {
    message(paste("Warning: Could not generate LIANA dot plot. Error:", e$message))
})


# Filter results for downstream plots (e.g., keep top 1% of interactions)
liana_trunc <- liana_results %>%
  filter(aggregate_rank <= 0.01)

# Generate and save the Frequency Heatmap using liana::heat_freq
message("Generating frequency heatmap with liana::heat_freq...")
tryCatch({
    # The heat_freq function returns a ggplot object
    p_heatmap <- heat_freq(liana_trunc)

    # Save the heatmap plot to a PNG file
    png(file.path(output_dir, "04h_frequency_heatmap.png"), width = 8, height = 7, units = 'in', res = 300)
    print(p_heatmap) # ggplot objects must be explicitly printed
    dev.off()
    
    message("Frequency heatmap saved.")
}, error = function(e) {
    message(paste("Warning: Could not generate frequency heatmap. Error:", e$message))
})


# Generate and save the Chord Diagram using liana::chord_freq
# This plot visualizes the flow and strength of interactions between specific source and target groups.
message("Generating chord diagram...")
tryCatch({
    # Note: Chord diagrams can become cluttered. The filtering above helps.
    
    # Open the PNG device
    png(file.path(output_dir, "04h_interaction_chord_diagram.png"), width = 10, height = 10, units = 'in', res = 300)
    # The chord_freq function prints directly to the graphics device
    chord_freq(liana_trunc, source_groups = source_groups_vec, target_groups = target_groups_vec)
    # Close the device
    dev.off()
    
    message("Chord diagram saved.")
}, error = function(e) {
    # In case of error, ensure the graphics device is closed
    if(dev.cur() != 1) dev.off()
    message(paste("Warning: Could not generate chord diagram. Error:", e$message))
})

message("Module 04h (Cell Chat) processing complete.")
