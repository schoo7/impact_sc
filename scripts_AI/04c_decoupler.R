# IMPACT-sc Script: 04c_decoupler_analysis.R
# Purpose: Add Pathway and TF Activity Scores with DecoupleR.
# Version: Network download via decoupleR, ggplot2 visualization, explicit BiocParallel SerialParam

# --- Git Bash Troubleshooting Notes ---
# If Git Bash hangs or has issues with character encoding/locale, consider the following:
# 1. Set Git Bash locale environment variables (in ~/.bashrc or ~/.bash_profile):
#    export LANG=en_US.UTF-8
#    export LC_ALL=en_US.UTF-8
# 2. Update Git for Windows to the latest version.
# 3. Use Windows Terminal as a host for Git Bash for potentially better Unicode support.
# --- End Git Bash Troubleshooting Notes ---

# --- Libraries ---
library(AnnotationDbi)
library(decoupleR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(patchwork)
library(Seurat)
library(stringr)
library(tibble)
library(tidyr)
library(viridis)
library(reshape2) # Added for melt
library(BiocParallel) # Explicitly load BiocParallel to set params


# --- IMPACT-sc Script Parameters (READ FROM ENVIRONMENT VARIABLES) ---
species <- Sys.getenv("IMPACT_SC_SPECIES", "human")
set.seed(1234) # Setting seed for reproducibility

# Base output path for results
base_output_path <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", "../output/")
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
  message(paste("Created output directory:", base_output_path))
} else {
  message(paste("Output directory set to:", base_output_path))
}

# Load species-specific annotation database
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    library(org.Hs.eg.db)
    org_db_object <- org.Hs.eg.db
    message("Loaded org.Hs.eg.db for human.")
  } else {
    stop("org.Hs.eg.db not installed. Please install it for human species analysis.")
  }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    org_db_object <- org.Mm.eg.db
    message("Loaded org.Mm.eg.db for mouse.")
  } else {
    stop("org.Mm.eg.db not installed. Please install it for mouse species analysis.")
  }
} else {
  stop(paste("Species '", species, "' is not supported. Please use 'human' or 'mouse'.", sep=""))
}

options(future.globals.maxSize = 4 * 1024^3) # 4GB
# --- End Script Parameters ---

#### Module 4.3: Add Pathway and TF Activity Scores with DecoupleR ####
message("Starting Module 4.3: DecoupleR Analysis (TF and Pathway Activities)")

# Define path for the input Seurat object from the previous module
obj_prev_module_path <- file.path(base_output_path, "04b_data_after_diffex.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
  message(paste("Input '04b_data_after_diffex.RDS' not found. Attempting to load from fallback:", obj_prev_module_path))
}

if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required input Seurat object not found at either primary or fallback path:", obj_prev_module_path))
}

data <- tryCatch({
  readRDS(obj_prev_module_path)
}, error = function(e) {
  message(paste("Error loading input RDS file '", obj_prev_module_path, "': ", e$message, sep=""))
  NULL
})

if(is.null(data)) {
  stop("Failed to load Seurat object for Module 4.3. Halting execution.")
}
message(paste("Successfully loaded Seurat object from:", obj_prev_module_path))

data[["RNA"]] <- split(data[["RNA"]], f = data$cell_type)  # Split RNA assay by sample
data<- JoinLayers(data) 
data$cell_type <- as.factor(data$cell_type)
DefaultAssay(data) <- "RNA"

if (!("data" %in% Layers(data, assay="RNA"))) {
    message("RNA assay 'data' layer is missing. Normalizing data using NormalizeData()...")
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
} else {
    message("RNA assay 'data' layer found. Using existing log-normalized data.")
}

mat_norm_decoupler <- GetAssayData(data, assay = "RNA", layer = "data")
if(is.null(mat_norm_decoupler) || ncol(mat_norm_decoupler) == 0 || nrow(mat_norm_decoupler) == 0){
    stop("Normalized data matrix for DecoupleR is empty or could not be retrieved. Halting.")
}
message(paste("Using species '", species, "' for DecoupleR network selection and analysis.", sep=""))

# --- TF activity (CollecTRI downloaded via decoupleR) ---
net_collectri <- NULL
message(paste("Attempting to download CollecTRI network for species:", species))
net_collectri <- tryCatch({
  decoupleR::get_collectri(organism = species, split_complexes = FALSE)
}, error = function(e) {
  message(paste("Error downloading CollecTRI network for '", species, "': ", e$message, sep=""))
  NULL
})

if (!is.null(net_collectri)) {
  message("CollecTRI network downloaded successfully.")
}

tf_acts_ulm <- NULL # Initialize before the if block
if(!is.null(net_collectri)){
  message("Running ULM (Univariate Linear Model) for TF activities using CollecTRI network.")
  message("Setting BiocParallel to SerialParam to avoid potential parallel backend issues.")
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

  current_mat_tf <- if (is(mat_norm_decoupler, "dgCMatrix")) mat_norm_decoupler else as.matrix(mat_norm_decoupler)
  
  genes_in_matrix_collectri <- rownames(current_mat_tf)
  genes_in_network_collectri <- unique(net_collectri$target)
  common_genes_collectri <- intersect(genes_in_matrix_collectri, genes_in_network_collectri)
  message(paste("Number of common genes between matrix and CollecTRI network:", length(common_genes_collectri)))
  if (length(common_genes_collectri) == 0) {
    message("WARNING: No common genes found between your expression data and the CollecTRI network. This will likely lead to issues or empty results.")
  }
  
  tf_acts_ulm <- decoupleR::run_ulm(current_mat_tf, net = net_collectri,
                     .source = 'source', .target = 'target', .mor = 'mor',
                     minsize = 5)

  message("run_ulm for TF activities finished.")


  if(!is.null(tf_acts_ulm) && nrow(tf_acts_ulm) > 0) {
    tf_scores_wide <- tf_acts_ulm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
      tibble::column_to_rownames('source')

    data[['tfs_ulm']] <- CreateAssayObject(counts = as.matrix(tf_scores_wide))
    DefaultAssay(data) <- "tfs_ulm"
    data <- ScaleData(data, assay = "tfs_ulm", verbose = FALSE)
    message("TF activities (ULM) added as 'tfs_ulm' assay and scaled.")

    n_top_tfs_decoupler <- 25
    avg_tf_activity_raw <- AverageExpression(data, assays = "tfs_ulm", features = rownames(data[["tfs_ulm"]]), group.by = "cell_type", layer = "counts", verbose = FALSE)$tfs_ulm

    if (nrow(avg_tf_activity_raw) > 0) {
        tf_variances_decoupler <- apply(avg_tf_activity_raw, 1, var, na.rm = TRUE)
        tf_variances_decoupler <- tf_variances_decoupler[!is.na(tf_variances_decoupler)]
        if (length(tf_variances_decoupler) > 0) {
            top_variable_tfs_decoupler <- names(sort(tf_variances_decoupler, decreasing = TRUE)[1:min(n_top_tfs_decoupler, length(tf_variances_decoupler))])
            
            # --- FIX: Added check to ensure top_variable_tfs_decoupler is not empty ---
            if (!is.null(top_variable_tfs_decoupler) && length(top_variable_tfs_decoupler) > 0) {
                avg_tf_activity_top_decoupler <- avg_tf_activity_raw[top_variable_tfs_decoupler, , drop=FALSE]

                if(nrow(avg_tf_activity_top_decoupler) > 0 && ncol(avg_tf_activity_top_decoupler) > 0) {
                    scaled_avg_tf_activity <- t(scale(t(avg_tf_activity_top_decoupler)))
                    scaled_avg_tf_activity[is.na(scaled_avg_tf_activity) | is.nan(scaled_avg_tf_activity) | is.infinite(scaled_avg_tf_activity)] <- 0
                    melted_tf_data <- reshape2::melt(scaled_avg_tf_activity, varnames = c("Gene", "CellType"), value.name = "Value")

                    tf_heatmap_plot <- ggplot(melted_tf_data, aes(x = Gene, y = CellType)) +
                      geom_tile(aes(fill = Value), color = "grey50") +
                      scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Activity") +
                      theme_minimal() +
                      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                      ggtitle("Mean TF Activity per Cell Type (Top Variable TFs)")

                    tf_heatmap_path_decoupler_gg <- file.path(base_output_path, "04c_tf_activity_heatmap_decoupler_ggplot.png")
                    ggsave(tf_heatmap_path_decoupler_gg, plot = tf_heatmap_plot, width = 10, height = max(6, 0.4*nlevels(as.factor(data$cell_type))), limitsize = FALSE)
                    message(paste("DecoupleR TF activity ggplot heatmap saved to:", tf_heatmap_path_decoupler_gg))
                } else { message("Not enough variable TFs or cell types to plot TF activity ggplot heatmap after filtering.")}
            } else {
                message("No top variable TFs could be identified after filtering. Skipping TF activity heatmap.")
            }
        } else { message("No TFs with variance found. Skipping TF activity heatmap.")}
    } else { message("Average TF activity matrix is empty. Skipping TF activity heatmap.")}
  } else {
      message("TF activity calculation (run_ulm) did not produce expected results. Skipping TF analysis downstream.")
  }
} else { message("Skipping TF activity analysis as CollecTRI network could not be downloaded.")}

# --- Pathway activity (PROGENy downloaded via decoupleR) ---
net_progeny <- NULL
message(paste("Attempting to download PROGENy network for species:", species))
net_progeny <- tryCatch({
  # PROGENy is available for human and mouse in decoupleR
  decoupleR::get_progeny(organism = species, top = 500)
}, error = function(e) {
  message(paste("Error downloading PROGENy network for '", species, "': ", e$message, sep=""))
  NULL
})

if (!is.null(net_progeny)) {
  message("PROGENy network downloaded successfully.")
}

pathway_acts_mlm <- NULL # Initialize before the if block
if(!is.null(net_progeny)){
  message("Running MLM (Multivariate Linear Model) for PROGENy pathway activities.")
  message("Setting BiocParallel to SerialParam to avoid potential parallel backend issues.")
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

  current_mat_pathway <- if (is(mat_norm_decoupler, "dgCMatrix")) mat_norm_decoupler else as.matrix(mat_norm_decoupler)

  genes_in_matrix_progeny <- rownames(current_mat_pathway)
  genes_in_network_progeny <- unique(net_progeny$target)
  common_genes_progeny <- intersect(genes_in_matrix_progeny, genes_in_network_progeny)
  message(paste("Number of common genes between matrix and PROGENy network:", length(common_genes_progeny)))
  if (length(common_genes_progeny) == 0 && nrow(net_progeny) > 0) {
    message("WARNING: No common genes found between your expression data and the PROGENy network. This will likely lead to issues or empty results.")
  }

  pathway_acts_mlm <- decoupleR::run_mlm(current_mat_pathway, net = net_progeny,
                     .source = 'source', .target = 'target', .mor = 'weight',
                     minsize = 5)

  message("run_mlm for pathway activities finished.")

  if(!is.null(pathway_acts_mlm) && nrow(pathway_acts_mlm) > 0) {
    pathway_scores_wide <- pathway_acts_mlm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
      tibble::column_to_rownames('source')

    data[['pathways_mlm']] <- CreateAssayObject(counts = as.matrix(pathway_scores_wide))
    DefaultAssay(data) <- "pathways_mlm"
    data <- ScaleData(data, assay = "pathways_mlm", verbose = FALSE)
    message("PROGENy pathway activities (MLM) added as 'pathways_mlm' assay and scaled.")

    avg_pathway_activity_raw <- AverageExpression(data, assays = "pathways_mlm", features = rownames(data[["pathways_mlm"]]), group.by = "cell_type", layer="counts", verbose=FALSE)$pathways_mlm

    if(nrow(avg_pathway_activity_raw) > 0 && ncol(avg_pathway_activity_raw) > 0) {
        scaled_avg_pathway_activity <- t(scale(t(avg_pathway_activity_raw)))
        scaled_avg_pathway_activity[is.na(scaled_avg_pathway_activity) | is.nan(scaled_avg_pathway_activity) | is.infinite(scaled_avg_pathway_activity)] <- 0
        melted_pathway_data <- reshape2::melt(scaled_avg_pathway_activity, varnames = c("Pathway", "CellType"), value.name = "Activity")

        pathway_heatmap_plot <- ggplot(melted_pathway_data, aes(x = Pathway, y = CellType, fill = Activity)) +
          geom_tile(color = "grey50") +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Activity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Mean Pathway Activity per Cell Type (PROGENy)")

        pathway_heatmap_path_decoupler_gg <- file.path(base_output_path, "04c_pathway_activity_heatmap_decoupler_ggplot.png")
        ggsave(pathway_heatmap_path_decoupler_gg, plot = pathway_heatmap_plot, width = 10, height = max(6, 0.4*nlevels(as.factor(data$cell_type))), limitsize=FALSE)
        message(paste("DecoupleR Pathway (PROGENy) activity ggplot heatmap saved to:", pathway_heatmap_path_decoupler_gg))
    } else { message("Not enough data (pathways or cell types) to plot PROGENy ggplot activity heatmap after processing.")}
  } else {
      message("Pathway activity calculation (run_mlm) did not produce expected results. Skipping PROGENy analysis downstream.")
  }
} else { message("Skipping PROGENy pathway activity analysis as network could not be downloaded.")}

DefaultAssay(data) <- "RNA"
message("Default assay reset to 'RNA'.")

intermediate_rds_path_04c <- file.path(base_output_path, "04c_data_after_decoupler.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04c)
  message(paste("Data object after DecoupleR analysis successfully saved to:", intermediate_rds_path_04c))
}, error = function(e) {
  warning(paste("Error saving Seurat data object after Module 4.3 to '", intermediate_rds_path_04c, "': ", e$message, sep=""))
})

message("Finished Module 04c: DecoupleR Analysis.")
