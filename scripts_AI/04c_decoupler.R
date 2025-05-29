# IMPACT-sc Script: 04c_decoupler_analysis.R
# Purpose: Add Pathway and TF Activity Scores with DecoupleR.
# Version: CSV input for networks, ggplot2 visualization, explicit BiocParallel SerialParam, no tryCatch on run_ulm/mlm

# --- Git Bash Troubleshooting Notes ---
# If Git Bash hangs or has issues with character encoding/locale, consider the following:
# 1. Set Git Bash locale environment variables (in ~/.bashrc or ~/.bash_profile):
#    export LANG=en_US.UTF-8
#    export LC_ALL=en_US.UTF-8
#    (or your preferred locale with .UTF-8, e.g., zh_CN.UTF-8)
#    Source the file (e.g., source ~/.bashrc) or restart Git Bash.
#
# 2. Configure Git's output encoding:
#    git config --global core.quotepath off
#
# 3. Set LESSCHARSET for the 'less' pager (in ~/.bashrc):
#    export LESSCHARSET=utf-8
#
# 4. Modify Git Bash terminal settings:
#    - Right-click Git Bash title bar -> Options...
#    - Text -> Character set: UTF-8
#    - Apply/Save.
#
# 5. Update Git for Windows to the latest version.
#
# 6. Check Windows system-level region/language settings:
#    - Search for "Region settings" or "Language settings".
#    - Administrative language settings -> "Language for non-Unicode programs".
#    - Consider "Beta: Use Unicode UTF-8 for worldwide language support" (use with caution).
#
# 7. Use Windows Terminal as a host for Git Bash for potentially better Unicode support.
# --- End Git Bash Troubleshooting Notes ---

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

# Get network file paths from environment variables
collectri_csv_path <- Sys.getenv("IMPACT_SC_COLLECTRI_CSV_PATH", "")
progeny_csv_path <- Sys.getenv("IMPACT_SC_PROGENY_CSV_PATH", "")


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

data <- tryCatch({ # Keep tryCatch for initial data loading
  readRDS(obj_prev_module_path)
}, error = function(e) {
  message(paste("Error loading input RDS file '", obj_prev_module_path, "': ", e$message, sep=""))
  NULL
})

if(is.null(data)) {
  stop("Failed to load Seurat object for Module 4.3. Halting execution.")
}
message(paste("Successfully loaded Seurat object from:", obj_prev_module_path))

DefaultAssay(data) <- "RNA"

if (!("data" %in% slotNames(data@assays$RNA))) {
    message("RNA assay 'data' slot is missing. Normalizing data using NormalizeData()...")
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
} else {
    message("RNA assay 'data' slot found. Using existing log-normalized data.")
}

mat_norm_decoupler <- GetAssayData(data, assay = "RNA", layer = "data")
if(is.null(mat_norm_decoupler) || ncol(mat_norm_decoupler) == 0 || nrow(mat_norm_decoupler) == 0){
    stop("Normalized data matrix for DecoupleR is empty or could not be retrieved. Halting.")
}
message(paste("Using species '", species, "' for DecoupleR network selection (if applicable) and analysis.", sep=""))

# --- TF activity (CollecTRI from CSV) ---
net_collectri <- NULL
if (collectri_csv_path != "" && file.exists(collectri_csv_path)) {
  message(paste("Attempting to load CollecTRI network from user-provided CSV:", collectri_csv_path))
  net_collectri <- tryCatch({ # Keep tryCatch for file reading
    read.csv(collectri_csv_path, stringsAsFactors = FALSE, header = TRUE)
  }, error = function(e) {
    message(paste("Error loading CollecTRI CSV file '", collectri_csv_path, "': ", e$message, sep=""))
    NULL
  })

  if (!is.null(net_collectri)) {
    required_collectri_cols <- c("source", "target", "mor")
    if (!all(required_collectri_cols %in% colnames(net_collectri))) {
      message(paste("Error: CollecTRI CSV ('", collectri_csv_path, "') must contain the columns: 'source', 'target', and 'mor'. Found columns: ", paste(colnames(net_collectri), collapse=", "), sep=""))
      net_collectri <- NULL
    } else {
      message("CollecTRI CSV loaded and columns validated successfully.")
    }
  }
} else {
  if (collectri_csv_path == "") {
    message("IMPACT_SC_COLLECTRI_CSV_PATH environment variable not set. Skipping TF activity analysis.")
  } else {
    message(paste("CollecTRI CSV file specified by IMPACT_SC_COLLECTRI_CSV_PATH ('", collectri_csv_path, "') does not exist. Skipping TF activity analysis.", sep=""))
  }
}

tf_acts_ulm <- NULL # Initialize before the if block
if(!is.null(net_collectri)){
  message("Running ULM (Univariate Linear Model) for TF activities using CollecTRI network.")
  message("Setting BiocParallel to SerialParam to avoid potential parallel backend issues.")
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

  current_mat_tf <- if (is(mat_norm_decoupler, "dgCMatrix")) mat_norm_decoupler else as.matrix(mat_norm_decoupler)
  message(paste("Dimensions of matrix for run_ulm: Rows =", nrow(current_mat_tf), "Cols =", ncol(current_mat_tf)))
  message(paste("Object size of current_mat_tf:", format(object.size(current_mat_tf), units = "Mb")))
  message(paste("Number of rows in net_collectri:", nrow(net_collectri)))
  message("First few rows of net_collectri:")
  print(head(net_collectri))
  
  genes_in_matrix_collectri <- rownames(current_mat_tf)
  genes_in_network_collectri <- unique(net_collectri$target)
  common_genes_collectri <- intersect(genes_in_matrix_collectri, genes_in_network_collectri)
  message(paste("Number of common genes between matrix and CollecTRI network:", length(common_genes_collectri)))
  if (length(common_genes_collectri) == 0) {
    message("WARNING: No common genes found between your expression data and the CollecTRI network. This will likely lead to issues or empty results.")
  }
  message(paste("Object size of net_collectri:", format(object.size(net_collectri), units = "Mb")))


  # Directly run without tryCatch for run_ulm
  tf_acts_ulm <- decoupleR::run_ulm(current_mat_tf, net = net_collectri,
                     .source = 'source', .target = 'target', .mor = 'mor',
                     minsize = 5)

  message("run_ulm finished.") # Add this message to see if it completes


  if(!is.null(tf_acts_ulm) && nrow(tf_acts_ulm) > 0) {
    tf_scores_wide <- tf_acts_ulm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
      tibble::column_to_rownames('source')

    data[['tfs_ulm']] <- CreateAssayObject(counts = as.matrix(tf_scores_wide))
    DefaultAssay(data) <- "tfs_ulm"
    data <- ScaleData(data, assay = "tfs_ulm", verbose = FALSE)
    message("TF activities (ULM) added as 'tfs_ulm' assay and scaled.")

    n_top_tfs_decoupler <- 25
    avg_tf_activity_raw <- AverageExpression(data, assays = "tfs_ulm", features = rownames(data[["tfs_ulm"]]), group.by = "cell_type", slot = "counts", verbose = FALSE)$tfs_ulm

    if (nrow(avg_tf_activity_raw) > 0) {
        tf_variances_decoupler <- apply(avg_tf_activity_raw, 1, var, na.rm = TRUE)
        tf_variances_decoupler <- tf_variances_decoupler[!is.na(tf_variances_decoupler)]
        if (length(tf_variances_decoupler) > 0) {
            top_variable_tfs_decoupler <- names(sort(tf_variances_decoupler, decreasing = TRUE)[1:min(n_top_tfs_decoupler, length(tf_variances_decoupler))])
            avg_tf_activity_top_decoupler <- avg_tf_activity_raw[top_variable_tfs_decoupler, , drop=FALSE]

            if(nrow(avg_tf_activity_top_decoupler) > 0 && ncol(avg_tf_activity_top_decoupler) > 0) {
                scaled_avg_tf_activity <- t(scale(t(avg_tf_activity_top_decoupler)))
                scaled_avg_tf_activity[is.na(scaled_avg_tf_activity) | is.nan(scaled_avg_tf_activity) | is.infinite(scaled_avg_tf_activity)] <- 0
                melted_tf_data <- reshape2::melt(scaled_avg_tf_activity, varnames = c("Gene", "CellType"), value.name = "Value")
                melted_tf_data$CellType <- factor(melted_tf_data$CellType, levels = rev(unique(as.character(melted_tf_data$CellType))))
                melted_tf_data$Gene <- factor(melted_tf_data$Gene, levels = unique(as.character(melted_tf_data$Gene)))
                max_tile_size <- 1

                tf_heatmap_plot <- ggplot(melted_tf_data, aes(x = Gene, y = CellType)) +
                  geom_tile(aes(fill = Value), color = "grey50", linewidth = 0.5,
                            width = pmin(abs(melted_tf_data$Value) / 2, max_tile_size),
                            height = pmin(abs(melted_tf_data$Value) / 2, max_tile_size)) +
                  scale_fill_gradient2(low = "#1E90FF", mid = "white", high = "#FF69B4", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Activity") +
                  theme_minimal(base_size = 12) +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "#374151", size=10),
                        axis.text.y = element_text(face = "bold", color = "#374151", size=10),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        legend.title = element_text(face = "bold", size=11),
                        legend.text = element_text(color = "#374151", size=10),
                        plot.background = element_rect(fill = "white", color = NA),
                        plot.margin = margin(15, 15, 15, 15),
                        plot.title = element_text(hjust = 0.5, face = "bold", size=14)) +
                  geom_tile(aes(x = Gene, y = CellType), color = "grey", linewidth = 0.2, fill = NA, width = 1, height = 1) +
                  coord_fixed() +
                  ggtitle("Mean TF Activity per Cell Type (Top Variable TFs)")

                tf_heatmap_path_decoupler_gg <- file.path(base_output_path, "04c_tf_activity_heatmap_decoupler_ggplot.png")
                ggsave(tf_heatmap_path_decoupler_gg, plot = tf_heatmap_plot, width = 10, height = max(6, 0.4*length(unique(melted_tf_data$CellType)) + 0.3*length(unique(melted_tf_data$Gene))), units = "in", dpi = 300, limitsize = FALSE)
                message(paste("DecoupleR TF activity ggplot heatmap saved to:", tf_heatmap_path_decoupler_gg))
            } else { message("Not enough variable TFs or cell types to plot TF activity ggplot heatmap after filtering.")}
        } else { message("No TFs with variance found. Skipping TF activity heatmap.")}
    } else { message("Average TF activity matrix is empty. Skipping TF activity heatmap.")}
  } else {
      if (is.null(tf_acts_ulm)) {
          message("TF activity calculation (run_ulm) resulted in NULL. Skipping TF analysis downstream.")
      } else if (nrow(tf_acts_ulm) == 0) {
          message("TF activity calculation (run_ulm) returned 0 rows. Skipping TF analysis downstream.")
      } else {
          message("TF activity calculation (run_ulm) did not produce expected results. Skipping TF analysis downstream.")
      }
  }
} else { message("Skipping TF activity analysis as CollecTRI network was not loaded successfully.")}

# --- Pathway activity (PROGENy from CSV, only for human) ---
net_progeny <- NULL
if (species == "human") {
  if (progeny_csv_path != "" && file.exists(progeny_csv_path)) {
    message(paste("Attempting to load PROGENy network from user-provided CSV for human:", progeny_csv_path))
    net_progeny <- tryCatch({ # Keep tryCatch for file reading
      read.csv(progeny_csv_path, stringsAsFactors = FALSE, header = TRUE)
    }, error = function(e) {
      message(paste("Error loading PROGENy CSV file '", progeny_csv_path, "': ", e$message, sep=""))
      NULL
    })

    if (!is.null(net_progeny)) {
      required_progeny_cols <- c("source", "target", "weight")
      if (!all(required_progeny_cols %in% colnames(net_progeny))) {
        message(paste("Error: PROGENy CSV ('", progeny_csv_path, "') must contain the columns: 'source', 'target', and 'weight'. Found columns: ", paste(colnames(net_progeny), collapse=", "), sep=""))
        net_progeny <- NULL
      } else {
        message("PROGENy CSV loaded and columns validated successfully for human.")
      }
    }
  } else {
    if (progeny_csv_path == "") {
      message("IMPACT_SC_PROGENY_CSV_PATH environment variable not set for human. Skipping PROGENy pathway analysis.")
    } else {
      message(paste("PROGENy CSV file specified by IMPACT_SC_PROGENY_CSV_PATH ('", progeny_csv_path, "') does not exist for human. Skipping PROGENy pathway analysis.", sep=""))
    }
  }
} else {
  message(paste("PROGENy pathway analysis is designed for human species. Skipping for '", species, "'.", sep=""))
}

pathway_acts_mlm <- NULL # Initialize before the if block
if(!is.null(net_progeny) && species == "human"){
  message("Running MLM (Multivariate Linear Model) for PROGENy pathway activities.")
  message("Setting BiocParallel to SerialParam to avoid potential parallel backend issues.")
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

  current_mat_pathway <- if (is(mat_norm_decoupler, "dgCMatrix")) mat_norm_decoupler else as.matrix(mat_norm_decoupler)
  message(paste("Dimensions of matrix for run_mlm: Rows =", nrow(current_mat_pathway), "Cols =", ncol(current_mat_pathway)))
  message(paste("Object size of current_mat_pathway:", format(object.size(current_mat_pathway), units = "Mb")))
  message(paste("Number of rows in net_progeny:", nrow(net_progeny)))
  message("First few rows of net_progeny:")
  print(head(net_progeny))

  genes_in_matrix_progeny <- rownames(current_mat_pathway)
  genes_in_network_progeny <- unique(net_progeny$target) # Assuming 'target' contains gene symbols
  common_genes_progeny <- intersect(genes_in_matrix_progeny, genes_in_network_progeny)
  message(paste("Number of common genes between matrix and PROGENy network:", length(common_genes_progeny)))
  if (length(common_genes_progeny) == 0 && nrow(net_progeny) > 0) { # only warn if network is not empty
    message("WARNING: No common genes found between your expression data and the PROGENy network for human. This will likely lead to issues or empty results.")
  }
  message(paste("Object size of net_progeny:", format(object.size(net_progeny), units = "Mb")))

  # Directly run without tryCatch for run_mlm
  pathway_acts_mlm <- decoupleR::run_mlm(current_mat_pathway, net = net_progeny,
                     .source = 'source', .target = 'target', .mor = 'weight',
                     minsize = 5)

  message("run_mlm finished.") # Add this message

  if(!is.null(pathway_acts_mlm) && nrow(pathway_acts_mlm) > 0) {
    pathway_scores_wide <- pathway_acts_mlm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
      tibble::column_to_rownames('source')

    data[['pathways_mlm']] <- CreateAssayObject(counts = as.matrix(pathway_scores_wide))
    DefaultAssay(data) <- "pathways_mlm"
    data <- ScaleData(data, assay = "pathways_mlm", verbose = FALSE)
    message("PROGENy pathway activities (MLM) added as 'pathways_mlm' assay and scaled.")

    avg_pathway_activity_raw <- AverageExpression(data, assays = "pathways_mlm", features = rownames(data[["pathways_mlm"]]), group.by = "cell_type", slot="counts", verbose=FALSE)$pathways_mlm

    if(nrow(avg_pathway_activity_raw) > 0 && ncol(avg_pathway_activity_raw) > 0) {
        scaled_avg_pathway_activity <- t(scale(t(avg_pathway_activity_raw)))
        scaled_avg_pathway_activity[is.na(scaled_avg_pathway_activity) | is.nan(scaled_avg_pathway_activity) | is.infinite(scaled_avg_pathway_activity)] <- 0
        melted_pathway_data <- reshape2::melt(scaled_avg_pathway_activity, varnames = c("Gene", "CellType"), value.name = "Value")
        melted_pathway_data$CellType <- factor(melted_pathway_data$CellType, levels = rev(unique(as.character(melted_pathway_data$CellType))))
        melted_pathway_data$Gene <- factor(melted_pathway_data$Gene, levels = unique(as.character(melted_pathway_data$Gene)))
        max_tile_size <- 1

        pathway_heatmap_plot <- ggplot(melted_pathway_data, aes(x = Gene, y = CellType)) +
          geom_tile(aes(fill = Value), color = "grey50", linewidth = 0.5,
                    width = pmin(abs(melted_pathway_data$Value) / 2, max_tile_size),
                    height = pmin(abs(melted_pathway_data$Value) / 2, max_tile_size)) +
          scale_fill_gradient2(low = "#1E90FF", mid = "white", high = "#FF69B4", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Activity") +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "#374151", size=10),
                axis.text.y = element_text(face = "bold", color = "#374151", size=10),
                panel.grid = element_blank(),
                axis.title = element_blank(),
                legend.title = element_text(face = "bold", size=11),
                legend.text = element_text(color = "#374151", size=10),
                plot.background = element_rect(fill = "white", color = NA),
                plot.margin = margin(15, 15, 15, 15),
                plot.title = element_text(hjust = 0.5, face = "bold", size=14)) +
          geom_tile(aes(x = Gene, y = CellType), color = "grey", linewidth = 0.2, fill = NA, width = 1, height = 1) +
          coord_fixed() +
          ggtitle("Mean Pathway Activity per Cell Type (PROGENy)")

        pathway_heatmap_path_decoupler_gg <- file.path(base_output_path, "04c_pathway_activity_heatmap_decoupler_ggplot.png")
        ggsave(pathway_heatmap_path_decoupler_gg, plot = pathway_heatmap_plot, width = 10, height = max(6, 0.4*length(unique(melted_pathway_data$CellType)) + 0.3*length(unique(melted_pathway_data$Gene))), units = "in", dpi = 300, limitsize = FALSE)
        message(paste("DecoupleR Pathway (PROGENy) activity ggplot heatmap saved to:", pathway_heatmap_path_decoupler_gg))
    } else { message("Not enough data (pathways or cell types) to plot PROGENy ggplot activity heatmap after processing.")}
  } else {
      if (is.null(pathway_acts_mlm)) {
          message("Pathway activity calculation (run_mlm) resulted in NULL. Skipping PROGENy analysis downstream.")
      } else if (nrow(pathway_acts_mlm) == 0) {
          message("Pathway activity calculation (run_mlm) returned 0 rows. Skipping PROGENy analysis downstream.")
      } else {
          message("Pathway activity calculation (run_mlm) did not produce expected results. Skipping PROGENy analysis downstream.")
      }
  }
} else { message("Skipping PROGENy pathway activity analysis as network was not loaded successfully or species is not human.")}

DefaultAssay(data) <- "RNA"
message("Default assay reset to 'RNA'.")

intermediate_rds_path_04c <- file.path(base_output_path, "04c_data_after_decoupler.RDS")
tryCatch({ # Keep tryCatch for saving data
  saveRDS(data, intermediate_rds_path_04c)
  message(paste("Data object after DecoupleR analysis successfully saved to:", intermediate_rds_path_04c))
}, error = function(e) {
  warning(paste("Error saving Seurat data object after Module 4.3 to '", intermediate_rds_path_04c, "': ", e$message, sep=""))
})

message("Finished Module 4.3: DecoupleR Analysis.")