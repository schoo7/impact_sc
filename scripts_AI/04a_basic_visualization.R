# IMPACT-sc Script: 04a_basic_visualization.R
# Purpose: Perform basic visualizations of the annotated Seurat object.

# --- Libraries ---
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(patchwork)
library(SCpubr)
library(Seurat)
library(stringr)
library(tibble)
library(viridis)
library(jsonlite) # For parsing JSON string from env var
library(SeuratExtend) # Added for DimPlot2, FeaturePlot3, DotPlot2
library(scRNAtoolVis) # Added for ClusterDistrBar, CalcStats, Heatmap
library(future) # Explicitly load future

# --- Set Future Plan to Sequential to Aid Debugging Potential Hangs ---
# This is a common troubleshooting step for hangs that occur in one environment but not another,
# as it disables parallelism which can sometimes be problematic.
plan(sequential)
message("Future plan set to sequential.")

# --- IMPACT-sc Script Parameters (READ FROM ENVIRONMENT VARIABLES) ---
species <- Sys.getenv("IMPACT_SC_SPECIES", "human")
set.seed(1234)
base_output_path <- Sys.getenv("IMPACT_SC_OUTPUT_DIR", "../output/")
if (!dir.exists(base_output_path)) dir.create(base_output_path, recursive = TRUE)
message(paste("Output directory set to:", base_output_path))

# --- [MODIFICATION] ---
# Read the user's chosen reduction method from the environment variable.
user_reduction_method <- Sys.getenv("IMPACT_SC_REDUCTION_METHOD", "") # Default to empty string

# Visualization specific parameters from environment
# These will now strictly depend on user input via environment variables.
# No R-script level defaults for the gene lists themselves.
featureplot_genes_env_var <- Sys.getenv("IMPACT_SC_FEATUREPLOT_GENES", "") # Default to empty string if not set
dotplot_genes_json_env_var <- Sys.getenv("IMPACT_SC_DOTPLOT_GENES_JSON", "[]") # Default to empty JSON array string if not set

if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { library(org.Hs.eg.db); org_db_object <- org.Hs.eg.db } else { stop("org.Hs.eg.db not installed.") }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { library(org.Mm.eg.db); org_db_object <- org.Mm.eg.db } else { stop("org.Mm.eg.db not installed.") }
} else { stop("Species not supported.") }
options(future.globals.maxSize = 4 * 1024^3) # 4GB
# --- End Script Parameters ---

#### Module 4.1: Basic Visualization ####
message("Starting Module 4.1: Basic Visualization")

# Input: '03_module3_final_annotated.RDS'
obj_module3_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
if (!file.exists(obj_module3_path)) {
  stop(paste("Annotated Seurat object from Module 3 not found at:", obj_module3_path,
             "\nPlease run previous modules first."))
}
data <- tryCatch(readRDS(obj_module3_path), error=function(e){
  message(paste("Error loading 03_module3_final_annotated.RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.1.")
message("Loaded final annotated object for Module 4.1.")

data$cell_type <- as.factor(data$cell_type)
DefaultAssay(data) <- "RNA"
celltype_counts <- table(data$cell_type)
valid_celltypes <- names(celltype_counts[celltype_counts >= 20])
data <- subset(data, subset = cell_type %in% valid_celltypes)

# --- [MODIFICATION] ---
# Determine reduction to use for plotting.
# Priority: 1. User's choice, 2. Fallback list (umap_c2s, umap, harmony), 3. Default to pca.
available_reductions <- Reductions(data)
reduction_to_use <- NULL

# 1. Check if user's choice is valid and available.
if (user_reduction_method != "" && user_reduction_method %in% available_reductions) {
    reduction_to_use <- user_reduction_method
    message(paste("Using user-selected reduction '", reduction_to_use, "' for visualization.", sep=""))
} else {
    if (user_reduction_method != "") {
        message(paste("Warning: User-selected reduction '", user_reduction_method, "' is not available in the Seurat object. Available reductions: ", paste(available_reductions, collapse=", "), ". Will proceed with fallback logic.", sep=""))
    }
    
    # 2. If user's choice is not used, go to fallback priority list.
    fallback_priority <- c("umap_c2s", "umap", "harmony", "pca")
    for (r_name in fallback_priority) {
        if (r_name %in% available_reductions) {
            reduction_to_use <- r_name
            message(paste("Using fallback reduction '", reduction_to_use, "' for visualization.", sep=""))
            break
        }
    }
}

# 3. If after all checks, no reduction is selected (highly unlikely if 'pca' is always present), stop.
if (is.null(reduction_to_use)) {
    stop("Could not find any suitable reduction (umap_c2s, umap, harmony, or pca) in the Seurat object. Stopping.")
}

message(paste("Final reduction selected for plotting:", reduction_to_use))


# --- DimPlot ---
# SeuratExtend::DimPlot2
p_dimplot_seuratextend <- tryCatch({
    SeuratExtend::DimPlot2(data, reduction = reduction_to_use, group.by = "cell_type", theme = theme_umap_arrows()) # "cell_type" is the correct string parameter
}, error = function(e) {
    message(paste("Error generating SeuratExtend::DimPlot2:", e$message)); NULL
})
if (!is.null(p_dimplot_seuratextend)) {
  ggsave(file.path(base_output_path, paste0("04a_basicviz_dimplot_seuratextend_", reduction_to_use, ".png")), plot = p_dimplot_seuratextend, width = 16, height = 14)
  message(paste0("SeuratExtend::DimPlot2 saved using ", reduction_to_use, "."))
}

# SCpubr::do_DimPlot
p_dimplot_scpubr <- tryCatch({
    SCpubr::do_DimPlot(sample = data, reduction = reduction_to_use, group.by = "cell_type", plot.axes = TRUE, label = F) # "cell_type" is the correct string parameter
}, error = function(e) {
    message(paste("Error generating SCpubr::do_DimPlot:", e$message, "- No Seurat basic fallback will be used.")); NULL
})
if (!is.null(p_dimplot_scpubr)) {
  ggsave(file.path(base_output_path, paste0("04a_basicviz_dimplot_scpubr_", reduction_to_use, ".png")), plot = p_dimplot_scpubr, width = 16, height = 14)
  message(paste0("SCpubr::do_DimPlot saved using ", reduction_to_use, "."))
}


# --- Proportion Plot ---
# scRNAtoolVis::ClusterDistrBar
if ("orig.ident" %in% colnames(data@meta.data) && nlevels(as.factor(data$orig.ident)) > 0 && "cell_type" %in% colnames(data@meta.data)) { # "orig.ident", "cell_type" are correct strings
  p_prop_clusterdistrbar <- tryCatch({
    ClusterDistrBar(data$orig.ident, data$cell_type) # Parameters are passed as vectors directly, no quotes needed
  }, error = function(e) {
    message(paste("Error generating scRNAtoolVis::ClusterDistrBar for proportions:", e$message)); NULL
  })
  if(!is.null(p_prop_clusterdistrbar)){
    ggsave(file.path(base_output_path, "04a_basicviz_proportion_clusterdistrbar.png"), plot = p_prop_clusterdistrbar, width = 10, height = 7)
    message("scRNAtoolVis::ClusterDistrBar proportion plot saved.")
  }
} else {message("Skipping proportion plot (ClusterDistrBar) as 'orig.ident' or 'cell_type' is missing or 'orig.ident' has only one level.")}


# --- Feature Plots ---
genes_for_featureplot_present <- c() # Initialize as empty character vector
if (is.null(featureplot_genes_env_var) || featureplot_genes_env_var == "") {
    message("IMPACT_SC_FEATUREPLOT_GENES environment variable not set or empty. Skipping FeaturePlot based on user input.")
} else {
    genes_for_featureplot_input <- strsplit(featureplot_genes_env_var, ",")[[1]]
    # Remove any empty strings that might result from trailing commas or multiple commas
    genes_for_featureplot_input <- genes_for_featureplot_input[genes_for_featureplot_input != "" & !is.na(genes_for_featureplot_input)]

    if (length(genes_for_featureplot_input) == 0) {
        message("No valid genes provided in IMPACT_SC_FEATUREPLOT_GENES after parsing. Skipping FeaturePlot.")
    } else {
        genes_for_featureplot_present <- genes_for_featureplot_input[genes_for_featureplot_input %in% rownames(data[["RNA"]])]
        if (length(genes_for_featureplot_present) == 0) {
            message(paste("None of the user-provided FeaturePlot genes were found in the data's RNA assay:", paste(genes_for_featureplot_input, collapse=", "), ". Skipping FeaturePlot."))
        }
    }
}

if(length(genes_for_featureplot_present) >= 1) {
    message(paste(Sys.time(), "- Generating FeaturePlots for user-provided genes:", paste(genes_for_featureplot_present, collapse=", ")))

    # SeuratExtend::FeaturePlot3
    message(paste(Sys.time(), "- Attempting SeuratExtend::FeaturePlot3..."))
    p_featureplot_seuratextend <- NULL # Initialize

    # Construct the call to FeaturePlot3 directly based on the number of available genes
    tryCatch({
        if (length(genes_for_featureplot_present) == 1) {
            p_featureplot_seuratextend <- SeuratExtend::FeaturePlot3(data, reduction = reduction_to_use, 
                                                                   feature.1 = genes_for_featureplot_present[1], 
                                                                   pt.size = 1)
        } else if (length(genes_for_featureplot_present) == 2) {
            p_featureplot_seuratextend <- SeuratExtend::FeaturePlot3(data, reduction = reduction_to_use, 
                                                                   feature.1 = genes_for_featureplot_present[1], 
                                                                   feature.2 = genes_for_featureplot_present[2], 
                                                                   pt.size = 1)
        } else if (length(genes_for_featureplot_present) >= 3) {
            # Plot up to 3 genes if more are provided
            p_featureplot_seuratextend <- SeuratExtend::FeaturePlot3(data, reduction = reduction_to_use, 
                                                                   feature.1 = genes_for_featureplot_present[1], 
                                                                   feature.2 = genes_for_featureplot_present[2], 
                                                                   feature.3 = genes_for_featureplot_present[3], 
                                                                   pt.size = 1)
        }
    }, error = function(e) {
        message(paste("Error generating SeuratExtend::FeaturePlot3:", e$message))
        # p_featureplot_seuratextend remains NULL
    })
    message(paste(Sys.time(), "- Finished SeuratExtend::FeaturePlot3 attempt."))

    if(!is.null(p_featureplot_seuratextend)){
        # Number of features actually plotted by FeaturePlot3 (up to 3)
        num_plotted_fp3 <- min(length(genes_for_featureplot_present), 3)
        
        plot_height_fp3 <- 7
        plot_width_fp3 <- 7 * num_plotted_fp3

        if(plot_width_fp3 > 0) { 
             message(paste(Sys.time(), "- Saving SeuratExtend::FeaturePlot3..."))
             ggsave(file.path(base_output_path, paste0("04a_basicviz_featureplot_seuratextend_", reduction_to_use, ".png")), plot = p_featureplot_seuratextend, width = plot_width_fp3, height = plot_height_fp3, limitsize = FALSE)
             message(paste0(Sys.time(), "- SeuratExtend::FeaturePlot3 saved using ", reduction_to_use, "."))
        }
    }

    # SCpubr::do_FeaturePlot
    message(paste(Sys.time(), "- Attempting SCpubr::do_FeaturePlot..."))
    p_featureplot_scpubr <- tryCatch({
        SCpubr::do_FeaturePlot(data, features = genes_for_featureplot_present, reduction = reduction_to_use, plot.axes = TRUE)
    }, error = function(e) {
        message(paste("Error generating SCpubr::do_FeaturePlot:", e$message)); NULL
    })
    message(paste(Sys.time(), "- Finished SCpubr::do_FeaturePlot attempt."))
    if(!is.null(p_featureplot_scpubr)){
        num_features_to_plot_scpubr <- length(genes_for_featureplot_present)
        plot_height_scpubr <- max(7, 4 * ceiling(num_features_to_plot_scpubr / 2))
        plot_width_scpubr <- 12
        if (num_features_to_plot_scpubr == 1) plot_width_scpubr <- 7
        message(paste(Sys.time(), "- Saving SCpubr::do_FeaturePlot..."))
        ggsave(file.path(base_output_path, paste0("04a_basicviz_featureplot_scpubr_", reduction_to_use, ".png")), plot = p_featureplot_scpubr, width = plot_width_scpubr, height = plot_height_scpubr, limitsize = FALSE) 
        message(paste0(Sys.time(), "- SCpubr::do_FeaturePlot saved using ", reduction_to_use, "."))
    }
} else {
    if (!is.null(featureplot_genes_env_var) && featureplot_genes_env_var != "") {
         message("FeaturePlot skipped as no valid or present genes were specified by the user via IMPACT_SC_FEATUREPLOT_GENES.")
    }
}


# --- Dot Plot ---
user_dotplot_groups <- list() 

if (is.null(dotplot_genes_json_env_var) || dotplot_genes_json_env_var == "" || dotplot_genes_json_env_var == "[]") {
    message("IMPACT_SC_DOTPLOT_GENES_JSON environment variable not set, empty, or '[]'. Skipping DotPlot based on user input.")
} else {
    user_dotplot_groups <- tryCatch({
        # simplifyDataFrame = TRUE is default for fromJSON, but explicit for clarity
        # This will parse '[{"name":"N", "genes":["G1","G2"]}]' into a data.frame
        parsed_json <- jsonlite::fromJSON(dotplot_genes_json_env_var, simplifyDataFrame = TRUE) 
        
        if (length(parsed_json) == 0) { 
             message("Parsed DotPlot JSON from IMPACT_SC_DOTPLOT_GENES_JSON is empty. Skipping user-defined DotPlot.")
             list()
        } else if (is.data.frame(parsed_json) && "name" %in% colnames(parsed_json) && "genes" %in% colnames(parsed_json)) {
            # Handle data.frame case (e.g., from '[{"name":"N1","genes":["G1","G2"]},{"name":"N2","genes":["G3"]}]')
            message("DotPlot JSON parsed as data.frame. Converting to named list.")
            # Ensure 'genes' column is a list of character vectors for setNames
            if (!is.list(parsed_json$genes) || !all(sapply(parsed_json$genes, is.character))) {
                 message("Warning: 'genes' column in parsed DotPlot JSON data.frame is not a list of character vectors as expected. Attempting to proceed.")
                 # Potentially try to coerce, or let it fail if structure is too different
            }
            # The 'genes' column in the data.frame is already a list of character vectors
            # e.g., parsed_json$genes might be list(c("CD3D", "CD19"))
            # parsed_json$name would be c("B")
            # setNames will correctly create list(B = c("CD3D", "CD19"))
            stats::setNames(parsed_json$genes, as.character(parsed_json$name))
        } else if (is.list(parsed_json) && !is.data.frame(parsed_json) && !is.null(names(parsed_json)) && all(sapply(parsed_json, function(x) is.character(x) || (is.list(x) && all(sapply(x, is.character)))))) {
            # This branch is for a simple named list format: e.g. '{"group1": ["gA", "gB"]}'
            # or where values could be single strings too (though less common for DotPlot genes)
            message("DotPlot input appears to be a simple named list of gene vectors. Using as is.")
            parsed_json 
        } else if (is.list(parsed_json) && !is.data.frame(parsed_json) && length(parsed_json) > 0 && all(sapply(parsed_json, function(x) is.list(x) && "name" %in% names(x) && "genes" %in% names(x)))) {
            # This branch is for the standard format: list of lists, e.g. from '[{"name":"group1", "genes":["gA","gB"]}]' if simplifyDataFrame=FALSE
            message("DotPlot JSON parsed as list of lists. Converting to named list.")
            stats::setNames(lapply(parsed_json, function(x) x$genes), lapply(parsed_json, function(x) x$name))
        } else {
            message("DotPlot JSON from IMPACT_SC_DOTPLOT_GENES_JSON is not in a recognized format. Skipping user-defined DotPlot. Expected format: '[{\"name\":\"group1\", \"genes\":[\"geneA\",\"geneB\"]},...]' or a simple named list in JSON.")
            print("Received JSON string for DotPlot:")
            print(dotplot_genes_json_env_var)
            print("Parsed structure:")
            print(str(parsed_json))
            list()
        }
    }, error = function(e) {
        message(paste("Error parsing DotPlot JSON from IMPACT_SC_DOTPLOT_GENES_JSON:", e$message, ". Skipping user-defined DotPlot."))
        list()
    })
}


dotplot_features_present <- list()
if (length(user_dotplot_groups) > 0) {
    dotplot_features_present <- lapply(user_dotplot_groups, function(g_list) {
        if(is.character(g_list)) { 
            g_list[g_list %in% rownames(data[["RNA"]])]
        } else {
            message(paste("Warning: A gene group in DotPlot input was not a character vector. Group name:", names(user_dotplot_groups)[which(sapply(user_dotplot_groups, identical, g_list))]))
            character(0) 
        }
    })
    dotplot_features_present <- Filter(function(x) length(x) > 0, dotplot_features_present)
}


if (length(dotplot_features_present) > 0) {
  message(paste(Sys.time(), "- Generating DotPlots with the following user-provided gene groups (after filtering for presence in data):"))
  print(dotplot_features_present)

  # SeuratExtend::DotPlot2
  message(paste(Sys.time(), "- Attempting SeuratExtend::DotPlot2..."))
  p_dotplot_seuratextend <- tryCatch({
    DotPlot2(data, features = dotplot_features_present, group.by = "cell_type")
  }, error = function(e) {
    message(paste("Error generating SeuratExtend::DotPlot2:", e$message)); NULL
  })
  message(paste(Sys.time(), "- Finished SeuratExtend::DotPlot2 attempt."))
  if(!is.null(p_dotplot_seuratextend)){
    num_total_genes_dp_se <- length(unlist(dotplot_features_present))
    plot_width_dp_se <- max(10, 1 + 0.5 * num_total_genes_dp_se)
    plot_height_dp_se <- max(7, 1 + 0.3 * nlevels(as.factor(data$cell_type)))
    message(paste(Sys.time(), "- Saving SeuratExtend::DotPlot2..."))
    ggsave(file.path(base_output_path, "04a_basicviz_dotplot_seuratextend.png"), plot = p_dotplot_seuratextend, width = plot_width_dp_se, height = plot_height_dp_se, limitsize = FALSE)
    message(paste0(Sys.time(), "- SeuratExtend::DotPlot2 saved."))
  }

  # SCpubr::do_DotPlot
  message(paste(Sys.time(), "- Attempting SCpubr::do_DotPlot..."))
  p_dotplot_scpubr <- tryCatch({
    SCpubr::do_DotPlot(data, features = dotplot_features_present)
  }, error = function(e) {
    message(paste("Error generating SCpubr::do_DotPlot:", e$message)); NULL
  })
  message(paste(Sys.time(), "- Finished SCpubr::do_DotPlot attempt."))
  if(!is.null(p_dotplot_scpubr)){
    num_total_genes_dp_sc <- length(unlist(dotplot_features_present))
    plot_width_dp_sc <- max(10, 1 + 0.5 * num_total_genes_dp_sc)
    plot_height_dp_sc <- max(7, 1 + 0.3 * nlevels(as.factor(data$cell_type)))
    message(paste(Sys.time(), "- Saving SCpubr::do_DotPlot..."))
    ggsave(file.path(base_output_path, "04a_basicviz_dotplot_scpubr.png"), plot = p_dotplot_scpubr, width = plot_width_dp_sc, height = plot_height_dp_sc, limitsize = FALSE)
    message(paste0(Sys.time(), "- SCpubr::do_DotPlot saved."))
  }
} else {
    if (!is.null(dotplot_genes_json_env_var) && dotplot_genes_json_env_var != "" && dotplot_genes_json_env_var != "[]") {
         message("DotPlot skipped as no valid gene groups were provided by the user via IMPACT_SC_DOTPLOT_GENES_JSON, or no genes within those groups were found in the data.")
    }
}


# --- Heatmap of top N variable genes per cluster (scRNAtoolVis) ---
DefaultAssay(data) <- "RNA"
if (length(VariableFeatures(data)) == 0) {
    message("Finding variable features for RNA assay before CalcStats.")
    data <- FindVariableFeatures(data, assay = "RNA", verbose = FALSE) 
}

if (length(VariableFeatures(data)) > 0 && "cell_type" %in% colnames(data@meta.data)) { 
    message("Calculating z-scores for top variable features per cluster using scRNAtoolVis::CalcStats.")

        message("Scaling data before CalcStats/Heatmap.")
        all_genes_for_scaling <- rownames(data[["RNA"]])
        data <- ScaleData(data, features = all_genes_for_scaling, assay = "RNA", verbose = FALSE)
    

    genes_zscore_scrnatoolvis <- tryCatch({
        CalcStats( # Added namespace
            data,
            features = VariableFeatures(data), 
            group.by = "cell_type",      
            order = "p",                    
            n = 4
        )
    }, error = function(e) {
        message(paste("Error in scRNAtoolVis::CalcStats:", e$message)); NULL
    })

    if (!is.null(genes_zscore_scrnatoolvis) && nrow(genes_zscore_scrnatoolvis) > 0) {
        p_heatmap_scrnatoolvis <- tryCatch({
            Heatmap(genes_zscore_scrnatoolvis, lab_fill = "zscore") # Added namespace
        }, error = function(e) {
            message(paste("Error generating scRNAtoolVis::Heatmap:", e$message)); NULL
        })
        if(!is.null(p_heatmap_scrnatoolvis)){
            ggsave(file.path(base_output_path, "04a_basicviz_heatmap_scrnatoolvis_zscore.png"), plot = p_heatmap_scrnatoolvis, width = 10, height = max(8, 0.2 * nrow(genes_zscore_scrnatoolvis)) , limitsize = FALSE)
            message("scRNAtoolVis::Heatmap of z-scores saved.")
        }
    } else {
        message("No data returned from CalcStats or an error occurred; skipping scRNAtoolVis::Heatmap.")
    }
} else {message("Not enough variable features or 'cell_type' missing for scRNAtoolVis::CalcStats based heatmap.")}


# Save the Seurat object (data) after this module if subsequent modules need its state
intermediate_rds_path_04a <- file.path(base_output_path, "04a_data_after_basicviz.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04a)
  message(paste("Data object after basic visualization saved to:", intermediate_rds_path_04a))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.1:", e$message))
})

message("Finished Module 4.1: Basic Visualization")