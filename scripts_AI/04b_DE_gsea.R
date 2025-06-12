# IMPACT-sc Script: 04b_diffex_gsea.R
# Purpose: Perform Differential Expression and Gene Set Enrichment Analysis.

# --- Libraries ---
library(AnnotationDbi)
library(CARD)
# library(celldex) # Celldex is no longer used for downloading reference data
library(decoupleR)
library(dplyr)
library(ensembldb)
library(ggplot2)
library(ggpubr) # Still loaded as it might be used elsewhere or by VlnPlot2 dependencies
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
library(SeuratExtend) # VlnPlot2 might come from here
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

# Gene for violin plot, read from environment variable
user_input_gene_for_vln_plot <- Sys.getenv("IMPACT_SC_DE_GENE", "") # Default to empty string if not set

if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) { library(org.Hs.eg.db); org_db_object <- org.Hs.eg.db } else { stop("org.Hs.eg.db not installed.") }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) { library(org.Mm.eg.db); org_db_object <- org.Mm.eg.db } else { stop("org.Mm.eg.db not installed.") }
} else { stop("Species not supported.") }
options(future.globals.maxSize = 4 * 1024^3)
# --- End Script Parameters ---

#### Module 4.2: Differential Expression & GSEA ####
message("Starting Module 4.2: Differential Expression & GSEA")

# Input: '04a_data_after_basicviz.RDS' (or '03_module3_final_annotated.RDS' if 04a is skipped)
# Output: DGE results, plots. Modifies 'data' object.
obj_prev_module_path <- file.path(base_output_path, "04a_data_after_basicviz.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS")
  message(paste("04a output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required input Seurat object not found at:", obj_prev_module_path))
}
data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading input RDS:",e$message));NULL
})
if(is.null(data)) stop("Failed to load data for Module 4.2.")
message("Loaded Seurat object for Module 4.2.")
DefaultAssay(data) <- "RNA"

# Gene expression across cell types (Violin/Box plots) using user-provided gene
if (is.null(user_input_gene_for_vln_plot) || user_input_gene_for_vln_plot == "") {
    message("IMPACT_SC_DE_GENE environment variable not set or empty. Skipping user-defined gene expression violin plot.")
} else {
    sanitized_gene_name <- gsub("[^A-Za-z0-9_.-]", "_", user_input_gene_for_vln_plot)
    if (user_input_gene_for_vln_plot %in% rownames(data[["RNA"]])) {
        p_vln_plot <- VlnPlot2(data, 
                               features = user_input_gene_for_vln_plot, 
                               group.by = "cell_type")
                               # stat.method = "wilcox.test" # Temporarily removed to avoid filter() namespace error
        
        # Adding title separately as it's a common ggplot2 modification
        p_vln_plot <- p_vln_plot + labs(title = paste("Expression of", user_input_gene_for_vln_plot))

        plot_filename <- paste0("04b_dge_vlnplot_", sanitized_gene_name, ".png") # Consistent naming
        ggsave(file.path(base_output_path, plot_filename), plot = p_vln_plot, width = max(7, 1.5*nlevels(as.factor(data$cell_type))), height = 6)
        message(paste("User-defined gene expression violin plot for '", user_input_gene_for_vln_plot, "' saved as ", plot_filename, ".", sep=""))
    } else {
        message(paste("Gene '", user_input_gene_for_vln_plot, "' provided via IMPACT_SC_DE_GENE was not found in the RNA assay rownames. Skipping violin plot.", sep=""))
    }
}

p_vln_plot2=SCpubr::do_ViolinPlot(data,features =user_input_gene_for_vln_plot,group.by ="cell_type" )
 plot_filename <- paste0("04b_dge_vlnplot2_", sanitized_gene_name, ".png") # Consistent naming
ggsave(file.path(base_output_path, plot_filename), plot = p_vln_plot2, width = max(7, 1.5*nlevels(as.factor(data$cell_type))), height = 6)
p_box_plot=SCpubr::do_BoxPlot(data,feature = user_input_gene_for_vln_plot,group.by = "cell_type")
 plot_filename <- paste0("04b_dge_boxplot_", sanitized_gene_name, ".png") # Consistent naming
ggsave(file.path(base_output_path, plot_filename), plot = p_box_plot, width = max(7, 1.5*nlevels(as.factor(data$cell_type))), height = 6)



# Identify marker genes
if (!("cell_type" %in% colnames(data@meta.data))) {
    warning("'cell_type' column not found in metadata. Skipping DGE and Volcano plot.")
    all_markers <- NULL # Ensure all_markers is NULL if cell_type is missing
} else {
    Idents(data) <- "cell_type" # Make sure Idents are set to cell_type for FindAllMarkers
    all_markers <- tryCatch({
        FindAllMarkers(
            data,
            min.cells.feature = 10, # Minimum number of cells a gene must be detected in
            min.cells.group = 10,   # Minimum number of cells in either of the two groups being compared
            min.pct = 0.25,         # Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations
            logfc.threshold = 0.25, # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
            verbose = FALSE         # Print progress bar
            # only.pos = TRUE # Consider if you only want positive markers, removed for more general use
        )
    }, error = function(e) {
        message(paste("FindAllMarkers failed:", e$message)); NULL
    })
}

if(!is.null(all_markers) && nrow(all_markers) > 0){
  write.csv(all_markers, file.path(base_output_path, "04b_all_celltype_markers.csv"), row.names = FALSE)
  message("Differential expression markers saved.")
  
  # Ensure 'cluster' column (or the one used for grouping by FindAllMarkers) and 'avg_log2FC' exist
  if ("cluster" %in% colnames(all_markers) && "avg_log2FC" %in% colnames(all_markers)) {
      top10_markers_m4_2 <- all_markers %>% # Corrected variable name from m4.2 to m4_2
                              group_by(cluster) %>%
                              top_n(n = 10, wt = avg_log2FC)
      print(head(top10_markers_m4_2))
  } else {
      message("Could not generate top 10 markers table: 'cluster' or 'avg_log2FC' column missing from FindAllMarkers output.")
  }

  # Volcano plot using scRNAtoolVis::jjVolcano with all_markers
  # Ensure all_markers has the required columns: p_val_adj and avg_log2FC
  if ("p_val_adj" %in% colnames(all_markers) && "avg_log2FC" %in% colnames(all_markers)) {
      volcano_plot_m4_2 <- jjVolcano( # Corrected variable name
          diffData = all_markers,       # Using the entire all_markers table as per user's correction
          log2FC.cutoff = 0.5,          # Cutoff for log2 fold change
          FDR.cutoff = 0.05,            # Cutoff for adjusted p-value
          topGeneN = 5,                 # Number of top genes to label
          fontface = 'italic'
      )
      ggsave(file.path(base_output_path, "04b_volcano_plot_all_markers.png"), plot = volcano_plot_m4_2, width = 10, height = 8) 
      message("Volcano plot for all markers saved as 04b_volcano_plot_all_markers.png")
  } else {
      message("Required columns (p_val_adj, avg_log2FC) not found in all_markers for jjVolcano.")
  }
} else {
  message("No markers found or FindAllMarkers failed. Skipping marker table and volcano plot.")
}


# GSEA with SeuratExtend (GeneSetAnalysisGO, RenameGO, GSEAplot) - Commented out due to dependency
message("GSEA with SeuratExtend is commented out due to specific dependencies. Standard GSEA tools like clusterProfiler or fgsea can be used as an alternative. These tools typically require a ranked list of genes from DGE results.")

intermediate_rds_path_04b <- file.path(base_output_path, "04b_data_after_diffex.RDS")
tryCatch({
  saveRDS(data, intermediate_rds_path_04b)
  message(paste("Data object after DGE saved to:", intermediate_rds_path_04b))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.2:", e$message))
})

message("Finished Module 4.2: Differential Expression & GSEA")
