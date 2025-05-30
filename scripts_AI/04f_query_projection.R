# IMPACT-sc Script: 04f_query_projection.R
# Purpose: Project a Query Dataset onto a Reference Embedding.
# THIS MODULE IS OPTIONAL AND REQUIRES A 'query.RDS' FILE.

# Increase expression limit to potentially avoid "evaluation nested too deeply" error
options(expressions = 10000)

# CRITICAL WARNING:
# The log files indicate significant R version and package version mismatches.
# For example, the script might be run with R 4.4.2, while packages like
# Seurat, ggplot2, harmony were built under R 4.4.3, and SingleR under R 4.5.0.
# Such incompatibilities are a VERY LIKELY CAUSE of errors like "node stack overflow"
# or other unexpected issues.
#
# STRONGLY RECOMMENDED:
# 1. Update your R environment to the latest stable version (or at least to match the
#    R version used to build your most recent packages, e.g., R 4.5.0 or higher).
# 2. After updating R, reinstall all key packages (Seurat, Bioconductor packages, etc.)
#    to ensure they are compiled correctly for your R version.
#
# The code changes below (e.g., compress=FALSE in saveRDS) are attempts to mitigate
# symptoms, but the underlying environment compatibility issue is paramount to resolve.

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
# Species: This is the species of the REFERENCE dataset ('data')
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

#### Module 4.7: Project Query onto Embedding ####
message("Starting Module 4.7: Project Query onto Embedding")

# Input:
#   - '04f_data_after_pseudotime.RDS' (or previous) as REFERENCE ('data')
#   - 'query.RDS' (USER-PROVIDED query dataset, path read from ENV variable)
# Output: Mapped query object, plots.

obj_prev_module_path <- file.path(base_output_path, "04f_data_after_pseudotime.RDS")
if (!file.exists(obj_prev_module_path)) {
  obj_prev_module_path <- file.path(base_output_path, "03_module3_final_annotated.RDS") # Fallback
  message(paste("04f output not found, loading from:", obj_prev_module_path))
}
if (!file.exists(obj_prev_module_path)) {
  stop(paste("Required reference Seurat object not found at:", obj_prev_module_path))
}
ref_data <- tryCatch(readRDS(obj_prev_module_path), error=function(e){
  message(paste("Error loading reference RDS:",e$message));NULL
})
if(is.null(ref_data)) stop("Failed to load reference data for Module 4.7.")
message("Loaded reference Seurat object (as 'ref_data') for Module 4.7.")
DefaultAssay(ref_data) <- "RNA"

# Path to query dataset, NOW READ FROM ENVIRONMENT VARIABLE
query_rds_path <- Sys.getenv("IMPACT_SC_QUERY_RDS_PATH", "") # Default to empty string if not set

if (query_rds_path == "" || !file.exists(query_rds_path)) {
  message(paste("Query dataset path 'IMPACT_SC_QUERY_RDS_PATH' not set, empty, or file '", query_rds_path, "' not found. Skipping Module 4.7.", sep=""))
} else {
  message(paste("Attempting to load query dataset from:", query_rds_path))
  query <- tryCatch(readRDS(query_rds_path), error=function(e){
    message(paste("Error loading query.RDS from", query_rds_path, ":", e$message)); NULL
  })
  if(is.null(query)) {
    message(paste("Failed to load query.RDS from", query_rds_path, ". Skipping Module 4.7."))
  } else {
    message("Loaded query dataset.")
    DefaultAssay(query) <- "RNA"

    query_species_param <- Sys.getenv("IMPACT_SC_QUERY_SPECIES", species)
    message(paste("Reference species is:", species, ". Query species set to:", query_species_param))

    if (query_species_param == "mouse" && species == "human") {
      message("Attempting to convert mouse genes in query to human orthologs.")
      mouse_genes <- rownames(query[["RNA"]])
      human_map <- tryCatch(homologene::mouse2human(mouse_genes), error=function(e){message(e);NULL})
      if(!is.null(human_map) && nrow(human_map) > 0){
        human_map <- human_map[!is.na(human_map$humanGene) & !duplicated(human_map$mouseGene), ]
        if(nrow(human_map) > 0){ # Check again after filtering
            query_counts_matrix <- GetAssayData(query, assay = "RNA", layer = "counts")[human_map$mouseGene, , drop=FALSE]
            rownames(query_counts_matrix) <- make.unique(human_map$humanGene)
            query <- CreateSeuratObject(counts = query_counts_matrix, meta.data = query@meta.data)
            message("Query genes converted to human orthologs.")
        } else { message("No valid orthologs found after filtering. Proceeding with original query genes.")}
      } else { message("homologene mouse2human failed or returned no results. Proceeding with original query genes.")}
    } else if (query_species_param == "human" && species == "mouse") {
      message("Attempting to convert human genes in query to mouse orthologs.")
      human_genes <- rownames(query[["RNA"]])
      mouse_map <- tryCatch(homologene::human2mouse(human_genes), error=function(e){message(e);NULL})
       if(!is.null(mouse_map) && nrow(mouse_map) > 0){
        mouse_map <- mouse_map[!is.na(mouse_map$mouseGene) & !duplicated(mouse_map$humanGene), ]
        if(nrow(mouse_map) > 0){ # Check again after filtering
            query_counts_matrix <- GetAssayData(query, assay = "RNA", layer = "counts")[mouse_map$humanGene, , drop=FALSE]
            rownames(query_counts_matrix) <- make.unique(mouse_map$mouseGene)
            query <- CreateSeuratObject(counts = query_counts_matrix, meta.data = query@meta.data)
            message("Query genes converted to mouse orthologs.")
        } else { message("No valid orthologs found after filtering. Proceeding with original query genes.")}
      } else { message("homologene human2mouse failed or returned no results. Proceeding with original query genes.")}
    }

    if(ncol(query) > 20){
        query_sce <- tryCatch(as.SingleCellExperiment(query, assay = "RNA"), error=function(e){message(e);NULL})
        if(!is.null(query_sce)){
            query_sce <- scDblFinder(query_sce, BPPARAM = BiocParallel::SerialParam())
            query <- as.Seurat(query_sce, counts = "counts", data = "logcounts")
            query <- subset(query, subset = scDblFinder.class == "singlet")
            message("Doublets removed from query dataset.")
        } else {message("Could not convert query to SCE for doublet removal.")}
    } else {message("Query has too few cells for scDblFinder, skipping doublet removal.")}


    query <- NormalizeData(query, verbose = FALSE)
    query <- FindVariableFeatures(query, verbose = FALSE)
    query <- ScaleData(query, verbose=FALSE) %>% RunPCA(verbose=FALSE)


    if (is.list(ref_data@assays$RNA) || (length(ref_data@assays$RNA@layers) > 0 && !is.null(names(ref_data@assays$RNA@layers)))) {
        ref_data <- JoinLayers(ref_data, assay="RNA")
    }
    DefaultAssay(ref_data) <- "RNA"
    if (length(VariableFeatures(ref_data))==0) ref_data <- FindVariableFeatures(ref_data, verbose=FALSE)
    if (!("pca" %in% Reductions(ref_data))) ref_data <- ScaleData(ref_data, verbose=FALSE) %>% RunPCA(verbose=FALSE)

    anchors <- tryCatch(FindTransferAnchors(
      reference = ref_data, query = query, dims = 1:min(30, ncol(ref_data[['pca']]), ncol(query[['pca']])),
      reference.reduction = "pca", verbose = FALSE
    ), error = function(e){message(paste("FindTransferAnchors failed:", e$message));NULL})

    # MODIFIED SECTION:
    # The following block directly executes the data transfer, UMAP, and mapping steps.
    # The original conditional check `if(!is.null(anchors))` and its `else` branch have been removed
    # as per user request. This means the script assumes `anchors` is non-NULL and valid.
    # If `FindTransferAnchors` failed and `anchors` is NULL, the following lines will likely error.

    predictions_transfer <- TransferData(anchorset = anchors, refdata = ref_data$cell_type, dims = 1:min(30, ncol(ref_data[['pca']])), verbose = FALSE)
    query <- AddMetaData(query, metadata = predictions_transfer)

    ref_reduction_for_umap_model <- "pca"
    ref_umap_name <- paste0("umap.", ref_reduction_for_umap_model)

    # Unconditionally run UMAP on reference data as per user request
    ref_data <- RunUMAP(ref_data, dims = 1:min(30, ncol(ref_data[[ref_reduction_for_umap_model]])),
                         reduction = ref_reduction_for_umap_model, return.model = TRUE,
                         reduction.name = ref_umap_name, verbose = FALSE)
    message(paste("Ran UMAP on reference (unconditionally) using reduction:", ref_reduction_for_umap_model, "and saved model as", ref_umap_name))


    query <- MapQuery(
      anchorset = anchors, reference = ref_data, query = query,
      refdata = list(celltype_ref = "cell_type"), reference.reduction = "pca",
      reduction.model = ref_umap_name, verbose = FALSE
    )
    message("Query mapped to reference UMAP space. New reduction in query: 'ref.umap'")

    p1_ref <- DimPlot(ref_data, reduction = ref_umap_name, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference Annotations")
    p2_query <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype_ref", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query Mapped (Predicted)")
    combined_map_plot <- p1_ref + p2_query
    ggsave(file.path(base_output_path, "04g_ref_query_umap_comparison.png"), plot = combined_map_plot, width = 16, height = 7)

    # Attempting to save RDS with compress = FALSE as a potential workaround for "node stack overflow"
    # However, resolving R/package version incompatibilities is the primary solution.
    saveRDS(query, file.path(base_output_path, "04g_query_mapped.RDS"), compress = FALSE)
    message("Mapped query object saved.")
    # END OF MODIFIED SECTION

  } # Closes if(!is.null(query))
} # Closes if (query_rds_path == "" || !file.exists(query_rds_path))

intermediate_rds_path_04g <- file.path(base_output_path, "04g_data_after_queryproj.RDS")
tryCatch({
  # Attempting to save RDS with compress = FALSE as a potential workaround
  saveRDS(ref_data, intermediate_rds_path_04g, compress = FALSE)
  message(paste("Reference data object after (attempted) Query Projection saved to:", intermediate_rds_path_04g))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.7:", e$message))
  warning("This error might be related to R/package version incompatibilities or object complexity.")
})

message("Finished Module 4.7: Project Query onto Embedding")
