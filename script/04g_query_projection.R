# IMPACT-sc Script: 04g_query_projection.R
# Purpose: Project a Query Dataset onto a Reference Embedding.
# THIS MODULE IS OPTIONAL AND REQUIRES A 'query.RDS' FILE.

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
species <- "human" # This is the species of the REFERENCE dataset ('data')
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

#### Module 4.7: Project Query onto Embedding ####
message("Starting Module 4.7: Project Query onto Embedding")

# Input: 
#   - '04f_data_after_pseudotime.RDS' (or previous) as REFERENCE ('data')
#   - 'query.RDS' (USER-PROVIDED query dataset, path relative to project root)
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

query_rds_path <- "../query.RDS" # Path to query dataset, relative to 'scripts/'
if (!file.exists(query_rds_path)) {
  message(paste("Query dataset 'query.RDS' not found at:", query_rds_path, ". Skipping Module 4.7."))
} else {
  query <- tryCatch(readRDS(query_rds_path), error=function(e){
    message(paste("Error loading query.RDS:", e$message)); NULL
  })
  if(is.null(query)) {
    message("Failed to load query.RDS. Skipping Module 4.7.")
  } else {
    message("Loaded query dataset.")
    DefaultAssay(query) <- "RNA" # Assume RNA assay for query

    # **USER ACTION: Define query_species if different from reference**
    query_species_param <- "mouse" # Example: if your query is mouse
    message(paste("Reference species is:", species, ". Query species assumed to be:", query_species_param))

    # Gene name conversion if species differ (example: mouse query, human reference)
    if (query_species_param == "mouse" && species == "human") {
      message("Attempting to convert mouse genes in query to human orthologs.")
      mouse_genes <- rownames(query[["RNA"]])
      human_map <- tryCatch(homologene::mouse2human(mouse_genes), error=function(e){message(e);NULL})
      if(!is.null(human_map)){
        human_map <- human_map[!is.na(human_map$humanGene) & !duplicated(human_map$mouseGene), ]
        if(nrow(human_map) > 0){
            query_counts_matrix <- GetAssayData(query, assay = "RNA", layer = "counts")[human_map$mouseGene, , drop=FALSE]
            rownames(query_counts_matrix) <- make.unique(human_map$humanGene)
            query <- CreateSeuratObject(counts = query_counts_matrix, meta.data = query@meta.data)
            message("Query genes converted to human orthologs.")
        } else { message("No orthologs found or mapping failed. Proceeding with original query genes.")}
      } else { message("homologene mouse2human failed. Proceeding with original query genes.")}
    } else if (query_species_param == "human" && species == "mouse") {
      message("Attempting to convert human genes in query to mouse orthologs.")
      human_genes <- rownames(query[["RNA"]])
      mouse_map <- tryCatch(homologene::human2mouse(human_genes), error=function(e){message(e);NULL})
       if(!is.null(mouse_map)){
        mouse_map <- mouse_map[!is.na(mouse_map$mouseGene) & !duplicated(mouse_map$humanGene), ]
        if(nrow(mouse_map) > 0){
            query_counts_matrix <- GetAssayData(query, assay = "RNA", layer = "counts")[mouse_map$humanGene, , drop=FALSE]
            rownames(query_counts_matrix) <- make.unique(mouse_map$mouseGene)
            query <- CreateSeuratObject(counts = query_counts_matrix, meta.data = query@meta.data)
            message("Query genes converted to mouse orthologs.")
        } else { message("No orthologs found or mapping failed. Proceeding with original query genes.")}
      } else { message("homologene human2mouse failed. Proceeding with original query genes.")}
    }
    
    # Doublet removal for query
    if(ncol(query) > 20){ # scDblFinder needs sufficient cells
        query_sce <- tryCatch(as.SingleCellExperiment(query, assay = "RNA"), error=function(e){message(e);NULL})
        if(!is.null(query_sce)){
            query_sce <- scDblFinder(query_sce, BPPARAM = BiocParallel::SerialParam()) # SerialParam for non-parallel
            query <- as.Seurat(query_sce, counts = "counts", data = "logcounts") # or data=NULL
            query <- subset(query, subset = scDblFinder.class == "singlet")
            message("Doublets removed from query dataset.")
        } else {message("Could not convert query to SCE for doublet removal.")}
    } else {message("Query has too few cells for scDblFinder, skipping doublet removal.")}


    query <- NormalizeData(query, verbose = FALSE)
    query <- FindVariableFeatures(query, verbose = FALSE)
    # Ensure query has PCA if FindTransferAnchors needs it (often does implicitly)
    query <- ScaleData(query, verbose=FALSE) %>% RunPCA(verbose=FALSE)


    # Prepare reference data
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

    if(!is.null(anchors)){
      predictions_transfer <- TransferData(anchorset = anchors, refdata = ref_data$cell_type, dims = 1:min(30, ncol(ref_data[['pca']])), verbose = FALSE)
      query <- AddMetaData(query, metadata = predictions_transfer)
      
      ref_reduction_for_umap_model <- "pca" # UMAP model should be built on the same reduction used for anchors
      ref_umap_name <- paste0("umap.", ref_reduction_for_umap_model)

      if (!(ref_umap_name %in% Reductions(ref_data)) || is.null(ref_data[[ref_umap_name]]@misc$model) ) {
          ref_data <- RunUMAP(ref_data, dims = 1:min(30, ncol(ref_data[[ref_reduction_for_umap_model]])), 
                           reduction = ref_reduction_for_umap_model, return.model = TRUE, 
                           reduction.name = ref_umap_name, verbose = FALSE)
          message(paste("Ran UMAP on reference using reduction:", ref_reduction_for_umap_model, "and saved model as", ref_umap_name))
      }

      query <- MapQuery(
        anchorset = anchors, reference = ref_data, query = query,
        refdata = list(celltype_ref = "cell_type"), reference.reduction = "pca", 
        reduction.model = ref_umap_name, verbose = FALSE
      )
      message("Query mapped to reference UMAP space. New reduction in query: 'ref.umap'")
      
      # Visualization
      p1_ref <- DimPlot(ref_data, reduction = ref_umap_name, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference Annotations")
      p2_query <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype_ref", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query Mapped (Predicted)")
      combined_map_plot <- p1_ref + p2_query
      ggsave(file.path(base_output_path, "04g_ref_query_umap_comparison.png"), plot = combined_map_plot, width = 16, height = 7)
      
      saveRDS(query, file.path(base_output_path, "04g_query_mapped.RDS"))
      message("Mapped query object saved.")
    } else { message("Skipping data transfer and mapping due to anchor finding failure.")}
    # Update the main 'data' object if this script is the last one modifying it before final save.
    # For now, query projection modifies 'query', not 'data'.
  }
}
# Save the reference 'data' object if it was modified (e.g., UMAP model added)
intermediate_rds_path_04g <- file.path(base_output_path, "04g_data_after_queryproj.RDS")
tryCatch({
  saveRDS(ref_data, intermediate_rds_path_04g) # Save ref_data as it might have new UMAP model
  message(paste("Reference data object after (attempted) Query Projection saved to:", intermediate_rds_path_04g))
}, error = function(e) {
  warning(paste("Error saving data object after Module 4.7:", e$message))
})

message("Finished Module 4.7: Project Query onto Embedding")
