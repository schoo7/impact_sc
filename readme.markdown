![Icon](impact_sc_logo.png)
# IMPACT-sc: Integrated and Modular Pipeline for Analysis of Cellular Transcriptomics

## 1. Introduction

Welcome to **IMPACT-sc (Integrated and Modular Pipeline for Analysis of Cellular Transcriptomics)**. This document outlines a comprehensive workflow for single-cell RNA sequencing (scRNA-seq) data analysis. IMPACT-sc is organized into several modules, guiding users from raw data processing, quality control, normalization, and dimensionality reduction, through batch correction, sophisticated cell type annotation, and a wide array of downstream analyses. These downstream modules include differential expression, pathway and transcription factor activity analysis, copy number variation prediction, spatial data integration, and more.

This version of IMPACT-sc incorporates improved error handling and species specificity, making it adaptable for human or mouse datasets.

**Key Features of IMPACT-sc:**
* **I**ntegrated: Combines R and Python tools (e.g., Seurat, Cell2Sentence) within a cohesive workflow.
* **M**odular: Structured for clarity, allowing step-by-step execution or selection of specific analytical modules.
* **P**ipeline for **A**nalysis of **C**ellular **T**ranscriptomics: Covers a broad spectrum of common and advanced scRNA-seq analysis tasks.
* Adaptable: Designed for human or mouse species with enhanced error handling.
* Comprehensive: Includes modules for advanced visualizations, pseudotime analysis (via noted dependencies), and GSEA (via noted dependencies).

## 2. Setup and Installation

### 2.1. Create a Virtual Environment

It's highly recommended to use a virtual environment to manage dependencies. We'll use `conda` to create an environment named `impact_sc` that can handle both R and Python packages.

```{bash create-conda-env, eval=FALSE}
# Create a conda environment
conda create -n impact_sc python=3.9 r-base=4.2 -y

# Activate the environment
conda activate impact_sc

# Install R and Python package managers within conda
conda install -c conda-forge mamba -y # Faster package installer
mamba install -c conda-forge r-essentials r-reticulate r-devtools -y
mamba install -c conda-forge pip -y
```

### 2.2. Install R Packages

Once the `impact_sc` environment is activated and R is set up, install the required R packages.

```{r install-r-packages, eval=FALSE}
# CRAN packages
install.packages(c(
  "AnnotationDbi", "dplyr", "ggplot2", "ggpubr", "homologene", "Matrix", 
  "patchwork", "presto", "reticulate", "rjags", "scater", "stringr", 
  "tibble", "tidyr", "viridis", "reshape2", "remotes", "pheatmap"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Common Bioconductor packages
BiocManager::install(c(
  "celldex", "decoupleR", "ensembldb", "infercnv", "msigdbr",
  "scDblFinder", "scRNAtoolVis", "SCpubr", "Seurat", "SeuratDisk", 
  "SingleCellExperiment", "SingleR", "SpatialExperiment", "scran", "UCell"
))

# Species-specific annotation databases (install both, script will load one)
BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db"))

# GitHub packages
# Ensure JAGS is installed for CARD: https://sourceforge.net/projects/mcmc-jags/files/
remotes::install_github("YingMa0107/CARD")

# SeuratExtend: This package seems to be custom or from a specific source not on CRAN/Bioc.
# If you have the GitHub repository URL, install it like this:
# remotes::install_github("developer/SeuratExtend")
```

**Note on SeuratExtend**: IMPACT-sc references functions like `DimPlot2`, `FeaturePlot3`, `DotPlot2`, `VlnPlot2`, `CalcStats`, `Heatmap`, `GeneSetAnalysisGO`, `RenameGO`, `GSEAplot`, `Palantir.RunDM`, `Palantir.Pseudotime`, `GeneTrendCurve.Palantir`, `GeneTrendHeatmap.Palantir`, and `ClusterDistrBar`. These appear to be from a package like `SeuratExtend` or a custom script collection. Ensure this package or the source of these functions is correctly installed. If `SeuratExtend` is not publicly available, these parts of the script will not work without its source code.

### 2.3. Install Python Packages

Install the required Python packages using `pip` within the `impact_sc` conda environment.

```{bash install-python-packages, eval=FALSE}
# Terminal within the 'impact_sc' conda environment
pip install anndata pandas scanpy cell2sentence numpy
```

**Note on cell2sentence**: Ensure the `cell2sentence` package is available via `pip`. The model `C2S-Pythia-410m-cell-type-prediction` will be downloaded by the `cell2sentence` library when first used.

### 2.4. External Files

- **Gene Order File for InferCNV**: The script uses a placeholder like `"hg38_gencode_v27.txt"`. This file is species-specific (e.g., human GRCh38 or mouse GRCm39). Ensure you have the correct file for your chosen species and update the path.
- **Input Data**: `ori.RDS` (your pre-processed Seurat object or object with raw counts), `query.RDS` (for Module 4.7), `Acinar_Cell_Carcinoma.RDS` (for Module 4.8).
- **Bulk Counts (Optional)**: `counts.csv` if running the bulk data TF analysis part in Module 4.3.

### 2.5. General Recommendations

- **Parameterization**: Key parameters like species are now at the top.
- **File Paths**: Use relative paths or define a base path variable.
- **Memory Management**: `options(future.globals.maxSize = ...)` is used.

## 3. Script Overview

IMPACT-sc is divided into the following main modules:
- **Module 1**: Data Processing - Loads raw data, performs QC, normalization, scaling, and PCA.
- **Module 2**: Dimensionality Reduction & Integration - Includes Harmony batch correction and Cell2Sentence AI-driven cell embeddings.
- **Module 3**: Cell Type Annotation - Employs Seurat clustering, Cell2Sentence predictions, and SingleR.
- **Module 4**: Downstream Analysis - A suite of analyses including:
  - Basic Visualization
  - Differential Expression & GSEA (GSEA via noted dependencies)
  - Pathway and TF Activity Scores (DecoupleR)
  - Gene Scores (UCell)
  - CNV Prediction (InferCNV)
  - Pseudotime Analysis (via noted SeuratExtend/Palantir dependencies)
  - Query Projection
  - Spatial Data Annotation (CARD)

## 4. Detailed Script Breakdown

### 4.0. Load Libraries and Set Parameters

```{r load-libraries, eval=FALSE}
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
# library(SeuratExtend) # Ensure this is available
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
```

```{r set-parameters, eval=FALSE}
# Specify the species: "human" or "mouse"
species <- "human" # MODIFY THIS PARAMETER AS NEEDED

# Set a seed for reproducibility
set.seed(1234)

# Define base paths
base_output_path <- "./output/"
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}

# Load species-specific annotation database
if (species == "human") {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    library(org.Hs.eg.db)
    org_db_object <- org.Hs.eg.db
    message("Loaded org.Hs.eg.db for human.")
  } else {
    stop("org.Hs.eg.db not installed. Please install it via BiocManager.")
  }
} else if (species == "mouse") {
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    org_db_object <- org.Mm.eg.db
    message("Loaded org.Mm.eg.db for mouse.")
  } else {
    stop("org.Mm.eg.db not installed. Please install it via BiocManager.")
  }
} else {
  stop("Species not supported by IMPACT-sc. Please choose 'human' or 'mouse'.")
}
```

### 4.1. Module 1: Data Processing

**Purpose**: Load raw count data, create a Seurat object, and perform initial quality control and preprocessing.

```{r module1-data-processing, eval=FALSE}
# Option 1: Load Demo Data from 10X Genomics (commented out)
# ...

# Option 2: Load Your Own Data
obj_rds_path <- "ori.RDS"
obj <- tryCatch({
  readRDS(obj_rds_path)
}, error = function(e) {
  stop(paste("Error loading ori.RDS:", e$message))
  return(NULL)
})
if (is.null(obj)) stop("Failed to load input data for IMPACT-sc.")

if (length(unique(obj$orig.ident)) > 1) {
    message("Multiple samples detected by 'orig.ident'.")
}
if (inherits(obj[["RNA"]], "list") || (is.list(obj@assays$RNA) && inherits(obj@assays$RNA[[1]], "Assay"))) {
    message("RNA assay is split. Joining layers for initial QC and processing.")
    obj <- JoinLayers(obj, assay="RNA")
}

# Perform basic quality control
mt_pattern <- if (species == "human") "^MT-" else if (species == "mouse") "^mt-" else "^MT-" # Default to human
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
message(paste("Calculated mitochondrial percentage using pattern:", mt_pattern))

# QC violin plots
qc_plot_path <- file.path(base_output_path, "qc_violin_plot.png")
png(qc_plot_path, width=10, height=6, units="in", res=300)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
message("QC violin plot saved to: ", qc_plot_path)

# Filter low-quality cells
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Standard single-cell processing
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Save processed Seurat object
module1_output_path <- file.path(base_output_path, "module1_processed.RDS")
tryCatch({
  saveRDS(obj, module1_output_path)
  message("Module 1 output saved to: ", module1_output_path)
}, error = function(e) {
  warning(paste("Error saving Module 1 output:", e$message))
})
```

### 4.2. Module 2.1: Harmony Batch Correction

```{r module2-harmony, eval=FALSE}
obj_module1_path <- file.path(base_output_path, "module1_processed.RDS")
if (file.exists(obj_module1_path)) {
  obj <- tryCatch(readRDS(obj_module1_path), error = function(e) {
    message(paste("Error loading Module 1 output:", e$message))
    NULL
  })
  if (is.null(obj)) stop("Failed to load data for Module 2.1")
} else {
  stop("Module 1 output (module1_processed.RDS) not found.")
}

if (length(unique(obj$orig.ident)) > 1) {
    message("Splitting RNA assay by orig.ident for Harmony integration if not already split.")
    is_v5_layered <- !is.null(obj@assays$RNA$counts) && inherits(obj@assays$RNA$counts, "dgCMatrix") &&
                     length(obj@assays$RNA@layers) > 0 && inherits(obj@assays$RNA@layers[[1]], "dgCMatrix")
    if (!is_v5_layered && !is.list(obj@assays$RNA)) {
         obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    }
    obj <- IntegrateLayers(
      object = obj, method = HarmonyIntegration, orig.reduction = "pca",
      new.reduction = "harmony", group.by.vars = "orig.ident", verbose = TRUE
    )
} else {
    message("Only one 'orig.ident' found. Harmony integration skipped or run on single batch.")
    if (!("harmony" %in% Reductions(obj))) {
      obj[["harmony"]] <- obj[["pca"]]
      message("'harmony' reduction set to be same as 'pca'.")
    }
}

module2_harmony_path <- file.path(base_output_path, "module2_harmony.RDS")
tryCatch({
  saveRDS(obj, module2_harmony_path)
  message("Module 2.1 (Harmony) output saved to: ", module2_harmony_path)
}, error = function(e) {
  warning(paste("Error saving Module 2.1 output:", e$message))
})
```

### 4.3. Module 2.2: Cell2Sentence (R & Python)

#### R Environment (Part 1: Preparation)

```{r module2-c2s-preparation, eval=FALSE}
obj_harmony_path <- file.path(base_output_path, "module2_harmony.RDS")
if (file.exists(obj_harmony_path)) {
  obj <- tryCatch(readRDS(obj_harmony_path), error = function(e) {
    message(paste("Error loading Module 2.1 output:", e$message))
    NULL
  })
  if (is.null(obj)) stop("Failed to load data for Module 2.2")
  reduction_for_clustering <- "harmony"
  message("Loaded Harmony corrected object for Cell2Sentence preparation.")
} else {
  stop("Harmony processed Seurat object (module2_harmony.RDS) not found for Module 2.2.")
}

is_rna_split_c2s <- (is.list(obj@assays$RNA) && inherits(obj@assays$RNA[[1]], "Assay")) ||
                (length(obj@assays$RNA@layers) > 0 && !is.null(names(obj@assays$RNA@layers)))
if (is_rna_split_c2s) {
    message("Joining RNA assay layers before clustering for Cell2Sentence.")
    obj <- JoinLayers(obj, assay = "RNA")
}
DefaultAssay(obj) <- "RNA"
obj <- FindNeighbors(obj, dims = 1:min(50, ncol(obj[[reduction_for_clustering]])), reduction = reduction_for_clustering)
obj <- FindClusters(obj, resolution = 0.1, cluster.name = "rna_snn_res.0.1")
obj <- FindClusters(obj, resolution = 2, cluster.name = "rna_snn_res.2")
obj$cell_type_low_res <- obj$rna_snn_res.0.1

h5seurat_path <- file.path(base_output_path, "module2_for_c2s.h5seurat")
h5ad_path <- file.path(base_output_path, "module2_for_c2s.h5ad")
SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE, assay = "RNA")
Convert(h5seurat_path, dest = "h5ad", assay = "RNA", overwrite = TRUE)
message("H5Seurat and H5AD files saved for Cell2Sentence.")

Sys.setenv(H5AD_FILE_PATH = h5ad_path)
c2s_output_dir_env <- file.path(base_output_path, "C2S_output")
Sys.setenv(C2S_OUTPUT_DIR = c2s_output_dir_env)
Sys.setenv(C2S_EMBEDDINGS_CSV = file.path(c2s_output_dir_env, "c2s_cell_embeddings.csv"))
Sys.setenv(C2S_PREDICTED_CSV = file.path(c2s_output_dir_env, "predicted_cell_types.csv"))
Sys.setenv(SPECIES_FOR_C2S = species)

if (!dir.exists(c2s_output_dir_env)) dir.create(c2s_output_dir_env, recursive = TRUE)
```

#### Python Environment (Cell2Sentence Processing)

```{python module2-c2s-python, eval=FALSE}
# Save this as run_cell2sentence.py and execute via reticulate::py_run_file()
import os
current_species = os.getenv("SPECIES_FOR_C2S", "human")
print(f"Cell2Sentence running for species: {current_species}")
# When setting adata.obs["organism"]:
# adata.obs["organism"] = adata.obs.get("organism", current_species)
# ... (rest of your Cell2Sentence Python code)
```

#### R Environment (Part 2: Load C2S Results)

```{r module2-c2s-load-results, eval=FALSE}
# Ensure run_cell2sentence.py has been executed
# reticulate::py_run_file("run_cell2sentence.py")

c2s_embeddings_csv_path <- Sys.getenv("C2S_EMBEDDINGS_CSV")
c2s_predicted_csv_path <- Sys.getenv("C2S_PREDICTED_CSV")

if (!file.exists(c2s_embeddings_csv_path) || !file.exists(c2s_predicted_csv_path)) {
    warning("Cell2Sentence output files not found. Skipping loading C2S results.")
} else {
    embedding_df <- read.csv(c2s_embeddings_csv_path, row.names = 1)
    embedding_matrix <- as.matrix(embedding_df)
    if (!all(rownames(embedding_matrix) %in% colnames(obj))) {
        stop("Mismatch in cell IDs between C2S embeddings and Seurat object.")
    }
    embedding_matrix <- embedding_matrix[colnames(obj), , drop = FALSE]
    colnames(embedding_matrix) <- paste0("C2S_", seq_len(ncol(embedding_matrix)))
    c2s_reduction <- CreateDimReducObject(embeddings = embedding_matrix, key = "C2S_", assay = DefaultAssay(obj))
    obj[["c2s_embeddings"]] <- c2s_reduction
    num_c2s_dims <- ncol(embedding_matrix)
    obj <- RunUMAP(obj, reduction = "c2s_embeddings", dims = 1:num_c2s_dims, reduction.name = "umap_c2s", min.dist = 0.1)
}

module2_c2s_path <- file.path(base_output_path, "module2_c2s_processed.RDS")
tryCatch({
  saveRDS(obj, module2_c2s_path)
  message("Module 2.2 (Cell2Sentence) output saved to: ", module2_c2s_path)
}, error = function(e) {
  warning(paste("Error saving Module 2.2 output:", e$message))
})
```

### 4.4. Module 3: Cell Type Annotation

```{r module3-cell-type-annotation, eval=FALSE}
obj_module2_path <- file.path(base_output_path, "module2_c2s_processed.RDS")
if (!file.exists(obj_module2_path)) {
    obj_module2_path <- file.path(base_output_path, "module2_harmony.RDS")
}
if (file.exists(obj_module2_path)) {
  obj <- tryCatch(readRDS(obj_module2_path), error = function(e) {
    message(e)
    NULL
  })
  if (is.null(obj)) stop("Failed to load data for Module 3.")
} else {
  stop("Processed Seurat object from Module 2 not found.")
}

# Module 3.1: Seurat Clustering (species-agnostic, code omitted for brevity)
# Module 3.2: Annotation with Cell2Sentence Prediction (species-agnostic, code omitted)

# Module 3.3: Annotation with SingleR Prediction
ref_data_celldex <- NULL
celldex_ref_name <- ""

if (species == "human") {
  celldex_ref_name <- "HumanPrimaryCellAtlasData"
  ref_data_celldex <- tryCatch({
    celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
  }, error = function(e) {
    message(paste("Failed to download", celldex_ref_name, "with symbols:", e$message, "\nTrying with ENSEMBL IDs."))
    hpca_ensembl <- tryCatch(celldex::HumanPrimaryCellAtlasData(ensembl = TRUE), error = function(e2) {
      message(e2)
      NULL
    })
    if (is.null(hpca_ensembl)) return(NULL)
    ensembl_ids <- rownames(hpca_ensembl)
    gene_symbols <- mapIds(org_db_object, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    valid_map_indices <- !is.na(gene_symbols)
    hpca_filtered <- hpca_ensembl[valid_map_indices, ]
    unique_gene_symbols <- make.unique(gene_symbols[valid_map_indices])
    rownames(hpca_filtered) <- unique_gene_symbols
    hpca_filtered <- hpca_filtered[!duplicated(rownames(hpca_filtered)), ]
    return(hpca_filtered)
  })
} else if (species == "mouse") {
  celldex_ref_name <- "MouseRNAseqData"
  ref_data_celldex <- tryCatch({
    celldex::MouseRNAseqData(ensembl = FALSE)
  }, error = function(e) {
    message(paste("Failed to download", celldex_ref_name, "with symbols:", e$message, "\nTrying with ENSEMBL IDs."))
    mrna_ensembl <- tryCatch(celldex::MouseRNAseqData(ensembl = TRUE), error = function(e2) {
      message(e2)
      NULL
    })
    if (is.null(mrna_ensembl)) return(NULL)
    ensembl_ids <- rownames(mrna_ensembl)
    gene_symbols <- mapIds(org_db_object, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    valid_map_indices <- !is.na(gene_symbols)
    mrna_filtered <- mrna_ensembl[valid_map_indices, ]
    unique_gene_symbols <- make.unique(gene_symbols[valid_map_indices])
    rownames(mrna_filtered) <- unique_gene_symbols
    mrna_filtered <- mrna_filtered[!duplicated(rownames(mrna_filtered)), ]
    return(mrna_filtered)
  })
} else {
  stop("SingleR reference data not configured for the specified species.")
}

if (is.null(ref_data_celldex) || nrow(ref_data_celldex) == 0) {
    stop(paste("Failed to load or process", celldex_ref_name, "for SingleR."))
}
message(paste("Using", celldex_ref_name, "for SingleR."))
message("Reference cell type distribution (Celldex):")
print(table(ref_data_celldex$label.main))

if (!("logcounts" %in% SummarizedExperiment::assayNames(as.SingleCellExperiment(obj, assay = "RNA")))) {
    if ("data" %in% slotNames(obj@assays$RNA)) {
        message("Logcounts not found, attempting to use 'data' slot from RNA assay for SingleR.")
    } else {
        stop("Suitable log-normalized counts not found for SingleR.")
    }
}

obj_sce <- as.SingleCellExperiment(obj, assay = "RNA")
predictions_singler <- tryCatch({
  SingleR(
    test = obj_sce,
    assay.type.test = "logcounts",
    ref = ref_data_celldex,
    labels = ref_data_celldex$label.main
  )
}, error = function(e) {
  message(paste("SingleR annotation failed:", e$message))
  NULL
})

if (!is.null(predictions_singler)) {
  message("SingleR predicted cell type distribution:")
  print(table(predictions_singler$labels))
  obj$singler_cell_type_main <- predictions_singler$labels
} else {
  obj$singler_cell_type_main <- NA
  warning("SingleR annotation was not successful.")
}

module3_final_annotated_path <- file.path(base_output_path, "module3_final_annotated.RDS")
tryCatch({
  saveRDS(obj, module3_final_annotated_path)
  message("Module 3 (All annotations) output saved to: ", module3_final_annotated_path)
}, error = function(e) {
  warning(paste("Error saving Module 3 output:", e$message))
})
```

### 4.5. Module 4: Downstream Analysis

#### 4.5.1. Load Data and Set Cell Types

```{r module4-load-data, eval=FALSE}
obj_module3_path <- file.path(base_output_path, "module3_final_annotated.RDS")
if (file.exists(obj_module3_path)) {
  data <- tryCatch(readRDS(obj_module3_path), error = function(e) {
    message(paste("Error loading Module 3 output:", e$message))
    NULL
  })
  if (is.null(data)) stop("Failed to load data for Module 4.")
  message("Loaded final annotated object for Module 4.")
} else {
  stop("Annotated Seurat object from Module 3 not found.")
}

if ("singler_cell_type_main" %in% colnames(data@meta.data) && sum(!is.na(data$singler_cell_type_main)) > 0) {
  data$cell_type <- data$singler_cell_type_main
  message("Using 'singler_cell_type_main' as 'cell_type'.")
} else if ("c2s_predicted_cell_type" %in% colnames(data@meta.data) && sum(!is.na(data$c2s_predicted_cell_type)) > 0) {
  data$cell_type <- data$c2s_predicted_cell_type
  message("Using 'c2s_predicted_cell_type' as 'cell_type'.")
} else if ("cell_type_seurat_clustering" %in% colnames(data@meta.data)) {
  data$cell_type <- data$cell_type_seurat_clustering
  message("Using 'cell_type_seurat_clustering' as 'cell_type'.")
} else {
  warning("No suitable 'cell_type' column found. Using 'seurat_clusters' if available.")
  if ("seurat_clusters" %in% colnames(data@meta.data)) data$cell_type <- data$seurat_clusters else stop("No cell type column available.")
}
data$cell_type <- as.factor(data$cell_type)
if (is.list(data@assays$RNA) || (length(data@assays$RNA@layers) > 0 && !is.null(names(data@assays$RNA@layers)))) {
    data <- JoinLayers(data, assay = "RNA")
}
DefaultAssay(data) <- "RNA"
```

#### 4.5.2. Basic Visualization

```{r module4-visualization, eval=FALSE}
# Visualization code using data$cell_type
# Ensure functions like DimPlot2 are available if used
```

#### 4.5.3. Differential Expression & GSEA

```{r module4-dge-gsea, eval=FALSE}
# DGE is species-agnostic
message("GSEA with SeuratExtend is commented out; use clusterProfiler with species-specific gene sets (e.g., GO.db, KEGG.db, msigdbr).")
# Example: options(spe = species)
# if (exists("GeneSetAnalysisGO") && exists("GO_Data")) {
#   go_data_species_specific <- if (species == "human") GO_Data$human else if (species == "mouse") GO_Data$mouse else NULL
#   if (!is.null(go_data_species_specific)) { ... }
# }
```

#### 4.5.4. Pathway and TF Activity Scores with DecoupleR

```{r module4-decoupler, eval=FALSE}
decoupler_organism <- if (species == "human") "human" else if (species == "mouse") "mouse" else "human"
message(paste("Using organism '", decoupler_organism, "' for DecoupleR.", sep = ""))

if (!("data" %in% slotNames(data@assays$RNA))) {
  data <- NormalizeData(data)
}
mat_norm_decoupler <- GetAssayData(data, assay = "RNA", layer = "data")

net_collectri <- tryCatch(decoupleR::get_collectri(organism = decoupler_organism, split_complexes = FALSE), error = function(e) {
  message(paste("Error getting CollecTRI:", e$message))
  NULL
})
if (!is.null(net_collectri)) {
  tf_acts_ulm <- decoupleR::run_ulm(as.matrix(mat_norm_decoupler), net = net_collectri, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
  tf_scores_wide <- tf_acts_ulm %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source')
  data[['tfs_ulm']] <- CreateAssayObject(counts = tf_scores_wide)
  DefaultAssay(data) <- "tfs_ulm"
  data <- ScaleData(data)
} else {
  message("Skipping TF activity due to network loading failure.")
}

net_progeny <- tryCatch(decoupleR::get_progeny(organism = decoupler_organism, top = 500), error = function(e) {
  message(paste("Error getting PROGENy:", e$message))
  NULL
})
if (!is.null(net_progeny)) {
  pathway_acts_mlm <- decoupleR::run_mlm(as.matrix(mat_norm_decoupler), net = net_progeny, .source = 'source', .target = 'target', .mor = 'weight', minsize = 5)
  pathway_scores_wide <- pathway_acts_mlm %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source')
  data[['pathways_mlm']] <- CreateAssayObject(counts = pathway_scores_wide)
  DefaultAssay(data) <- "pathways_mlm"
  data <- ScaleData(data)
} else {
  message("Skipping PROGENy pathway activity due to network loading failure.")
}
DefaultAssay(data) <- "RNA"
```

#### 4.5.5. Gene Scores with UCell

```{r module4-ucell, eval=FALSE}
msigdbr_species_param <- if (species == "human") "Homo sapiens" else if (species == "mouse") "Mus musculus" else "Homo sapiens"
message(paste("Using species '", msigdbr_species_param, "' for MSigDB (UCell).", sep = ""))
hallmark_sets_ucell <- tryCatch(msigdbr(species = msigdbr_species_param, category = "H"), error = function(e) {
  message(paste("Error getting MSigDB Hallmark sets:", e$message))
  NULL
})
if (!is.null(hallmark_sets_ucell) && nrow(hallmark_sets_ucell) > 0) {
  hallmark_list_ucell <- hallmark_sets_ucell %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(gene_symbols = list(gene_symbol), .groups = "drop") %>%
    tibble::deframe()
  data <- AddModuleScore_UCell(data, features = hallmark_list_ucell, assay = "RNA", slot = "data", name = "_UCell")
} else {
  message("Failed to get Hallmark gene sets from MSigDB. Skipping UCell.")
}
```

#### 4.5.6. Predict CNV Activity with InferCNV

```{r module4-infercnv, eval=FALSE}
ref_group_names_infercnv <- "Normal_Epithelial" # MODIFY THIS
message("For InferCNV, ensure 'ref_group_names_infercnv' is correctly set: ", paste(ref_group_names_infercnv, collapse = ", "))

gene_order_file_path_m4.5 <- if (species == "human") {
  "hg38_gencode_v27.txt"
} else if (species == "mouse") {
  "mm10_gencode_vM25.txt"
} else {
  stop("Species not supported for InferCNV gene order file example.")
}
message(paste("IMPORTANT: Gene order file for InferCNV is '", gene_order_file_path_m4.5, "'. Ensure it is correct for species '", species, "'.", sep = ""))

if (!file.exists(gene_order_file_path_m4.5)) {
  warning(paste("Gene order file '", gene_order_file_path_m4.5, "' not found. InferCNV may not run correctly."))
}

if (!(any(ref_group_names_infercnv %in% unique(data$cell_type)))) {
  warning(paste("Reference group(s) for InferCNV ('", paste(ref_group_names_infercnv, collapse = ","), "') not found. Skipping InferCNV.", sep = ""))
} else if (!file.exists(gene_order_file_path_m4.5)) {
  warning("Gene order file not found for InferCNV. Skipping InferCNV.")
} else {
  infercnv_out_dir_m4.5 <- file.path(base_output_path, "InferCNV_output_m4.5/")
  if (!dir.exists(infercnv_out_dir_m4.5)) dir.create(infercnv_out_dir_m4.5, recursive = TRUE)

  counts_matrix_infercnv_m4.5 <- GetAssayData(data, assay = "RNA", layer = "counts")
  counts_file_path_m4.5 <- file.path(infercnv_out_dir_m4.5, "counts_matrix.txt")
  tryCatch(write.table(as.matrix(counts_matrix_infercnv_m4.5), file = counts_file_path_m4.5, quote = FALSE, sep = "\t", col.names = NA),
           error = function(e) warning("Error writing InferCNV counts matrix."))

  annotations_df_infercnv_m4.5 <- data.frame(CellID = colnames(data), CellType = data$cell_type)
  annotations_file_path_m4.5 <- file.path(infercnv_out_dir_m4.5, "cellAnnotations.txt")
  tryCatch(write.table(annotations_df_infercnv_m4.5, file = annotations_file_path_m4.5, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE),
           error = function(e) warning("Error writing InferCNV annotations file."))

  infercnv_obj_m4.5 <- tryCatch(CreateInfercnvObject(
    raw_counts_matrix = counts_file_path_m4.5, annotations_file = annotations_file_path_m4.5, delim = "\t",
    gene_order_file = gene_order_file_path_m4.5, ref_group_names = ref_group_names_infercnv
  ), error = function(e) {
    message(paste("Error creating InferCNV object:", e$message))
    NULL
  })

  if (!is.null(infercnv_obj_m4.5)) {
      num_threads_infercnv <- max(1, floor(parallel::detectCores() / 2))
      infercnv_run_obj_m4.5 <- tryCatch({
        infercnv::run(
          infercnv_obj_m4.5, cutoff = 0.1, out_dir = infercnv_out_dir_m4.5,
          cluster_by_groups = TRUE, denoise = TRUE, HMM = TRUE,
          num_threads = num_threads_infercnv, analysis_mode = "cells"
        )
      }, error = function(e) {
        message(paste("InferCNV run failed:", e$message))
        return(NULL)
      })

      if (!is.null(infercnv_run_obj_m4.5)) {
        message("InferCNV analysis completed. Output in: ", infercnv_out_dir_m4.5)
      } else {
        message("InferCNV run was not successful.")
      }
  } else {
    message("Skipping InferCNV run due to object creation failure.")
  }
}
```

#### 4.5.7. Pseudotime Analysis with SeuratExtend (Palantir)

```{r module4-pseudotime, eval=FALSE}
message("Module 4.6 (Pseudotime Analysis with Palantir/SeuratExtend) is commented out...")
# Requires SeuratExtend or alternative tools (Monocle3, Slingshot)
```

#### 4.5.8. Project Query onto Embedding

```{r module4-query-projection, eval=FALSE}
message("Module 4.7 (Project Query): Species handling for query projection is complex. Adjust homologene function (mouse2human or human2mouse) based on your query and reference species ('", species, "').")
# Example: query_species <- "mouse"
# ref_species <- species
# if (query_species == "mouse" && ref_species == "human") {
#   mouse_genes <- rownames(query[["RNA"]])
#   human_map <- homologene::mouse2human(mouse_genes)
# }
```

#### 4.5.9. Annotate Spatial Data with CARD

```{r module4-card, eval=FALSE}
# CARD analysis is species-agnostic once gene names are harmonized
# Ensure 'Acinar_Cell_Carcinoma.RDS' is available
```

#### 4.5.10. Save Final Output

```{r module4-save-final, eval=FALSE}
module4_final_path <- file.path(base_output_path, "module4_fully_processed_final.RDS")
tryCatch({
  saveRDS(data, module4_final_path)
  message("IMPACT-sc Module 4 processing complete. Final object saved to: ", module4_final_path)
}, error = function(e) {
  warning(paste("Error saving final Module 4 output:", e$message))
})

message("End of IMPACT-sc analysis workflow.")
```
