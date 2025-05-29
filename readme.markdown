# IMPACT-sc: Integrated Multi-Pipeline Analysis and Characterization of Single-Cell Data

<p align="center">
  <img src="impact_sc_logo.png" alt="IMPACT-sc Logo" width="800"/>
</p>

IMPACT-sc is a modular pipeline designed for the comprehensive analysis of single-cell RNA sequencing (scRNA-seq) data. It integrates various R and Python scripts to perform steps from data processing and normalization to advanced downstream analyses like cell type annotation, differential expression, trajectory inference, and query data projection. The pipeline is configured through an interactive setup script and orchestrated by a main Python script, ensuring flexibility and reproducibility.

<p align="center">
  <img src="overview_image.png" alt="IMPACT-sc Pipeline Overview" width="800"/>
  <br/>
  <em>Figure 1: Overview of the IMPACT-sc pipeline workflow. </em>
</p>

## Table of Contents

1.  [Overview](#overview)
2.  [Features](#features)
3.  [Prerequisites](#prerequisites)
4.  [Installation](#installation)
5.  [Configuration](#configuration)
6.  [Running the Pipeline](#running-the-pipeline)
7.  [Available Modules](#available-modules)
8.  [Key Environment Variables](#key-environment-variables)
9.  [Directory Structure](#directory-structure)
10. [Troubleshooting & Important Notes](#troubleshooting--important-notes)

## Overview

The IMPACT-sc pipeline consists of three main components:
1.  **Dependency Installer (`install_dependencies.sh`):** A shell script to set up the R and Python environments with all necessary packages.
2.  **Interactive Setup (`interactive_setup.py`):** A Python script that guides the user through a series of questions to generate a `impact_sc_params.json` configuration file.
3.  **Pipeline Orchestrator (`run_impact_sc_pipeline.py`):** The main Python script that reads the `params.json` file, selects the specified modules, sets up the appropriate environment variables, and executes the R and Python scripts for each module in sequence.

## Features

* **Modular Design:** Select and run only the analysis modules you need.
* **Interactive Configuration:** Easy setup of parameters via a guided questionnaire.
* **Reproducibility:** All parameters are saved in a JSON file, and module execution is logged.
* **Environment Management:** Uses Conda for Python environments and provides guidance for system R package installation.
* **Comprehensive Analysis:** Covers a wide range of scRNA-seq analysis tasks including:
    * Data processing, normalization, and QC.
    * Batch correction (Harmony).
    * Cell embedding and semantic interpretation (Cell2Sentence).
    * Automated and manual cell type annotation (Seurat, SingleR, Cell2Sentence).
    * Visualization (UMAP, tSNE, feature plots, dot plots).
    * Differential Gene Expression (DGE) and Gene Set Enrichment Analysis (GSEA).
    * Pathway and transcription factor activity analysis (DecoupleR, PROGENy).
    * Gene signature scoring (UCell).
    * Pseudotime/trajectory analysis (Palantir).
    * Query data projection onto a reference.
* **Species Support:** Primarily designed for 'human' and 'mouse' data.

## Prerequisites

* **Operating System:** Primarily developed and tested on Windows (with Git Bash), but should be adaptable to Linux/macOS with appropriate path adjustments.
* **Git Bash (for Windows users):** Required for running the `install_dependencies.sh` script. **Must be run as Administrator.**
* **Conda / Anaconda / Miniconda:** Required for Python environment management. Ensure `conda` is in your system PATH.
* **R Installation:** A working R installation is required.
    * For Windows, **Rtools** (e.g., Rtools42 for R 4.2.x) must be installed, and its `bin` directories added to the system PATH for package compilation.
* **(Windows Specific) Locale Settings:**
    * Go to Control Panel > Clock and Region > Region > Administrative tab.
    * Under 'Language for non-Unicode programs', click 'Change system locale...'
    * Ensure **'Beta: Use Unicode UTF-8 for worldwide language support' is UNCHECKED.**
    * Set the system locale to **'English (United States)'** for best compatibility.
    * **RESTART** your computer after making this change. Failure to do so can cause R package installations to fail.

## Installation

1.  **Clone the Repository:**
    ```bash
    git clone <your-repository-url>
    cd <repository-name>
    ```

2.  **Run the Dependency Installer:**
    This script will attempt to install R packages into your system R library and create a Conda environment named `impact_sc` for Python packages.
    * **On Windows:** Open Git Bash **as Administrator**.
    * Navigate to the cloned repository directory.
    * Modify `USER_R_BIN_PATH` in `install_dependencies.sh` if your R `bin` directory is not `F:/R-4.2.3/bin` or not in the system PATH.
    * Execute the script:
        ```bash
        bash install_dependencies.sh
        ```
    * This process can take a very long time. Monitor the output and check the log files:
        * `r_package_install_system_r_only.log` for R package installation details.
        * `python_env_install.log` for Python environment and package installation details.
    * **Address any errors reported in the logs.** Common issues include permission errors (if not run as admin), missing Rtools, or incorrect locale settings.

## Configuration

Once dependencies are installed, you need to generate a parameters file for your specific analysis run.

1.  **Activate the Conda environment (if not already active for running the setup):**
    ```bash
    conda activate impact_sc
    ```
    This is recommended in the `install_dependencies.sh` and `interactive_setup.py` output.

2.  **Run the Interactive Setup Script:**
    ```bash
    python interactive_setup.py
    ```
    This script will ask you a series of questions about:
    * Paths to your R and Python module scripts (e.g., the `scripts.AI` directory).
    * Path to your `Rscript` executable.
    * Input data file(s) and species.
    * Desired output directory.
    * Modules you wish to run.
    * Specific parameters required for selected modules (e.g., paths to reference files, gene lists).

3.  **Parameters File:**
    The script will save your answers into a JSON file named `impact_sc_params.json` inside your specified output directory. This file will be used by the pipeline orchestrator.

## Running the Pipeline

1.  **Activate the Conda Environment:**
    ```bash
    conda activate impact_sc
    ```
    This is specified in the `run_impact_sc_pipeline.py` script's warning and general Conda best practices.

2.  **Execute the Pipeline Orchestrator:**
    Provide the path to the `impact_sc_params.json` file generated during the configuration step.
    ```bash
    python run_impact_sc_pipeline.py /path/to/your/output_directory/impact_sc_params.json
    ```
    * The orchestrator will execute the selected modules in sequence.
    * Log files for each module (e.g., `01_data_processing_log.txt`) will be saved in the output directory.
    * If a module fails, the pipeline will stop. Check the corresponding log file for error messages.

## Available Modules

The following modules can be selected during the interactive setup. The orchestrator (`run_impact_sc_pipeline.py`) determines the script (R or Python) to run for each module.

* **`01_data_processing`** (R): Initial data loading, QC, filtering, and normalization.
    * _Input:_ Raw scRNA-seq data path(s) from `params:input_data_paths`.
* **`02a_harmony_c2s_prep`** (R): Prepares data for Harmony batch correction and subsequent Cell2Sentence analysis (e.g., generates an H5AD file).
* **`02b_c2s`** (Python): Runs Cell2Sentence for semantic embedding and cell type prediction.
    * _Requires:_
        * `params:h5ad_path_for_c2s`: Path to H5AD input file (typically output of module 02a).
        * `params:c2s_model_path_or_name`: Path to a local Cell2Sentence model or a Hugging Face model name.
* **`02c_load_c2s_result`** (R): Loads Cell2Sentence results back into the Seurat object.
* **`03_cell_type_annotation`** (R): Performs cell type annotation using Seurat clustering, SingleR, and integrates Cell2Sentence predictions.
    * _Requires if module selected:_
        * `params:local_singler_ref_path`: Path to a local SingleR reference RDS file.
    * _Optional:_
        * `params:final_cell_type_source`: Determines which annotation source (`auto`, `seurat`, `c2s`, `singler`) is used for the final `cell_type` column.
* **`04a_basic_visulization`** (R): Generates basic visualizations like UMAPs, tSNEs, feature plots, and dot plots.
    * _Optional:_
        * `params:featureplot_genes`: Comma-separated list of genes for FeaturePlots.
        * `params:dotplot_gene_groups`: List of gene groups for DotPlots.
* **`04b_DE_gsea`** (R): Performs Differential Gene Expression analysis and Gene Set Enrichment Analysis.
    * _Optional:_
        * `params:de_gsea_plot_gene`: A single gene name for generating an expression violin plot.
* **`04c_decoupler`** (R): Infers pathway and transcription factor activities using DecoupleR.
    * _Requires if module selected:_
        * `params:collectri_csv_path`: Path to CollecTRI network CSV file for the specified species.
    * _Optional (for human species):_
        * `params:progeny_csv_path`: Path to PROGENy pathway activity network CSV file.
* **`04d_ucell_scores`** (R): Calculates gene signature scores using UCell.
    * _Optional:_
        * `params:msigdb_category`: MSigDB category for UCell gene sets (e.g., "H", "C2").
        * `params:ucell_plot_pathway_name`: Specific pathway/gene set name to plot.
* **`04e_pseudotime`** (R): Performs pseudotime/trajectory analysis using Palantir. (Note: `interactive_setup.py` collects `palantir_start_cell` under a key that seems to align with this module functionality).
    * _Optional:_
        * `params:conditional_paths:palantir_start_cell`: Barcode of a known progenitor/start cell.
* **`04f_query_projection`** (R): Projects a query dataset onto the processed reference dataset. (Note: `interactive_setup.py` collects query RDS info under a key that seems to align with this module functionality).
    * _Requires if module selected:_
        * `params:conditional_paths:query_rds_path`: Path to the query Seurat object RDS file.
        * `params:conditional_paths:query_species`: Species of the query dataset.

## Key Environment Variables

The `run_impact_sc_pipeline.py` script sets various environment variables for the individual R and Python scripts based on the `impact_sc_params.json` file. Some important ones include:

* `IMPACT_SC_SPECIES`: e.g., "human", "mouse"
* `IMPACT_SC_OUTPUT_DIR`: Path to the main output directory.
* `IMPACT_SC_INPUT_DATA_PATH`: Path to the primary input data file (first from the list).
* `IMPACT_SC_INPUT_DATA_PATHS`: Semicolon-separated list of all input data paths.
* `H5AD_FILE_PATH` (for module `02b_c2s`): Path to the H5AD input file.
* `C2S_MODEL_PATH_OR_NAME` (for module `02b_c2s`): Model name or path for Cell2Sentence.
* `C2S_OUTPUT_DIR`, `C2S_EMBEDDINGS_CSV`, `C2S_PREDICTED_CSV` (for module `02b_c2s`).
* `IMPACT_SC_LOCAL_SINGLER_REF_PATH` (for module `03_cell_type_annotation`).
* `IMPACT_SC_FINAL_CELL_TYPE_SOURCE` (for module `03_cell_type_annotation`).
* `IMPACT_SC_FEATUREPLOT_GENES`, `IMPACT_SC_DOTPLOT_GENES_JSON` (for module `04a_basic_visulization`).
* `IMPACT_SC_DE_GENE` (for module `04b_DE_gsea`).
* `IMPACT_SC_COLLECTRI_CSV_PATH`, `IMPACT_SC_PROGENY_CSV_PATH` (for module `04c_decoupler`).
* `IMPACT_SC_MSIGDB_CATEGORY`, `IMPACT_SC_UCELL_PLOT_PATHWAY_NAME` (for module `04d_ucell_scores`).
* `IMPACT_SC_PALANTIR_START_CELL` (for module `04e_pseudotime` or similar, based on `04e_pseudotime` in orchestrator).
* `IMPACT_SC_QUERY_RDS_PATH`, `IMPACT_SC_QUERY_SPECIES` (for module `04f_query_projection` or similar, based on `04f_query_projection` in orchestrator).

## Directory Structure (Example)
