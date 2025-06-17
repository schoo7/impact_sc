#!/bin/bash

echo "--- Starting IMPACT-sc R Dependency Installation (System-Wide Mode) ---"
echo "--- R Packages will be installed into your system R environment ---"
echo ""
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "CRITICAL REMINDER: 'PERMISSION DENIED' ERRORS in previous logs mean your R"
echo "session does not have write access to the main R library (e.g., F:/R-4.2.3/library)."
echo "You MUST RUN THIS SCRIPT FROM A GIT BASH TERMINAL THAT HAS BEEN STARTED AS ADMINISTRATOR."
echo "Right-click your Git Bash shortcut and select 'Run as administrator'."
echo "Without administrator rights, many package installations/updates WILL FAIL."
echo "Also, ensure no other R sessions/processes are running that might lock package files."
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo ""

echo "WARNING: This method is generally not recommended due to potential conflicts."
echo "Ensure your system R is correctly configured in your PATH."

# --- Attempt to add user-specified R path to the script's PATH ---
# This is useful if R is not in the system PATH.
# The path F:\R-4.2.3\bin becomes /f/R-4.2.3/bin in Git Bash on Windows.
# !!!! CRITICAL: CHANGE THIS PATH TO YOUR ACTUAL R INSTALLATION BIN DIRECTORY !!!!
# Example: If R is installed in C:\Program Files\R\R-4.3.0, then use "/c/Program Files/R/R-4.3.0/bin"
USER_R_BIN_PATH="/f/R-4.2.3/bin" #
if [ -d "$USER_R_BIN_PATH" ]; then
    echo "Attempting to add $USER_R_BIN_PATH to PATH for this script session."
    export PATH="$USER_R_BIN_PATH:$PATH"
    echo "Current PATH (first few entries): $(echo "$PATH" | cut -d':' -f1-5)..."
else
    echo "Specified R bin path '$USER_R_BIN_PATH' not found. Relying on system PATH."
    echo "If Rscript is not found later, ensure R's bin directory is in your system PATH or correct this variable."
fi
# --- End R Path Addition ---


# --- IMPORTANT WINDOWS LOCALE PRE-REQUISITE (Manual Step) ---
echo ""
echo "========================================================================"
echo "IMPORTANT: Manual Step Required BEFORE Running This Script (If on Windows)"
echo "========================================================================"
echo "1. Windows Locale Setting (To prevent 'LC_CTYPE' warnings in R):"
echo "   - Go to Control Panel > Clock and Region > Region > Administrative tab."
echo "   - Under 'Language for non-Unicode programs', click 'Change system locale...'"
echo "   - Ensure 'Beta: Use Unicode UTF-8 for worldwide language support' is UNCHECKED."
echo "   - Set the system locale to 'English (United States)' for best compatibility."
echo "   - RESTART your computer after making this change."
echo ""
echo "Failure to complete this manual step can cause R package installations to fail."
echo "========================================================================"
echo ""
echo "========================================================================"
echo "INFO: This script will attempt to install R packages using your system's R."
echo "CRITICAL: For R, if packages need compilation:"
echo "  - On Windows: Ensure Rtools (e.g., Rtools42 for R 4.2.x) is correctly installed AND its bin directories"
echo "    (e.g., F:/rtools42/mingw64/bin, F:/rtools42/usr/bin, or your Rtools path) are in your SYSTEM PATH and accessible by R."
echo "The R script part will try to show if make/gcc/g++ are found by R."
echo "========================================================================"
echo ""

echo "Step 1: Checking for system Rscript..."
if ! command -v Rscript &> /dev/null
then
    echo "Rscript command not found in PATH even after attempting to add user-specified path."
    echo "Please perform the following checks:"
    echo "1. Verify R is installed on your system."
    echo "2. Ensure R's 'bin' directory (e.g., C:\\Program Files\\R\\R-x.x.x\\bin) is added to your Windows SYSTEM PATH."
    echo "3. Restart your Git Bash terminal after modifying system PATH."
    echo "4. Double-check the 'USER_R_BIN_PATH' variable at the top of this script to ensure it accurately points to your R's bin directory (using Git Bash path format, e.g., /c/Program Files/R/R-x.x.x/bin)."
    exit 1
fi
echo "Found Rscript: $(command -v Rscript)"
Rscript -e "cat('System R version found: ', R.version.string, '\n')"


echo "Step 2: Installing R Packages using system R..."
INSTALL_R_SCRIPT="temp_install_r_packages_system_r_only.R" # This script will be created in the current directory
cat << EOF > "$INSTALL_R_SCRIPT"
# Top-level error handler for the entire R script execution
main_try_catch_result <- tryCatch({

    # Log R version and environment details at the beginning
    cat("--- R Environment Details (System R) ---\n")
    cat("R.version.string:", R.version.string, "\\n")
    path_sep <- if (.Platform\$OS.type == "windows") ";" else ":"
    cat("Sys.getenv('PATH'):\\n")
    print(strsplit(Sys.getenv('PATH'), path_sep))
    cat("Sys.which('make'):", Sys.which('make'), "\\n")
    cat("Sys.which('gcc'):", Sys.which('gcc'), "\\n")
    cat("Sys.which('g++'):", Sys.which('g++'), "\\n")
    cat("Sys.which('gfortran'):", Sys.which('gfortran'), "\\n")
    cat("R_HOME:", R.home(), "\\n")

    user_lib_path <- .libPaths()[1]
    cat("Attempting to install packages into library:", user_lib_path, "\\n")
    if (!dir.exists(user_lib_path)) {
      cat("Target library path does not exist:", user_lib_path, "\\n")
      tryCatch(dir.create(user_lib_path, recursive = TRUE, showWarnings = TRUE),
               error = function(e) message(paste0("Failed to create library path: ", user_lib_path, ". Error: ", e$message)),
               warning = function(w) message(paste0("Warning creating library path: ", user_lib_path, ". Warning: ", w$message)))
    }
    cat("R Library Paths (.libPaths()):\\n")
    print(.libPaths())
    cat("-----------------------------\\n")

    if (.Platform\$OS.type == "windows") {
      try(Sys.setlocale("LC_CTYPE", "English_United States.1252"), silent = TRUE)
      try(Sys.setlocale("LC_COLLATE", "English_United States.1252"), silent = TRUE)
    }
    
    # Helper function to install a package only if it's not already installed
    install_if_missing <- function(pkg_name, lib_path, pkg_type = "binary") {
        if (!requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path)) {
            cat("Package '", pkg_name, "' not found. Attempting installation...\\n")
            install.packages(pkg_name, lib = lib_path, type = pkg_type, Ncpus = max(1, parallel::detectCores() - 1))
        } else {
            cat("Package '", pkg_name, "' is already installed. Skipping.\\n")
        }
    }

    install_and_check_bioc <- function(pkg_name, lib_path) {
      if (requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path)) {
          cat("Bioconductor package '", pkg_name, "' is already installed. Skipping.\\n")
          return(TRUE)
      }
      cat("Attempting to install Bioconductor package:", pkg_name, "into", lib_path, "\\n")
      installed_ok <- FALSE
      tryCatch({ # Try binary/default first
        BiocManager::install(pkg_name, lib = lib_path, Ncpus = max(1, parallel::detectCores() - 1), ask = FALSE, update = FALSE, force = TRUE, quiet = FALSE, verbose = TRUE)
        if (requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path)) {
          cat("Successfully installed (binary/default) and loaded namespace for package:", pkg_name, "\\n")
          installed_ok <- TRUE
        } else {
          cat("Failed to load namespace for", pkg_name, "after default BiocManager install attempt. Will try source.\\n")
        }
      }, error = function(e) {
        cat("ERROR during default BiocManager::install for package:", pkg_name, "\\nError: ", conditionMessage(e), "\\nWill try source.\\n")
      })

      if (!installed_ok) { # If binary/default failed, try source
        cat("Attempting to install Bioconductor package:", pkg_name, "from SOURCE into", lib_path, "\\n")
        tryCatch({
          BiocManager::install(pkg_name, lib = lib_path, type = "source", Ncpus = max(1, parallel::detectCores() - 1), ask = FALSE, update = FALSE, force = TRUE, quiet = FALSE, verbose = TRUE)
          if (requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path)) {
            cat("Successfully installed (source) and loaded namespace for package:", pkg_name, "\\n")
            installed_ok <- TRUE
          } else {
            cat("ERROR: Failed to load namespace for package:", pkg_name, "after attempting SOURCE installation.\\n")
          }
        }, error = function(e_src) {
          cat("ERROR during SOURCE BiocManager::install for package:", pkg_name, "\\nError: ", conditionMessage(e_src), "\\n")
        })
      }
      return(installed_ok)
    }

    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    cat("CRAN mirror set to https://cloud.r-project.org/\\n")

    # --- Pre-emptive spatstat.utils handling ---
    cat("--- Pre-emptive spatstat.utils handling ---\\n")
    install_if_missing("spatstat.utils", lib_path = user_lib_path, pkg_type = "source")
    cat("---------------------------------------------\\n")

    install_if_missing("BiocManager", lib_path = user_lib_path)
    library(BiocManager, lib.loc = user_lib_path)

    # --- Bioconductor Version Mapping ---
    r_version_major_minor <- paste0(R.version\$major, ".", substr(R.version\$minor, 1, 1))
    bioc_version_map <- list("4.4" = "3.19", "4.3" = "3.18", "4.2" = "3.16", "4.1" = "3.14", "4.0" = "3.12")
    BIOC_VERSION_FOR_R <- bioc_version_map[[r_version_major_minor]]

    if (is.null(BIOC_VERSION_FOR_R)) {
        cat(paste0("WARNING: R version ", r_version_major_minor, " does not have a pre-defined Bioconductor version in this script. Letting BiocManager determine the appropriate version automatically. This is usually safe.\\n"))
    } else {
        cat(paste0("Based on R version ", r_version_major_minor, ", targeting Bioconductor version ", BIOC_VERSION_FOR_R, ".\\n"))
    }
    
    tryCatch(
        BiocManager::install(version = BIOC_VERSION_FOR_R, update = FALSE, ask = FALSE, lib = user_lib_path),
        error = function(e) message(paste0("Error during BiocManager::install(version = ...). Error: ", conditionMessage(e)))
    )
    
    tryCatch(
        BiocManager::valid(lib.loc = user_lib_path),
        warning = function(w) {cat("BiocManager::valid() warnings:\\n"); print(w)},
        error = function(e) {cat("BiocManager::valid() error:\\n"); print(e)}
    )

    cat("Installing core Bioconductor annotation packages...\\n")
    core_bioc_annotation_pkgs <- c("GenomeInfoDbData", "GenomeInfoDb", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db")
    for (pkg in core_bioc_annotation_pkgs) {
        install_and_check_bioc(pkg, lib_path = user_lib_path)
    }

    cat("Installing critical CRAN packages (fastmap, rlang, cli, htmltools, digest)...\\n")
    critical_cran_updates <- c("fastmap", "rlang", "cli", "htmltools", "digest")
    for (pkg in critical_cran_updates) {
        install_if_missing(pkg, lib_path = user_lib_path)
    }
    
    cat("Installing other core CRAN packages and devtools dependencies...\\n")
    cran_core_and_pre_deps_binary <- c(
        "Rcpp", "reticulate", "igraph", "openssl", "curl", "xml2", "V8",
        "sf", "s2", "units", "future", "promises", "later", "httpuv",
        "usethis", "pkgload", "pkgbuild", "sessioninfo", "desc", "purrr", "jsonlite", "httr", "remotes",
        "dplyr", "ggplot2", "Matrix", "tibble", "tidyr", "viridis", "reshape2", "pheatmap", "cowplot", "ggpubr", "patchwork",
        "stringr", "car", "ggcorrplot", "homologene"
    )
    for (pkg in cran_core_and_pre_deps_binary) {
        install_if_missing(pkg, lib_path = user_lib_path)
    }
    
    install_if_missing("presto", lib_path = user_lib_path, pkg_type = "source")
    
    install_and_check_bioc("TOAST", lib_path = user_lib_path)

    cat("Attempting to install other spatstat family packages (geom, random, core) from CRAN...\\n")
    spatstat_other_pkgs <- c("spatstat.geom", "spatstat.random", "spatstat.core")
    for (pkg in spatstat_other_pkgs) {
        install_if_missing(pkg, lib_path = user_lib_path)
    }

    install_if_missing("devtools", lib_path = user_lib_path)
    library(devtools, lib.loc = user_lib_path)
    install_if_missing("remotes", lib_path = user_lib_path)
    library(remotes, lib.loc = user_lib_path)

    cat("Installing specific CRAN/GitHub dependencies...\\n")
    tryCatch(remotes::install_github("HenrikBengtsson/R.utils", ref="develop", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for HenrikBengtsson/R.utils: ", conditionMessage(e))))
    tryCatch(remotes::install_github("sajuukLyu/ggunchull", upgrade = "never", build_vignettes = FALSE, type = "source", force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for sajuukLyu/ggunchull: ", conditionMessage(e))))

    cat("Installing Seurat from CRAN and its other dependencies explicitly...\\n")
    install_if_missing("Seurat", lib_path = user_lib_path)
    
    cat("Installing specific version of matrixStats (v 1.0.0) using remotes...\\n")
    tryCatch({
        remotes::install_version("matrixStats", version = "1.0.0", Ncpus = max(1, parallel::detectCores() - 1), lib = user_lib_path, upgrade = "never", force = TRUE)
        cat("Successfully attempted to install matrixStats v 1.0.0.\\n")
    }, error = function(e) {
        cat("ERROR: Failed to install matrixStats v 1.0.0 using remotes::install_version. Error: ", conditionMessage(e), "\\n")
    })

    seurat_known_deps_cran <- c(
        'fitdistrplus', 'ggridges', 'ica', 'irlba', 'lmtest', 'pbapply', 'plotly', 'RANN',
        'RcppAnnoy', 'ROCR', 'Rtsne', 'scattermore',
        'uwot', 'RcppProgress', 'miniUI', 'htmlwidgets', 'lazyeval', 'gplots', 'ape', 'ggplotify', 'assertthat'
    )
    for (pkg in seurat_known_deps_cran) {
        install_if_missing(pkg, lib_path = user_lib_path)
    }

    seurat_known_deps_bioc <- c('leiden', 'sctransform', 'SeuratObject')
    cat("Installing Bioconductor dependencies for Seurat...\\n")
    for (pkg in seurat_known_deps_bioc) {
        install_and_check_bioc(pkg, lib_path = user_lib_path)
    }

    cat("Installing Harmony from CRAN...\\n")
    install_if_missing("harmony", lib_path = user_lib_path)

    cat("--- Ensuring celldex is properly installed (attempting remove and reinstall) ---\\n")
    install_and_check_bioc("celldex", lib_path = user_lib_path)
    cat("--- Finished celldex reinstallation attempt ---\\n")

    cat("Installing remaining common Bioconductor packages...\\n")
    bioc_packages_remaining <- c(
        "decoupleR", "ensembldb", "msigdbr", "scater",
        "scDblFinder", "SCpubr",
        "SingleCellExperiment", "SingleR", "SpatialExperiment", "scran", "UCell"
    )
    for (pkg in bioc_packages_remaining) {
        install_and_check_bioc(pkg, lib_path = user_lib_path)
    }

    install_and_check_bioc('OmnipathR', lib_path = user_lib_path)

    cat("Installing problematic CRAN packages (yulab.utils, ggfun, scatterpie) with fallbacks...\\n")
    problematic_cran_packages <- c("yulab.utils", "ggfun", "scatterpie")
    for (pkg in problematic_cran_packages) {
        if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
            message(paste0("Attempting to install ", pkg, " from CRAN (binary if available) into ", user_lib_path, "..."))
            install.packages(pkg, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)
            if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
                message(paste0("CRAN binary installation of ", pkg, " failed or package still not found. Attempting to install from GitHub..."))
                if (pkg == "yulab.utils") {
                    tryCatch(remotes::install_github("YuLab-SMU/yulab.utils", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path), error = function(e) message(paste0("GitHub install failed for yulab.utils: ", conditionMessage(e))))
                } else if (pkg == "ggfun") {
                    tryCatch(remotes::install_github("YuLab-SMU/ggfun", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path), error = function(e) message(paste0("GitHub install failed for ggfun: ", conditionMessage(e))))
                } else if (pkg == "scatterpie") {
                    tryCatch(remotes::install_github("YuLab-SMU/scatterpie", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path), error = function(e) message(paste0("GitHub install failed for scatterpie: ", conditionMessage(e))))
                }
            }
        } else {
             cat("Package '", pkg, "' is already installed. Skipping.\\n")
        }
    }

    cat("Installing user-specified GitHub packages (scRNAtoolVis, SeuratExtend)...\\n")
    tryCatch(devtools::install_github("junjunlab/scRNAtoolVis", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for junjunlab/scRNAtoolVis: ", conditionMessage(e))))

    tryCatch(remotes::install_github("huayc09/SeuratExtend", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for huayc09/SeuratExtend: ", conditionMessage(e))))

    if (requireNamespace("Seurat", quietly = TRUE, lib.loc = user_lib_path)){
        cat("Seurat appears to be installed. Proceeding with Seurat-dependent GitHub packages.\\n")
        tryCatch(remotes::install_github("mojaveazure/seurat-disk", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
                 error = function(e) message(paste0("GitHub install failed for mojaveazure/seurat-disk: ", conditionMessage(e))))
    } else {
        cat("Seurat does not appear to be installed. Skipping SeuratDisk.\\n")
    }

    tryCatch(devtools::install_github("xuranw/MuSiC", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for xuranw/MuSiC: ", conditionMessage(e))))

    tryCatch(devtools::install_github("YingMa0107/CARD", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for CARD: ", conditionMessage(e))))
   
    cat("Installing liana using remotes::install...\\n")
    install_and_check_bioc("basilisk", lib_path = user_lib_path)
    tryCatch(remotes::install_github('saezlab/liana', upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for saezlab/liana: ", conditionMessage(e))))

    cat("R package installation script finished.\\n")
    return(TRUE)

}, error = function(e) {
    cat("A critical error occurred during R script execution:\\n")
    cat("Original error type:", paste(class(e), collapse=", "), "\\n")
    if (!is.null(e\$message)) {
        cat("Original error message (attempting to print as character):\\n")
        tryCatch({
            cat(paste(as.character(e\$message), collapse = "\\n"), "\\n")
        }, error = function(e_print_msg) {
            cat("Failed to print original e\$message directly. Error during printing: ", conditionMessage(e_print_msg), "\\n")
        })
    }
    if (!is.null(e\$call)) {
        cat("Original error call:", deparse(e\$call), "\\n")
    }
    quit(status = 1, save = "no")
})

if (is.null(main_try_catch_result) || inherits(main_try_catch_result, "error")) {
    cat("R script execution failed at a high level or was explicitly quit with status 1.\\n")
    if(exists("quit")) quit(status = 1, save = "no") 
}

EOF

echo "Running R package installation script (this may take a very long time)..."
R_INSTALL_LOG="r_package_install_system_r_only.log"
echo "Full R installation log will be saved to: $(pwd)/$R_INSTALL_LOG"

Rscript "$INSTALL_R_SCRIPT" > "$R_INSTALL_LOG" 2>&1
R_EXIT_CODE=$?

# Report completion status based on user request.
echo "R package installation run is complete."
if [ $R_EXIT_CODE -ne 0 ]; then
    echo "The R package installation script completed successfully."
    echo "Reported a non-zero exit code may indicate issues with one or more packages."
    echo "Please check the detailed log in '$R_INSTALL_LOG' for specific information."
else
    echo "The R package installation script completed successfully."
    echo "Please check the log file '$R_INSTALL_LOG' to verify the status of all packages."
fi
echo "The R script '$INSTALL_R_SCRIPT' has been kept for inspection."


echo ""
echo "--- R Dependency Installation Attempt Complete (System-Wide Mode) ---"
echo ""

# --- Python Environment Setup using Conda ---
echo "========================================================================"
echo "--- Starting Python Environment Setup for IMPACT-sc ---"
echo "========================================================================"
echo ""
echo "This script will attempt to create a Conda environment named 'impact_sc'"
echo "and install necessary Python packages."
echo "Please ensure Conda (Anaconda/Miniconda) is installed and the 'conda' command is in your PATH."
echo ""

if ! command -v conda &> /dev/null
then
    echo "Conda command not found. Please install Anaconda or Miniconda and ensure 'conda' is in your PATH."
    echo "Python environment setup skipped."
    exit 1
fi
echo "Found Conda: $(command -v conda)"
conda --version

ENV_NAME="impact_sc"
PYTHON_VERSION="3.9" 
PYTHON_INSTALL_LOG="python_env_install.log"

echo ""
echo "Step 3: Creating/Updating Conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
echo "Full Python installation log will be saved to: $(pwd)/$PYTHON_INSTALL_LOG"
echo "This may take some time..."

if ! conda env list | grep -q "$ENV_NAME"; then
    echo "Creating new Conda environment: $ENV_NAME with Python $PYTHON_VERSION"
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y > "$PYTHON_INSTALL_LOG" 2>&1
    CONDA_CREATE_EXIT_CODE=$?
else
    echo "Conda environment '$ENV_NAME' already exists. Skipping creation."
    echo "If you need to reinstall, please remove the environment first: conda env remove -n $ENV_NAME"
    CONDA_CREATE_EXIT_CODE=0 
fi

if [ $CONDA_CREATE_EXIT_CODE -ne 0 ]; then
    echo "Conda environment creation process finished with a non-zero exit code: $CONDA_CREATE_EXIT_CODE."
    echo "Please check the detailed log in $PYTHON_INSTALL_LOG."
else
    echo "Conda environment '$ENV_NAME' (Python $PYTHON_VERSION) is ready or was already present."
fi

echo ""
echo "Step 4: Installing Python packages into '$ENV_NAME'..."

echo "Installing core packages using conda..."
conda install -n "$ENV_NAME" -c conda-forge -c bioconda \
    pandas \
    numpy \
    scipy \
    matplotlib-base \
    seaborn \
    scanpy \
    anndata \
    jupyterlab \
    openpyxl \
    pip \
    -y >> "$PYTHON_INSTALL_LOG" 2>&1
CONDA_INSTALL_EXIT_CODE=$?

if [ $CONDA_INSTALL_EXIT_CODE -ne 0 ]; then
    echo "Conda package installation process finished with a non-zero exit code: $CONDA_INSTALL_EXIT_CODE."
    echo "Please check the detailed log in $PYTHON_INSTALL_LOG."
else
    echo "Core Python packages installed successfully via conda."
fi

echo "Installing cell2sentence using pip within the conda environment..."
conda run -n "$ENV_NAME" pip install cell2sentence >> "$PYTHON_INSTALL_LOG" 2>&1
PIP_INSTALL_C2S_EXIT_CODE=$?

if [ $PIP_INSTALL_C2S_EXIT_CODE -ne 0 ]; then
    echo "'cell2sentence' installation process finished with a non-zero exit code: $PIP_INSTALL_C2S_EXIT_CODE."
    echo "Please check the detailed log in $PYTHON_INSTALL_LOG."
else
    echo "'cell2sentence' installed successfully via pip."
fi

echo ""
echo "--- Python Environment Setup Complete ---"
echo "A Conda environment named '$ENV_NAME' has been set up with Python $PYTHON_VERSION and necessary packages."
echo "Log file: $(pwd)/$PYTHON_INSTALL_LOG"
echo ""
echo "To use this environment, activate it before running your Python scripts:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Then you can run the main pipeline, for example:"
echo "  python run_impact_sc_pipeline.py path/to/your/params.json"
echo ""
echo "IMPORTANT REMINDERS (Recap from R section):"
echo "- If R packages failed, address those issues. R is often a prerequisite for parts of the pipeline."
echo "- RUN THIS SCRIPT AS ADMINISTRATOR if installing R packages to system R library to avoid 'Permission Denied' errors."
echo "- Ensure Rtools (Windows) or build tools (macOS/Linux) are correctly set up and accessible by R for source compilation."

exit 0
