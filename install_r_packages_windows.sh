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
# Removed JAGS/infercnv specific warnings from here
echo "WARNING: This method is generally not recommended due to potential conflicts."
echo "Ensure your system R is correctly configured in your PATH."

# --- Attempt to add user-specified R path to the script's PATH ---
# This is useful if R is not in the system PATH.
# The path F:\R-4.2.3\bin becomes /f/R-4.2.3/bin in Git Bash on Windows.
USER_R_BIN_PATH="/f/R-4.2.3/bin" #
if [ -d "$USER_R_BIN_PATH" ]; then
    echo "Attempting to add $USER_R_BIN_PATH to PATH for this script session."
    export PATH="$USER_R_BIN_PATH:$PATH"
    echo "Current PATH (first few entries): $(echo $PATH | cut -d':' -f1-5)..."
else
    echo "Specified R bin path $USER_R_BIN_PATH not found, relying on system PATH."
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
    echo "ERROR: Rscript command not found in PATH even after attempting to add user-specified path. Please install R and add its bin directory to your system PATH, or verify the USER_R_BIN_PATH in this script."
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
    cat("--- R Environment Details (System R) ---\\n")
    cat("R.version.string:", R.version.string, "\\n")
    cat("Sys.getenv('PATH'):\\n")
    path_sep <- if (.Platform\$OS.type == "windows") ";" else ":"
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

    install_and_check_bioc <- function(pkg_name, lib_path) {
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
        cat("ERROR during default BiocManager::install for package:", pkg_name, "\\nError: ", e$message, "\\nWill try source.\\n")
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
            tryCatch(library(pkg_name, lib.loc = lib_path), error = function(e_lib) {
                cat("Error during library(", pkg_name, "):\\n")
                print(e_lib)
            })
          }
        }, error = function(e_src) {
          cat("ERROR during SOURCE BiocManager::install for package:", pkg_name, "\\nError: ", e_src$message, "\\n")
        })
      }
      return(installed_ok)
    }

    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    cat("CRAN mirror set to https://cloud.r-project.org/\\n")

    # --- Pre-emptive spatstat.utils handling ---
    cat("--- Pre-emptive spatstat.utils handling ---\\n")
    cat("Attempting to ensure spatstat.utils is version >= 3.1.0 by installing from source.\\n")
    cat("This may involve removing an older version if present and R has permission.\\n")
    tryCatch({
      if (requireNamespace("spatstat.utils", quietly = TRUE, lib.loc = user_lib_path)) {
        current_spatstat_version <- packageVersion("spatstat.utils", lib.loc=user_lib_path)
        cat("spatstat.utils is currently installed. Version:", as.character(current_spatstat_version), "\\n")
        if (current_spatstat_version < "3.1.0") {
          cat("Attempting to unload and remove old spatstat.utils...\\n")
          try(detach("package:spatstat.utils", unload=TRUE, character.only=TRUE, force=TRUE), silent=TRUE)
          remove.packages("spatstat.utils", lib = user_lib_path)
          cat("Old spatstat.utils removed (or attempt made). Now installing from source.\\n")
          install.packages("spatstat.utils", Ncpus = max(1, parallel::detectCores() - 1), type = "source", lib = user_lib_path)
        } else {
          cat("spatstat.utils is already version >= 3.1.0.\\n")
        }
      } else {
        cat("spatstat.utils not found. Installing from source...\\n")
        install.packages("spatstat.utils", Ncpus = max(1, parallel::detectCores() - 1), type = "source", lib = user_lib_path)
      }
    }, error = function(e) {
      cat("Error during pre-emptive spatstat.utils handling: ", e$message, "\\n")
      cat("If spatstat.utils is not correctly installed/updated, downstream packages might fail.\\n")
    })
    cat("Finished pre-emptive spatstat.utils handling. Current version (if installed):\\n")
    tryCatch({
        if (requireNamespace("spatstat.utils", quietly = TRUE, lib.loc = user_lib_path)) {
            print(packageVersion("spatstat.utils", lib.loc=user_lib_path))
        } else {
            cat("spatstat.utils still not found after handling attempt.\\n")
        }
    }, error = function(e_pv) {
        cat("Error checking spatstat.utils version after install attempt: ", e_pv$message, "\\n")
    })
    cat("---------------------------------------------\\n")

    cat("Attempting to update other installed packages from CRAN in library:", user_lib_path, "...\\n")
    update.packages(ask = FALSE, checkBuilt = TRUE, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)

    cat("Installing BiocManager if not present...\\n")
    if (!requireNamespace("BiocManager", quietly = TRUE, lib.loc = user_lib_path)) {
        install.packages("BiocManager", Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)
    }
    library(BiocManager, lib.loc = user_lib_path)

    r_version_major_minor <- paste0(R.version\$major, ".", substr(R.version\$minor, 1, 1))
    BIOC_VERSION_FOR_R <- "3.16" # For R 4.2.x
    cat(paste0("Targeting Bioconductor version ", BIOC_VERSION_FOR_R, " for R ", R.version.string, "\\n"))
    tryCatch(
        BiocManager::install(version = BIOC_VERSION_FOR_R, update = FALSE, ask = FALSE, lib = user_lib_path),
        error = function(e) message(paste0("Error ensuring BiocManager version ", BIOC_VERSION_FOR_R, ". Error: ", e$message))
    )
    tryCatch(
        BiocManager::valid(lib.loc = user_lib_path),
        warning = function(w) {cat("BiocManager::valid() warnings:\\n"); print(w)},
        error = function(e) {cat("BiocManager::valid() error:\\n"); print(e)}
    )

    cat("Installing core Bioconductor annotation packages (will attempt source if binary fails)...\\n")
    core_bioc_annotation_pkgs <- c("GenomeInfoDbData", "GenomeInfoDb", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db")
    for (pkg in core_bioc_annotation_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
        install_and_check_bioc(pkg, lib_path = user_lib_path)
      } else {
        cat("Package", pkg, "already installed.\\n")
      }
    }

    cat("Attempting to install/update critical CRAN packages (fastmap, rlang, cli, htmltools, digest)...\\n")
    critical_cran_updates <- c("fastmap", "rlang", "cli", "htmltools", "digest")
    install.packages(critical_cran_updates, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)

    cat("Installing other core CRAN packages and devtools dependencies...\\n")
    cran_core_and_pre_deps_binary <- c(
        "Rcpp", "reticulate", "igraph", "openssl", "curl", "xml2", "V8",
        "sf", "s2", "units", "future", "future.apply", "promises", "later", "httpuv", "shiny",
        "usethis", "pkgload", "pkgbuild", "sessioninfo", "desc", "purrr", "jsonlite", "httr", "remotes",
        "dplyr", "ggplot2", "Matrix", "tibble", "tidyr", "viridis", "reshape2", "pheatmap", "cowplot", "ggpubr", "patchwork",
        "stringr", "car", "ggcorrplot", "homologene"
    )
    install.packages(cran_core_and_pre_deps_binary, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)
    
    # MODIFIED: Install presto from source, TOAST will be installed via BiocManager
    cat("Attempting to install presto from source...\\n")
    install.packages("presto", Ncpus = max(1, parallel::detectCores() - 1), type = "source", lib = user_lib_path)

    # MODIFIED: Install TOAST using BiocManager
    BiocManager::install("TOAST")

    cat("Attempting to install other spatstat family packages (geom, random, core) from CRAN...\\n")
    spatstat_other_pkgs <- c("spatstat.geom", "spatstat.random", "spatstat.core")
    install.packages(spatstat_other_pkgs, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)

    cat("Installing devtools if not present...\\n")
    if (!requireNamespace("devtools", quietly = TRUE, lib.loc = user_lib_path)) {
        install.packages("devtools", Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)
    }
    library(devtools, lib.loc = user_lib_path)
    if (!requireNamespace("remotes", quietly = TRUE, lib.loc = user_lib_path)) {
        install.packages("remotes", Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)
    }
    library(remotes, lib.loc = user_lib_path)

    cat("Installing specific CRAN/GitHub dependencies...\\n")
    tryCatch(remotes::install_github("HenrikBengtsson/R.utils", ref="develop", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for HenrikBengtsson/R.utils: ", e)))
    tryCatch(remotes::install_github("sajuukLyu/ggunchull", upgrade = "never", build_vignettes = FALSE, type = "source", force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for sajuukLyu/ggunchull: ", e)))

    cat("Installing Seurat from CRAN and its other dependencies explicitly...\\n")
    install.packages("Seurat", Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)

    seurat_known_deps_cran <- c(
        'fitdistrplus', 'ggridges', 'ica', 'irlba', 'lmtest', 'matrixStats', 'pbapply', 'plotly', 'RANN',
        'RcppAnnoy', 'ROCR', 'Rtsne', 'scattermore',
        'uwot', 'RcppProgress', 'miniUI', 'htmlwidgets', 'lazyeval', 'gplots', 'ape', 'ggplotify', 'assertthat'
    )
    install.packages(seurat_known_deps_cran, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)

    seurat_known_deps_bioc <- c('leiden', 'sctransform', 'SeuratObject')
    cat("Installing Bioconductor dependencies for Seurat (will attempt source if binary fails)...\\n")
    for (pkg in seurat_known_deps_bioc) {
      if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
        install_and_check_bioc(pkg, lib_path = user_lib_path)
      } else {
        cat("Seurat Bioc dependency", pkg, "already installed or install attempted.\\n")
      }
    }

    cat("Installing Harmony from CRAN...\\n")
    tryCatch(install.packages("harmony", Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path),
             error = function(e) {
                cat("CRAN install failed for harmony: ", e$message, "\\nAttempting source install.\\n")
                tryCatch(install.packages("harmony", Ncpus = max(1, parallel::detectCores() - 1), type = "source", lib = user_lib_path),
                         error = function(e_src) cat("Source install also failed for harmony: ", e_src$message, "\\n"))
             })


    cat("Installing remaining common Bioconductor packages (will attempt source if binary fails)...\\n")
    bioc_packages_remaining <- c(
        "celldex", "decoupleR", "ensembldb", "msigdbr", "scater",
        "scDblFinder", "SCpubr",
        "SingleCellExperiment", "SingleR", "SpatialExperiment", "scran", "UCell"
    )
    for (pkg in bioc_packages_remaining) {
      if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
        install_and_check_bioc(pkg, lib_path = user_lib_path)
      } else {
        cat("Bioconductor package", pkg, "already installed or previously attempted.\\n")
      }
    }

    cat("Installing OmnipathR using BiocManager (will attempt source if binary fails)...\\n")
    install_and_check_bioc('OmnipathR', lib_path = user_lib_path)

    # Removed rjags and infercnv installation sections

    cat("Installing problematic CRAN packages (yulab.utils, ggfun, scatterpie) with fallbacks...\\n")
    problematic_cran_packages <- c("yulab.utils", "ggfun", "scatterpie")
    for (pkg in problematic_cran_packages) {
        if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
            message(paste0("Attempting to install ", pkg, " from CRAN (binary if available) into ", user_lib_path, "..."))
            install.packages(pkg, Ncpus = max(1, parallel::detectCores() - 1), type = "binary", lib = user_lib_path)
            if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
                message(paste0("CRAN binary installation of ", pkg, " failed or package still not found. Attempting to install from GitHub (may require compilation)..."))
                if(!requireNamespace("remotes", quietly=TRUE, lib.loc=user_lib_path)) install.packages("remotes", lib=user_lib_path, type="binary")
                library(remotes, lib.loc=user_lib_path)

                if (pkg == "yulab.utils") {
                    tryCatch(remotes::install_github("YuLab-SMU/yulab.utils", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path), error = function(e) message(paste0("GitHub install failed for yulab.utils: ", e)))
                } else if (pkg == "ggfun") {
                    tryCatch(remotes::install_github("YuLab-SMU/ggfun", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path), error = function(e) message(paste0("GitHub install failed for ggfun: ", e)))
                } else if (pkg == "scatterpie") {
                    tryCatch(remotes::install_github("YuLab-SMU/scatterpie", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path), error = function(e) message(paste0("GitHub install failed for scatterpie: ", e)))
                }
            }
        }
    }

    cat("Installing user-specified GitHub packages (scRNAtoolVis, SeuratExtend)...\\n")
    cat("Installing scRNAtoolVis from GitHub (junjunlab/scRNAtoolVis)...\\n")
    tryCatch(devtools::install_github("junjunlab/scRNAtoolVis", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for junjunlab/scRNAtoolVis: ", e$message)))

    cat("Installing SeuratExtend from GitHub (huayc09/SeuratExtend)...\\n")
    tryCatch(remotes::install_github("huayc09/SeuratExtend", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for huayc09/SeuratExtend: ", e$message)))

    if (requireNamespace("Seurat", quietly = TRUE, lib.loc = user_lib_path)){
        cat("Seurat appears to be installed. Proceeding with Seurat-dependent GitHub packages.\\n")
        cat("Installing SeuratDisk from GitHub (mojaveazure/seurat-disk)...\\n")
        tryCatch(remotes::install_github("mojaveazure/seurat-disk", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
                 error = function(e) message(paste0("GitHub install failed for mojaveazure/seurat-disk: ", e)))
    } else {
        cat("Seurat does not appear to be installed. Skipping SeuratDisk.\\n")
    }

    cat("Installing MuSiC from GitHub (xuranw/MuSiC) using devtools...\\n")
    tryCatch(devtools::install_github("xuranw/MuSiC", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for xuranw/MuSiC: ", e)))

    cat("Installing CARD from GitHub (YingMa0107/CARD) using devtools...\\n")
    tryCatch(devtools::install_github("YingMa0107/CARD", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for CARD: ", e)))

    cat("R package installation script finished.\\n")
    # Return a success indicator
    return(TRUE)

}, error = function(e) {
    # This will catch any error that occurs anywhere in the R script
    cat("A critical error occurred during R script execution:\\n")
    # MODIFIED: More robust error message printing
    cat("Original error type:", paste(class(e), collapse=", "), "\\n")
    if (!is.null(e$message)) {
        cat("Original error message (attempting to print as character):\\n")
        tryCatch({
            cat(paste(as.character(e$message), collapse = "\\n"), "\\n")
        }, error = function(e_print_msg) {
            cat("Failed to print original e$message directly using as.character(). Error during printing: ", e_print_msg$message, "\\n")
            cat("Attempting to print structure of original e$message:\\n")
            try(print(str(e$message)), silent = TRUE)
            cat("Attempting to print original e$message using print():\\n")
            try(print(e$message), silent = TRUE)
        })
    } else {
        cat("No e$message component found in the error object.\\n")
    }
    if (!is.null(e$call)) {
        cat("Original error call:", deparse(e$call), "\\n")
    }
    cat("Full original error object structure (attempting str(e)):\\n")
    try(str(e), silent = TRUE)
    # Force a non-zero exit status for the Rscript process
    quit(status = 1, save = "no")
})

# If main_try_catch_result is NULL (due to quit()), or if it's an error object itself.
if (is.null(main_try_catch_result) || inherits(main_try_catch_result, "error")) {
    # This part might not be reached if quit() was effective
    cat("R script execution failed at a high level or was explicitly quit with status 1.\\n")
    if(exists("quit")) quit(status = 1, save = "no") # Ensure exit if not already
}

EOF

echo "Running R package installation script (this may take a very long time)..."
R_INSTALL_LOG="r_package_install_system_r_only.log" # Log file will be overwritten
echo "Full R installation log will be saved to: $(pwd)/$R_INSTALL_LOG"

Rscript "$INSTALL_R_SCRIPT" > "$R_INSTALL_LOG" 2>&1
R_EXIT_CODE=$?

if [ $R_EXIT_CODE -ne 0 ]; then
    echo "ERROR: R package installation script had errors. Exit code: $R_EXIT_CODE."
    echo "Please check the output above and the detailed log in $R_INSTALL_LOG."
    echo "The R script '$INSTALL_R_SCRIPT' was kept for debugging."
    # Do not exit yet, allow Python setup to proceed if user wants
else
    echo "The R script '$INSTALL_R_SCRIPT' was kept for inspection."
    echo "R package installation script completed. Check $R_INSTALL_LOG for details and any warnings/errors."
fi

echo ""
echo "--- R Dependency Installation Attempt Complete (System-Wide Mode) ---"
echo "R Packages have been attempted to be installed into your system R library."
echo "Please carefully review the log: $R_INSTALL_LOG and the console output for any errors."
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

# Check if conda is available
if ! command -v conda &> /dev/null
then
    echo "ERROR: conda command not found. Please install Anaconda or Miniconda and ensure 'conda' is in your PATH."
    echo "Python environment setup skipped."
    exit 1
fi
echo "Found Conda: $(command -v conda)"
conda --version

ENV_NAME="impact_sc"
PYTHON_VERSION="3.9" # Specify Python version
PYTHON_INSTALL_LOG="python_env_install.log"

echo ""
echo "Step 3: Creating/Updating Conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
echo "Full Python installation log will be saved to: $(pwd)/$PYTHON_INSTALL_LOG"
echo "This may take some time..."

# Create the conda environment if it doesn't exist
if ! conda env list | grep -q "$ENV_NAME"; then
    echo "Creating new Conda environment: $ENV_NAME with Python $PYTHON_VERSION"
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y > "$PYTHON_INSTALL_LOG" 2>&1
    CONDA_CREATE_EXIT_CODE=$?
else
    echo "Conda environment '$ENV_NAME' already exists. Skipping creation."
    echo "If you need to reinstall, please remove the environment first: conda env remove -n $ENV_NAME"
    CONDA_CREATE_EXIT_CODE=0 # Treat as success for this step
fi

if [ $CONDA_CREATE_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Failed to create Conda environment '$ENV_NAME'. Exit code: $CONDA_CREATE_EXIT_CODE."
    echo "Please check the detailed log in $PYTHON_INSTALL_LOG."
    exit 1
else
    echo "Conda environment '$ENV_NAME' (Python $PYTHON_VERSION) is ready or was already present."
fi

echo ""
echo "Step 4: Installing Python packages into '$ENV_NAME'..."

# Activate the environment for subsequent commands (best effort for scripts, conda run is more robust)
# For direct installation, using `conda install -n ENV_NAME` or `conda run -n ENV_NAME pip install` is safer.

echo "Installing core packages using conda..."
# Install packages using conda, specifying channels
# scanpy often brings in anndata, numpy, scipy, pandas, matplotlib, seaborn
# Explicitly listing them for clarity and to ensure they are from preferred channels if needed.
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
    echo "ERROR: Failed to install some core Python packages using conda. Exit code: $CONDA_INSTALL_EXIT_CODE."
    echo "Please check the detailed log in $PYTHON_INSTALL_LOG."
    # Consider exiting or allowing pip installs to proceed
else
    echo "Core Python packages installed successfully via conda."
fi

echo "Installing cell2sentence using pip within the conda environment..."
conda run -n "$ENV_NAME" pip install cell2sentence >> "$PYTHON_INSTALL_LOG" 2>&1
PIP_INSTALL_C2S_EXIT_CODE=$?

if [ $PIP_INSTALL_C2S_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Failed to install 'cell2sentence' using pip. Exit code: $PIP_INSTALL_C2S_EXIT_CODE."
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
# Removed JAGS/infercnv specific reminder from here
echo "- For 'infercnv' to work with rjags, JAGS 4.x must be manually installed first AND R must be able to find it for 'rjags' compilation."

exit 0
