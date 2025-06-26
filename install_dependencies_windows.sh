#!/bin/bash

# Ollama Model Configuration (if needed for future script enhancements)
OLLAMA_MODEL_NAME_DEFAULT="gemma3:12b-it-qat" # You can change the default Ollama model if needed
OLLAMA_BASE_URL_DEFAULT="http://localhost:11434" # Base URL for Ollama API

# Detect OS
OS_TYPE=""
case "$(uname -s)" in
    Linux*)     OS_TYPE="linux";;
    Darwin*)    OS_TYPE="mac";;
    CYGWIN*|MSYS*|MINGW*) OS_TYPE="windows";;
    *)          OS_TYPE="unknown";;
esac

echo "--- Starting IMPACT-sc R & Python Dependency Installation ($OS_TYPE Mode) ---"
echo "--- R Packages will be installed into your system R environment ---"
echo ""

if [ "$OS_TYPE" == "windows" ]; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "CRITICAL REMINDER (Windows):"
    echo "- To avoid 'Permission Denied' or 'cannot rename file' errors, you MUST RUN THIS SCRIPT FROM A GIT BASH TERMINAL THAT HAS BEEN STARTED AS ADMINISTRATOR."
    echo "- ENSURE NO OTHER R SESSIONS OR PROCESSES ARE RUNNING that might lock package files (e.g., RStudio, other R consoles, background R processes)."
    echo "- If you still encounter issues, try closing all applications and restarting your computer before rerunning the script."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
elif [ "$OS_TYPE" == "mac" ]; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "CRITICAL REMINDER (macOS):"
    echo "- If you encounter 'Permission Denied' errors, try running this script with 'sudo' (e.g., 'sudo bash install_deps.sh')."
    echo "- Ensure no other R sessions or processes are running that might lock package files."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
else
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "CRITICAL REMINDER:"
    echo "- If you encounter 'Permission Denied' errors, try running this script with 'sudo' (e.g., 'sudo bash install_deps.sh')."
    echo "- Ensure no other R sessions or processes are running that might lock package files."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
fi
echo ""

echo "WARNING: This method is generally not recommended due to potential conflicts."
echo "Ensure your system R is correctly configured in your PATH."

# IMPORTANT OS-SPECIFIC PRE-REQUISITES for R Installation (Manual Step)
echo ""
echo "========================================================================"
echo "IMPORTANT: Manual R 4.2.3 Installation (CRITICAL FIRST STEP)"
echo "========================================================================"
echo "This script CANNOT automatically install R itself. You MUST manually download and install R 4.2.3."
echo "Please follow the instructions below to download R 4.2.3:"
echo ""

if [ "$OS_TYPE" == "windows" ]; then
    echo "  For Windows, download R 4.2.3 from CRAN: https://cran.r-project.org/bin/windows/base/old/4.2.3/R-4.2.3-win.exe"
    echo "  Run the downloaded installer (.exe file) and follow the prompts. Ensure R's 'bin' directory is added to your system PATH during installation."
elif [ "$OS_TYPE" == "mac" ]; then
    echo "  For macOS, download R 4.2.3 from CRAN: https://cran.r-project.org/bin/macosx/old/R-4.2.3.pkg"
    echo "  Run the downloaded installer (.pkg file) and follow the prompts. Ensure R's path is correctly set."
elif [ "$OS_TYPE" == "linux" ]; then
    echo "  For Linux, follow the official CRAN instructions for your specific distribution to install R 4.2.3:"
    echo "  https://cran.r-project.org/bin/linux/"
fi
echo ""
echo "After manually installing R 4.2.3, you MUST provide its 'bin' directory path to this script."
echo "This path contains the 'Rscript' command, which this script needs to install R packages."
echo "========================================================================"
echo ""

# Get R bin path from user interactively
USER_R_BIN_PATH=""
echo "Please enter the full path to your installed R 4.2.3 'bin' directory, for example:"
if [ "$OS_TYPE" == "windows" ]; then
    echo "  Windows (Git Bash format): /c/Program Files/R/R-4.2.3/bin or /f/R-4.2.3/bin"
elif [ "$OS_TYPE" == "mac" ]; then
    echo "  macOS: /Library/Frameworks/R.framework/Versions/4.2/Resources/bin or /opt/homebrew/bin"
elif [ "$OS_TYPE" == "linux" ]; then
    echo "  Linux: /usr/bin/R or /usr/local/bin/R"
fi
read -p "R 4.2.3 'bin' directory path: " USER_R_BIN_PATH

# Attempt to add user-specified R path to the script's PATH
if [ -z "$USER_R_BIN_PATH" ]; then
    echo "WARNING: R 'bin' directory path not provided. Script will rely on 'Rscript' being in system PATH."
    echo "If 'Rscript command not found' error occurs later, please ensure R's 'bin' directory is in your system PATH, or re-run the script and provide the correct path."
else
    if [ -d "$USER_R_BIN_PATH" ]; then
        echo "Attempting to add $USER_R_BIN_PATH to PATH for this script session."
        export PATH="$USER_R_BIN_PATH:$PATH"
        echo "Current PATH (first few entries): $(echo "$PATH" | cut -d':' -f1-5)..."
    else
        echo "ERROR: Specified R 'bin' directory path '$USER_R_BIN_PATH' not found. Please ensure the path is correct."
        echo "Script will rely on 'Rscript' being in system PATH."
    fi
fi
# End R Path Addition


# IMPORTANT OS-SPECIFIC PRE-REQUISITES for Development Tools
echo ""
echo "========================================================================"
echo "IMPORTANT: Additional Manual Steps BEFORE Running This Script (OS Specific)"
echo "========================================================================"

if [ "$OS_TYPE" == "windows" ]; then
    echo "1. Windows Locale Setting (To prevent 'LC_CTYPE' warnings in R):"
    echo "   - Go to Control Panel > Clock and Region > Region > Administrative tab."
    echo "   - Under 'Language for non-Unicode programs', click 'Change system locale...'"
    echo "   - Ensure 'Beta: Use Unicode UTF-8 for worldwide language support' is UNCHECKED."
    echo "   - Set the system locale to 'English (United States)' for best compatibility."
    echo ""
    echo "2. Rtools Installation (For compiling R packages from source):"
    echo "   - Ensure Rtools (e.g., Rtools42 for R 4.2.x, Rtools43 for R 4.3.x) is correctly installed."
    echo "   - Crucially, ensure Rtools's 'mingw64/bin' and 'usr/bin' directories (e.g., F:/rtools42/mingw64/bin, F:/rtools42/usr/bin) "
    echo "     are added to your SYSTEM PATH and accessible by R."
    echo "   - You can check if R finds them by running 'Sys.which('make')' in R."

elif [ "$OS_TYPE" == "mac" ]; then
    echo "1. Xcode Command Line Tools (For compiling R packages from source):"
    echo "   - Open Terminal and run: 'xcode-select --install'"
    echo "   - This provides 'make', 'gcc', 'g++', which are necessary for compiling R packages."
    echo "2. System R Locale (Usually handled by R.app):"
    echo "   - R.app typically handles locale settings well. If you encounter warnings, consult R documentation."
    
elif [ "$OS_TYPE" == "linux" ]; then
    echo "1. Build Essentials (For compiling R packages from source):"
    echo "   - Ensure you have build tools installed. For Debian/Ubuntu: 'sudo apt-get install build-essential'"
    echo "   - For Fedora/RHEL: 'sudo yum groupinstall \"Development Tools\"'"
    echo "2. System R Locale (Usually configured correctly):"
    echo "   - If you encounter locale warnings, consult your distribution's documentation on R locale setup."
fi
echo ""
echo "Failure to complete these manual steps can cause R package installations to fail."
echo "========================================================================"
echo ""

echo "Step 1: Checking for system Rscript..."
if ! command -v Rscript &> /dev/null
then
    echo "Rscript command not found in PATH even after attempting to add user-specified path."
    echo "Please perform the following checks:"
    echo "1. Verify R 4.2.3 is installed on your system."
    echo "2. Ensure R's 'bin' directory (e.g., C:\\Program Files\\R\\R-x.x.x\\bin or /Library/Frameworks/R.framework/Versions/R-x.x.x/Resources/bin) "
    echo "   is added to your system's PATH."
    echo "3. Restart your terminal after modifying system PATH."
    echo "4. Double-check the 'R 4.2.3 bin directory path' input you provided earlier, and ensure it accurately points to your R's bin directory "
    echo "   (using Git Bash path format for Windows, or standard Unix path for macOS/Linux)."
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

    # Define package versions to fix for stable installations.
    package_versions <- list(
        "spatstat.utils" = "3.1.4",
        "BiocManager" = "1.30.22",
        "GenomeInfoDbData" = "1.2.9",
        "GenomeInfoDb" = "1.34.9",
        "AnnotationDbi" = "1.60.2",
        "org.Hs.eg.db" = "3.16.0",
        "org.Mm.eg.db" = "3.16.0",
        "fastmap" = "1.2.0",
        "rlang" = "1.1.6",
        "cli" = "3.6.5",
        "htmltools" = "0.5.8.1",
        "digest" = "0.6.35",
        "Rcpp" = "1.0.12",
        "reticulate" = "1.36.1",
        "igraph" = "2.0.3",
        "openssl" = "2.1.2",
        "curl" = "5.2.1",
        "xml2" = "1.3.6",
        "V8" = "4.4.2",
        "sf" = "1.0.16",
        "s2" = "1.1.6",
        "units" = "0.8.5",
        "future" = "1.33.2",
        "promises" = "1.3.0",
        "later" = "1.3.2",
        "httpuv" = "1.6.15",
        "usethis" = "2.2.3",
        "pkgload" = "1.3.4",
        "pkgbuild" = "1.4.4",
        "sessioninfo" = "1.2.2",
        "desc" = "1.4.3",
        "purrr" = "1.0.2",
        "jsonlite" = "1.8.8",
        "httr" = "1.4.7",
        "remotes" = "2.5.0",
        "dplyr" = "1.1.4",
        "ggplot2" = "3.5.1",
        "Matrix" = "1.6.5",
        "tibble" = "3.2.1",
        "tidyr" = "1.3.1",
        "viridis" = "0.6.5",
        "reshape2" = "1.4.4",
        "pheatmap" = "1.0.12",
        "cowplot" = "1.1.3",
        "ggpubr" = "0.6.0",
        "patchwork" = "1.2.0",
        "stringr" = "1.5.1",
        "car" = "3.1.2",
        "ggcorrplot" = "0.1.4.1",
        "homologene" = "1.4.68.19.3.27",
        "TOAST" = "1.12.0",
        "spatstat.geom" = "3.2.9",
        "spatstat.random" = "3.2.3",
        "spatstat.explore" = "3.2.7",
        "spatstat.model" = "3.2.11",
        "devtools" = "2.4.5",
        "R.utils" = "2.12.3",
        "ggunchull" = "1.0.1",
        "Seurat" = "5.0.3",
        "matrixStats" = "1.0.0",
        "fitdistrplus" = "1.1.11",
        "ggridges" = "0.5.6",
        "ica" = "1.0.3",
        "irlba" = "2.3.5.1",
        "lmtest" = "0.9.40",
        "pbapply" = "1.7.2",
        "plotly" = "4.10.4",
        "RANN" = "2.6.1",
        "RcppAnnoy" = "0.0.22",
        "ROCR" = "1.0.11",
        "Rtsne" = "0.17",
        "scattermore" = "1.2",
        "uwot" = "0.2.2",
        "RcppProgress" = "0.4.2",
        "miniUI" = "0.1.1.1",
        "htmlwidgets" = "1.6.4",
        "lazyeval" = "0.2.2",
        "gplots" = "3.1.3.1",
        "ape" = "5.8",
        "ggplotify" = "0.1.2",
        "assertthat" = "0.2.1",
        "leiden" = "0.4.3.1",
        "sctransform" = "0.4.1",
        "SeuratObject" = "5.1.0",
        "harmony" = "1.2.0",
        "celldex" = "1.8.0",
        "decoupleR" = "2.4.0",
        "ensembldb" = "2.22.0",
        "msigdbr" = "24.1.0",
        "scater" = "1.26.1",
        "scDblFinder" = "1.12.0",
        "SCpubr" = "2.0.2",
        "SingleCellExperiment" = "1.20.1",
        "SingleR" = "2.0.0",
        "SpatialExperiment" = "1.8.1",
        "scran" = "1.26.2",
        "UCell" = "2.2.0",
        "OmnipathR" = "3.6.0",
        "scRNAtoolVis" = "0.1.0",
        "SeuratExtend" = "1.2.3",
        "SeuratDisk" = "0.0.0.9021",
        "MuSiC" = "1.0.0",
        "CARD" = "1.1",
        "liana" = "0.1.14",
        "basilisk" = NULL,  # Adjusted to match Bioconductor 3.16 default
        "basilisk.utils" = "1.21.2",
        "ceLLama" = "0.1.0",
        "thinkr" = "0.1.2"
    )

    # Determine a writable library path for package installation
    user_lib_path <- Sys.getenv("R_LIBS_USER", unset = NA)
    if (is.na(user_lib_path) || !file.exists(user_lib_path) || !file.access(user_lib_path, 2) == 0) {
        temp_lib_path <- .libPaths()[1]
        if (!file.exists(temp_lib_path) || !file.access(temp_lib_path, 2) == 0) {
            if (.Platform\$OS.type == "windows") {
                temp_lib_path <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", paste0(R.version\$major, ".", R.version\$minor))
            } else {
                temp_lib_path <- file.path(Sys.getenv("HOME"), "R", paste0(R.version\$major, ".", R.version\$minor))
            }
            cat("Default R library path not writable or missing. Attempting to use/create user-specific library: ", temp_lib_path, "\\n")
        }
        user_lib_path <- temp_lib_path
    }

    cat("Attempting to install packages into library:", user_lib_path, "\\n")
    if (!dir.exists(user_lib_path)) {
      cat("Target library path does not exist:", user_lib_path, "\\n")
      tryCatch(dir.create(user_lib_path, recursive = TRUE, showWarnings = TRUE),
               error = function(e) message(paste0("Failed to create library path: ", user_lib_path, ". Error: ", e$message)),
               warning = function(w) message(paste0("Warning creating library path: ", user_lib_path, ". Warning: ", w$message)))
    }
    .libPaths(c(user_lib_path, .libPaths()))

    cat("R Library Paths (.libPaths()):\\n")
    print(.libPaths())
    cat("-----------------------------\\n")

    # Helper function to install a package only if it's not already installed
    install_if_missing <- function(pkg_name, lib_path, pkg_type = NULL, version = NULL) {
        if (pkg_name %in% names(package_versions)) {
            fixed_version <- package_versions[[pkg_name]]
            if (!is.null(fixed_version)) {
                version <- fixed_version
                cat(paste0("  (Fixed version for '", pkg_name, "': ", version, ")\\n"))
            }
        }

        current_version <- tryCatch(as.character(packageVersion(pkg_name, lib.loc = lib_path)), error = function(e) NA)
        
        if (!is.na(current_version) && (is.null(version) || current_version == version)) {
            cat("Package '", pkg_name, "' (version ", current_version, ") is already installed and matches target. Skipping.\\n")
            return(TRUE)
        }

        cat("Package '", pkg_name, "' not found or version mismatch. Attempting installation...\\n")
        if (!is.na(current_version)) {
            cat(paste0("  (Currently installed version: ", current_version, ", Target version: ", ifelse(is.null(version), "latest", version), ")\\n"))
        }
        
        current_pkg_type <- pkg_type
        if (is.null(current_pkg_type)) {
            if (.Platform\$OS.type == "windows") {
                current_pkg_type <- "binary"
            } else {
                current_pkg_type <- "source"
            }
        }

        install_success <- FALSE
        if (!is.null(version) && pkg_name != "remotes" && requireNamespace("remotes", quietly = TRUE, lib.loc = lib_path)) {
            cat(paste0("  Attempting to install specific version '", version, "' of '", pkg_name, "' using remotes::install_version.\\n"))
            tryCatch({
                if (!isNamespaceLoaded("remotes")) {
                  library(remotes, lib.loc = lib_path)
                }
                remotes::install_version(pkg_name, version = version, Ncpus = max(1, parallel::detectCores() - 1), lib = lib_path, upgrade = "never", force = TRUE)
                if (requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path) && as.character(packageVersion(pkg_name, lib.loc = lib_path)) == version) {
                    cat("Successfully installed '", pkg_name, "' version ", version, ".\\n")
                    install_success <- TRUE
                } else {
                    cat("WARNING: Installed '", pkg_name, "' but version mismatch or failed to load after install_version.\\n")
                }
            }, error = function(e) {
                message(paste0("ERROR installing '", pkg_name, "' v", version, " with remotes::install_version: ", conditionMessage(e)))
            })
        }

        if (!install_success) {
            cat(paste0("  Attempting generic install.packages for '", pkg_name, "' (type=", current_pkg_type, ").\\n"))
            install_args <- list(pkgs = pkg_name, lib = lib_path, Ncpus = max(1, parallel::detectCores() - 1))
            if (!is.null(current_pkg_type)) {
                install_args\$type <- current_pkg_type
            }
            
            tryCatch({
                do.call(install.packages, install_args)
                if (requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path)) {
                    cat("Successfully installed (generic) and loaded namespace for package:", pkg_name, "\\n")
                    install_success <- TRUE
                } else {
                    cat("WARNING: Package '", pkg_name, "' still not found after generic installation attempt.\\n")
                }
            }, error = function(e) {
                message(paste0("ERROR installing '", pkg_name, "' (generic install): ", conditionMessage(e)))
                if (.Platform\$OS.type == "windows" && current_pkg_type == "binary") {
                    message(paste0("Retrying '", pkg_name, "' from source on Windows.\\n"))
                    tryCatch({
                        install_args\$type <- "source"
                        do.call(install.packages, install_args)
                        if (requireNamespace(pkg_name, quietly = TRUE, lib.loc = lib_path)) {
                            cat("Successfully installed (source fallback) and loaded namespace for package:", pkg_name, "\\n")
                            install_success <- TRUE
                        } else {
                            cat("WARNING: Package '", pkg_name, "' still not found after source fallback installation attempt.\\n")
                        }
                    }, error = function(e_src) {
                        message(paste0("ERROR installing '", pkg_name, "' (source fallback) on Windows: ", conditionMessage(e_src)))
                    })
                }
            })
        }
        return(install_success)
    }

    install_and_check_bioc <- function(pkg_name, lib_path, version = NULL) {
        if (pkg_name %in% names(package_versions)) {
            fixed_pkg_version_from_list <- package_versions[[pkg_name]]
            if (!is.null(fixed_pkg_version_from_list)) {
                cat(paste0("  (Fixed package version for '", pkg_name, "': ", fixed_pkg_version_from_list, ") will be checked for skipping.\\n"))
            }
        }

        current_version <- tryCatch(as.character(packageVersion(pkg_name, lib.loc = lib_path)), error = function(e) NA)
        if (!is.na(current_version) && (is.null(fixed_pkg_version_from_list) || current_version == fixed_pkg_version_from_list)) {
            cat("Bioconductor package '", pkg_name, "' (version ", current_version, ") is already installed and matches target. Skipping.\\n")
            return(TRUE)
        }
        cat("Attempting to install Bioconductor package:", pkg_name, "into", lib_path, ".\\n")
        if (!is.na(current_version)) {
            cat(paste0("  (Currently installed version: ", current_version, ", Target version: ", ifelse(is.null(fixed_pkg_version_from_list), "latest for Bioc release", fixed_pkg_version_from_list), ")\\n"))
        }

        installed_ok <- FALSE
        
        tryCatch({ 
            BiocManager::install(pkg_name, lib = lib_path, Ncpus = max(1, parallel::detectCores() - 1), ask = FALSE, update = FALSE, force = TRUE, quiet = FALSE, verbose = TRUE)
            installed_version_after_attempt <- tryCatch(as.character(packageVersion(pkg_name, lib.loc = lib_path)), error = function(e) NA)
            if (!is.na(installed_version_after_attempt) && (is.null(fixed_pkg_version_from_list) || installed_version_after_attempt == fixed_pkg_version_from_list)) {
                cat("Successfully installed Bioconductor package:", pkg_name, " (version ", installed_version_after_attempt, ").\\n")
                installed_ok <- TRUE
            } else {
                cat("Failed to install Bioconductor package:", pkg_name, "or version mismatch after install attempt. Installed: ", ifelse(is.na(installed_version_after_attempt), "None", installed_version_after_attempt), ". Targeted fixed: ", ifelse(is.null(fixed_pkg_version_from_list), "None", fixed_pkg_version_from_list), ".\\n")
            }
        }, error = function(e) {
            cat("ERROR during default BiocManager::install for package:", pkg_name, ". Error: ", conditionMessage(e), ".\\n")
        })

        return(installed_ok)
    }

    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    cat("CRAN mirror set to https://cloud.r-project.org/\\n")

    # Ensure Matrix is installed with the correct version first
    if ("Matrix" %in% names(package_versions)) {
        fixed_version <- package_versions[["Matrix"]]
        cat("Ensuring 'Matrix' v", fixed_version, " is installed...\\n")
        install_if_missing("Matrix", lib_path = user_lib_path, version = fixed_version)
    }

    # Special handling for matrixStats
    pkg_name_matrixStats <- "matrixStats"
    target_version_matrixStats <- package_versions[[pkg_name_matrixStats]]
    cat(paste0("Attempting to ensure '", pkg_name_matrixStats, "' v", target_version_matrixStats, " is cleanly installed...\\n"))

    if (requireNamespace(pkg_name_matrixStats, quietly = TRUE, lib.loc = user_lib_path)) {
        cat(paste0("Removing existing '", pkg_name_matrixStats, "' to ensure clean reinstallation...\\n"))
        tryCatch({
            remove.packages(pkg_name_matrixStats, lib = user_lib_path)
            cat(paste0("Successfully removed existing '", pkg_name_matrixStats, "'.\\n"))
        }, error = function(e) {
            cat(paste0("WARNING: Could not remove existing '", pkg_name_matrixStats, "'. Error: ", conditionMessage(e), ". This might lead to 'cannot rename file' errors or other installation failures if the package is in use.\\n"))
            cat("Please ensure no other R sessions or processes are running that might lock package files (e.g., RStudio, other R consoles, background R processes) and try again.\\n")
        })
    } else {
        cat(paste0("No existing '", pkg_name_matrixStats, "' package found for removal.\\n"))
    }

    tryCatch({
        remotes::install_version(pkg_name_matrixStats, version = target_version_matrixStats, 
                                 Ncpus = max(1, parallel::detectCores() - 1), lib = user_lib_path, 
                                 upgrade = "never", force = TRUE)
        if (requireNamespace(pkg_name_matrixStats, quietly = TRUE, lib.loc = user_lib_path) && 
            as.character(packageVersion(pkg_name_matrixStats, lib.loc = user_lib_path)) == target_version_matrixStats) {
            cat(paste0("Successfully installed '", pkg_name_matrixStats, "' v", target_version_matrixStats, ".\\n"))
        } else {
            cat(paste0("ERROR: Failed to install '", pkg_name_matrixStats, "' v", target_version_matrixStats, " or correct version not found after installation attempt.\\n"))
            message("This is a critical dependency. Please check the R log for specific errors. Often, ensuring no other R sessions are running and retrying solves this.\\n")
        }
    }, error = function(e) {
        cat(paste0("CRITICAL ERROR during installation of '", pkg_name_matrixStats, "' v", target_version_matrixStats, ": ", conditionMessage(e), "\\n"))
        message("This package is a common source of installation issues. Please ensure Rtools (Windows) or build essentials (macOS/Linux) are correctly configured and accessible by R, and that no R sessions are locking files.\\n")
    })

    # Pre-emptive spatstat.utils handling
    cat("--- Pre-emptive spatstat.utils handling ---\\n")
    install_if_missing("spatstat.utils", lib_path = user_lib_path, pkg_type = "source", version = package_versions[["spatstat.utils"]])
    cat("---------------------------------------------\\n")

    install_if_missing("BiocManager", lib_path = user_lib_path, version = package_versions[["BiocManager"]])
    library(BiocManager, lib.loc = user_lib_path)

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
        install_and_check_bioc(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }

    cat("Installing critical CRAN packages (fastmap, rlang, cli, htmltools, digest)...\\n")
    critical_cran_updates <- c("fastmap", "rlang", "cli", "htmltools", "digest")
    for (pkg in critical_cran_updates) {
        install_if_missing(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }
    
    cat("Installing other core CRAN packages and devtools dependencies...\\n")
    cran_core_and_pre_deps_binary <- c(
        "Rcpp", "reticulate", "igraph", "openssl", "curl", "xml2", "V8",
        "sf", "s2", "units", "future", "promises", "later", "httpuv",
        "usethis", "pkgload", "pkgbuild", "sessioninfo", "desc", "purrr", "jsonlite", "httr", "remotes",
        "dplyr", "ggplot2", "tibble", "tidyr", "viridis", "reshape2", "pheatmap", "cowplot", "ggpubr", "patchwork",
        "stringr", "car", "ggcorrplot", "homologene"
    )
    for (pkg in cran_core_and_pre_deps_binary) {
        install_if_missing(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }
    
    install_and_check_bioc("TOAST", lib_path = user_lib_path, version = package_versions[["TOAST"]])

    cat("Attempting to install other spatstat family packages (geom, random, explore, model) from CRAN...\\n")
    spatstat_other_pkgs <- c("spatstat.geom", "spatstat.random", "spatstat.explore", "spatstat.model")
    for (pkg in spatstat_other_pkgs) {
        install_if_missing(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }

    install_if_missing("devtools", lib_path = user_lib_path, version = package_versions[["devtools"]])
    library(devtools, lib.loc = user_lib_path)
    install_if_missing("remotes", lib_path = user_lib_path, version = package_versions[["remotes"]])
    library(remotes, lib.loc = user_lib_path)

    cat("Installing specific CRAN/GitHub dependencies...\\n")
    tryCatch(remotes::install_github("HenrikBengtsson/R.utils", ref="develop", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for HenrikBengtsson/R.utils: ", conditionMessage(e))))
    tryCatch(remotes::install_github("sajuukLyu/ggunchull", upgrade = "never", build_vignettes = FALSE, type = "source", force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for sajuukLyu/ggunchull: ", conditionMessage(e))))

    cat("Installing Seurat from CRAN and its other dependencies explicitly...\\n")
    install_if_missing("Seurat", lib_path = user_lib_path, version = package_versions[["Seurat"]])
    
    seurat_known_deps_cran <- c(
        'fitdistrplus', 'ggridges', 'ica', 'irlba', 'lmtest', 'pbapply', 'plotly', 'RANN',
        'RcppAnnoy', 'ROCR', 'Rtsne', 'scattermore',
        'uwot', 'RcppProgress', 'miniUI', 'htmlwidgets', 'lazyeval', 'gplots', 'ape', 'ggplotify', 'assertthat'
    )
    for (pkg in seurat_known_deps_cran) {
        install_if_missing(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }

    seurat_known_deps_bioc <- c('leiden', 'sctransform', 'SeuratObject')
    cat("Installing Bioconductor dependencies for Seurat...\\n")
    for (pkg in seurat_known_deps_bioc) {
        install_and_check_bioc(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }

    cat("Installing Harmony from CRAN...\\n")
    install_if_missing("harmony", lib_path = user_lib_path, version = package_versions[["harmony"]])

    cat("--- Ensuring celldex is properly installed (attempting remove and reinstall) ---\\n")
    install_and_check_bioc("celldex", lib_path = user_lib_path, version = package_versions[["celldex"]])
    cat("--- Finished celldex reinstallation attempt ---\\n")

    cat("Installing remaining common Bioconductor packages...\\n")
    bioc_packages_remaining <- c(
        "decoupleR", "ensembldb", "msigdbr", "scater",
        "scDblFinder", "SCpubr",
        "SingleCellExperiment", "SingleR", "SpatialExperiment", "scran", "UCell"
    )
    for (pkg in bioc_packages_remaining) {
        install_and_check_bioc(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
    }

    install_and_check_bioc('OmnipathR', lib_path = user_lib_path, version = package_versions[["OmnipathR"]])

    cat("Installing problematic CRAN packages (yulab.utils, ggfun, scatterpie) with fallbacks...\\n")
    problematic_cran_packages <- c("yulab.utils", "ggfun", "scatterpie")
    for (pkg in problematic_cran_packages) {
        if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
            message(paste0("Attempting to install ", pkg, " from CRAN (binary if available) into ", user_lib_path, "...\\n"))
            install_if_missing(pkg, lib_path = user_lib_path, version = package_versions[[pkg]])
            
            if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib_path)) {
                message(paste0("CRAN installation of ", pkg, " failed or package still not found. Attempting to install from GitHub...\\n"))
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

    cat("Installing additional user-specified GitHub packages (ceLLama, thinkr)...\\n")
    tryCatch(devtools::install_github("CelVoxes/ceLLama", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
            error = function(e) message(paste0("GitHub install failed for CelVoxes/ceLLama: ", conditionMessage(e))))

    tryCatch(devtools::install_github("eonurk/thinkR", upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
            error = function(e) message(paste0("GitHub install failed for eonurk/thinkR: ", conditionMessage(e))))

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
    cat("Attempting to ensure 'basilisk' is cleanly installed...\\n")

    if (requireNamespace("basilisk", quietly = TRUE, lib.loc = user_lib_path)) {
      cat("Attempting to remove existing 'basilisk' package...\\n")
      tryCatch(
        remove.packages("basilisk", lib = user_lib_path),
        error = function(e) message(paste0("WARNING: Failed to remove 'basilisk' during pre-installation cleanup. This might indicate file locks. Error: ", conditionMessage(e), "\\n")),
        warning = function(w) message(paste0("WARNING: Problems encountered removing 'basilisk'. Warning: ", conditionMessage(w), "\\n"))
      )
    }
    if (requireNamespace("basilisk.utils", quietly = TRUE, lib.loc = user_lib_path)) {
        cat("Attempting to remove existing 'basilisk.utils' package...\\n")
        tryCatch(
            remove.packages("basilisk.utils", lib = user_lib_path),
            error = function(e) message(paste0("WARNING: Failed to remove 'basilisk.utils' during pre-installation cleanup. Error: ", conditionMessage(e), "\\n"))
        )
    }
      
    if (requireNamespace("basilisk", quietly = TRUE, lib.loc = user_lib_path)) {
      cat("WARNING: 'basilisk' still detected after removal attempt. Please ensure no other R sessions are using it, or try restarting R and rerunning the script.\\n")
    } else {
      cat("'basilisk' successfully removed (if it existed) or not found.\\n")
    }

    install_and_check_bioc("basilisk", lib_path = user_lib_path)
    if (requireNamespace("basilisk", quietly = TRUE, TRUE)) {
      cat("Successfully ensured 'basilisk' is installed and its namespace is available.\\n")
    } else {
      cat("ERROR: 'basilisk' did not install correctly or could not be loaded. 'liana' installation may fail.\\n")
      message("Please check the R log for 'basilisk' installation errors. You might need to manually inspect its dependencies or try installing it separately.\\n")
    }
    
    tryCatch(remotes::install_github('saezlab/liana', upgrade = "never", build_vignettes = FALSE, force = TRUE, lib = user_lib_path),
             error = function(e) message(paste0("GitHub install failed for saezlab/liana. This might be due to 'basilisk' issues. Error: ", conditionMessage(e))))

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

# --- Section to save installed package versions to CSV ---
cat("\n========================================================================\n")
cat("SAVING INSTALLED R PACKAGE VERSIONS TO CSV (For Reference / Pinning)\n")
cat("========================================================================\n")

output_csv_path <- "F:/R_PROJECT/impact_sc/version.csv" 
output_dir <- dirname(output_csv_path)

target_packages <- c(
    "spatstat.utils", "BiocManager",
    "GenomeInfoDbData", "GenomeInfoDb", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
    "fastmap", "rlang", "cli", "htmltools", "digest",
    "Rcpp", "reticulate", "igraph", "openssl", "curl", "xml2", "V8",
    "sf", "s2", "units", "future", "promises", "later", "httpuv",
    "usethis", "pkgload", "pkgbuild", "sessioninfo", "desc", "purrr", "jsonlite", "httr", "remotes",
    "dplyr", "ggplot2", "Matrix", "tibble", "tidyr", "viridis", "reshape2", "pheatmap", "cowplot", "ggpubr", "patchwork",
    "stringr", "car", "ggcorrplot", "homologene",
    "TOAST",
    "spatstat.geom", "spatstat.random", "spatstat.explore", "spatstat.model",
    "devtools",
    "R.utils", "ggunchull",
    "Seurat", "matrixStats",
    "leiden", "sctransform", "SeuratObject",
    "harmony", "celldex",
    "decoupleR", "ensembldb", "msigdbr", "scater",
    "scDblFinder", "SCpubr",
    "SingleCellExperiment", "SingleR", "SpatialExperiment", "scran", "UCell",
    "OmnipathR",
    "yulab.utils", "ggfun", "scatterpie",
    "scRNAtoolVis", "SeuratExtend", "SeuratDisk", "MuSiC", "CARD", "liana",
    "basilisk", "basilisk.utils", "ceLLama", "thinkr"
)

installed_pkg_versions <- data.frame(Package = character(), Version = character(), stringsAsFactors = FALSE)

for (pkg_name in unique(target_packages)) {
    pkg_version <- tryCatch(
        as.character(packageVersion(pkg_name, lib.loc = user_lib_path)),
        error = function(e) NA
    )
    if (!is.na(pkg_version)) {
        installed_pkg_versions <- rbind(installed_pkg_versions, data.frame(Package = pkg_name, Version = pkg_version, stringsAsFactors = FALSE))
    }
}

installed_pkg_versions <- installed_pkg_versions[order(installed_pkg_versions$Package), ]

df_pkgs_to_save <- as.data.frame(installed_pkg_versions, stringsAsFactors = FALSE)
colnames(df_pkgs_to_save) <- c("Package", "Version")

if (!dir.exists(output_dir)) {
    message(paste0("Creating directory for CSV: ", output_dir, "\\n"))
    tryCatch(dir.create(output_dir, recursive = TRUE, showWarnings = TRUE),
             error = function(e) message(paste0("Failed to create directory: ", output_dir, ". Error: ", e$message, "\\n")))
}

tryCatch({
    write.csv(df_pkgs_to_save, file = output_csv_path, row.names = FALSE)
    cat(paste0("Successfully saved package versions to: ", output_csv_path, "\\n"))
}, error = function(e) {
    cat(paste0("ERROR: Failed to save package versions to CSV. Error: ", conditionMessage(e), "\\n"))
    cat(paste0("Please check write permissions for the directory: ", output_dir, " and ensure the path is valid.\\n"))
})

cat("------------------------------------------------------------------------\n")
cat("This CSV contains versions of packages targeted by this script.\n")
cat("To fix specific package versions for future installations, copy the desired\n")
cat("version numbers from the CSV and add 'version=\"X.Y.Z\"' to the corresponding\n")
cat("'install_if_missing' or 'install_and_check_bioc' calls in this script.\n")
cat("Example: install_if_missing(\"dplyr\", lib_path = user_lib_path, version = \"1.0.0\")\n")
cat("========================================================================\n")

EOF

echo "Running R package installation script (this may take a very long time)..."
R_INSTALL_LOG="r_package_install_system_r_only.log"
echo "Full R installation log will be saved to: $(pwd)/$R_INSTALL_LOG"

Rscript "$INSTALL_R_SCRIPT" > "$R_INSTALL_LOG" 2>&1
R_EXIT_CODE=$?

echo "R package installation run is complete."
if [ $R_EXIT_CODE -ne 0 ]; then
    echo "The R package installation script reported a non-zero exit code: $R_EXIT_CODE."
    echo "This indicates issues with one or more packages."
    echo "Please check the detailed log in '$R_INSTALL_LOG' for specific information."
else
    echo "The R package installation script completed successfully (exit code 0)."
    echo "Please check the log file '$R_INSTALL_LOG' to verify the status of all packages."
fi
echo "The R script '$INSTALL_R_SCRIPT' has been kept for inspection."

echo ""
echo "--- R Dependency Installation Attempt Complete (System-Wide Mode) ---"
echo ""

# Python Environment Setup using Conda
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
echo "IMPORTANT REMINDERS (Recap):"
echo "- If R packages failed, address those issues. R is often a prerequisite for parts of the pipeline."
echo "- Ensure no R sessions or processes are running when installing R packages, especially if you encounter 'cannot remove earlier installation' errors."

if [ "$OS_TYPE" == "windows" ]; then
    echo "- Remember to RUN THIS SCRIPT AS ADMINISTRATOR on Windows if installing R packages to system R library to avoid 'Permission Denied' errors."
    echo "- Ensure Rtools (Windows) is correctly set up and accessible by R for source compilation."
elif [ "$OS_TYPE" == "mac" ]; then
    echo "- If R package installation fails with permission errors, try running this script with 'sudo'."
    echo "- Ensure Xcode Command Line Tools are installed (run 'xcode-select --install') for compiling R packages from source on macOS."
elif [ "$OS_TYPE" == "linux" ]; then
    echo "- If R package installation fails with permission errors, try running this script with 'sudo'."
    echo "- Ensure build-essential packages are installed for compiling R packages from source on Linux."
fi

exit 0
