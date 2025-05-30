#!/bin/bash

echo "==================================================================="
echo "IMPACT-sc R Package Installation for Windows"
echo "==================================================================="
echo "This script installs all required R packages for IMPACT-sc on Windows."
echo "Prerequisites: R and Rtools must be installed and in PATH"
echo ""

# Check if running as Administrator
if ! net session >/dev/null 2>&1; then
    echo "❌ ERROR: This script requires Administrator privileges"
    echo "Please run Git Bash as Administrator and try again"
    exit 1
fi

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "❌ ERROR: Rscript not found in PATH"
    echo "Please install R and Rtools, then add to PATH"
    echo "R Download: https://cran.r-project.org/bin/windows/base/"
    echo "Rtools Download: https://cran.r-project.org/bin/windows/Rtools/"
    exit 1
fi

echo "✅ Found R: $(Rscript --version 2>&1 | head -1)"

# Check for Rtools
echo "Checking for Rtools..."
if ! command -v make &> /dev/null; then
    echo "⚠️  WARNING: make command not found - Rtools may not be properly installed"
    echo "Some packages may fail to compile"
else
    echo "✅ Found make: $(make --version | head -1)"
fi

LOG_FILE="r_package_install_windows.log"
echo "Installation log will be saved to: $LOG_FILE"
echo ""

# Create R installation script for Windows
INSTALL_R_SCRIPT="temp_install_r_packages_windows.R"
cat << 'EOF' > "$INSTALL_R_SCRIPT"
# Windows-optimized R package installation script
main_result <- tryCatch({
    
    cat("--- R Environment Details (Windows) ---\n")
    cat("R.version.string:", R.version.string, "\n")
    cat("Platform:", R.version$platform, "\n") 
    cat("R_HOME:", R.home(), "\n")
    cat("Library paths:", paste(.libPaths(), collapse = "; "), "\n")
    
    # Check for Rtools
    has_rtools <- FALSE
    tryCatch({
        sys_info <- Sys.which("make")
        if (sys_info != "") {
            cat("Rtools make found at:", sys_info, "\n")
            has_rtools <- TRUE
        }
    }, error = function(e) {
        cat("Rtools check error:", e$message, "\n")
    })
    
    if (!has_rtools) {
        cat("⚠️  WARNING: Rtools not detected. Some packages may fail to compile.\n")
        cat("Install Rtools from: https://cran.r-project.org/bin/windows/Rtools/\n")
    }
    
    # Set up CRAN mirror
    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    
    # Increase timeout for downloads
    options(timeout = 600)
    
    # Use parallel installation
    ncpus <- max(1, parallel::detectCores() - 1)
    cat("Using", ncpus, "CPU cores for installation\n")
    
    # Function to install packages with retry
    install_with_retry <- function(pkg, max_attempts = 3, from_source = FALSE) {
        for (attempt in 1:max_attempts) {
            cat("Installing", pkg, "(attempt", attempt, "of", max_attempts, ")\n")
            success <- FALSE
            
            tryCatch({
                if (from_source) {
                    install.packages(pkg, type = "source", Ncpus = ncpus)
                } else {
                    install.packages(pkg, type = "binary", Ncpus = ncpus)
                }
                
                # Verify installation
                if (requireNamespace(pkg, quietly = TRUE)) {
                    cat("✅ Successfully installed:", pkg, "\n")
                    success <- TRUE
                    break
                }
            }, error = function(e) {
                cat("❌ Attempt", attempt, "failed for", pkg, ":", e$message, "\n")
                if (attempt == max_attempts) {
                    cat("⚠️  FINAL FAILURE for package:", pkg, "\n")
                }
            })
        }
        return(success)
    }
    
    # Install BiocManager first
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        cat("Installing BiocManager...\n")
        install_with_retry("BiocManager")
    }
    
    library(BiocManager)
    
    # Set Bioconductor version based on R version
    r_version <- paste0(R.version$major, ".", substr(R.version$minor, 1, 1))
    if (r_version >= "4.3") {
        bioc_version <- "3.18"
    } else if (r_version >= "4.2") {
        bioc_version <- "3.16"
    } else {
        bioc_version <- "3.14"
    }
    
    cat("Setting Bioconductor version:", bioc_version, "for R", r_version, "\n")
    BiocManager::install(version = bioc_version, ask = FALSE)
    
    # Core CRAN packages (binary installation preferred on Windows)
    cat("\n=== Installing Core CRAN Packages ===\n")
    core_packages <- c(
        "devtools", "remotes", "BiocManager",
        "dplyr", "ggplot2", "Matrix", "tibble", "tidyr",
        "viridis", "reshape2", "pheatmap", "cowplot", "ggpubr", "patchwork",
        "stringr", "reticulate", "future", "future.apply",
        "Rcpp", "RcppAnnoy", "igraph", "jsonlite",
        "openxlsx", "readxl", "data.table"
    )
    
    for (pkg in core_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install_with_retry(pkg)
        } else {
            cat("✅ Already installed:", pkg, "\n")
        }
    }
    
    # Seurat (critical package)
    cat("\n=== Installing Seurat ===\n")
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        install_with_retry("Seurat")
    }
    
    # Bioconductor packages
    cat("\n=== Installing Bioconductor Packages ===\n")
    bioc_packages <- c(
        "GenomeInfoDb", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
        "SingleCellExperiment", "SingleR", "scater", "scran",
        "celldex", "decoupleR", "UCell", "scDblFinder",
        "clusterProfiler", "enrichplot", "GSVA"
    )
    
    for (pkg in bioc_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat("Installing Bioconductor package:", pkg, "\n")
            tryCatch({
                BiocManager::install(pkg, ask = FALSE, update = FALSE)
                if (requireNamespace(pkg, quietly = TRUE)) {
                    cat("✅ Successfully installed:", pkg, "\n")
                } else {
                    cat("❌ Failed to install:", pkg, "\n")
                }
            }, error = function(e) {
                cat("❌ Error installing", pkg, ":", e$message, "\n")
            })
        } else {
            cat("✅ Already installed:", pkg, "\n")
        }
    }
    
    # GitHub packages (may require source compilation)
    cat("\n=== Installing GitHub Packages ===\n")
    library(devtools)
    library(remotes)
    
    github_packages <- list(
        "mojaveazure/seurat-disk" = "SeuratDisk",
        "huayc09/SeuratExtend" = "SeuratExtend", 
        "junjunlab/scRNAtoolVis" = "scRNAtoolVis",
        "YingMa0107/CARD" = "CARD"
    )
    
    for (repo in names(github_packages)) {
        pkg_name <- github_packages[[repo]]
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
            cat("Installing from GitHub:", repo, "\n")
            tryCatch({
                remotes::install_github(repo, upgrade = "never", force = TRUE)
                if (requireNamespace(pkg_name, quietly = TRUE)) {
                    cat("✅ Successfully installed:", pkg_name, "\n")
                } else {
                    cat("❌ Failed to install:", pkg_name, "\n")
                }
            }, error = function(e) {
                cat("❌ Error installing", pkg_name, "from GitHub:", e$message, "\n")
            })
        } else {
            cat("✅ Already installed:", pkg_name, "\n")
        }
    }
    
    # Special packages that may need source compilation
    cat("\n=== Installing Special Packages ===\n")
    special_packages <- c("harmony", "monocle3")
    
    for (pkg in special_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat("Installing special package:", pkg, "\n")
            # Try binary first, then source
            if (!install_with_retry(pkg, max_attempts = 2)) {
                cat("Binary install failed, trying source compilation...\n")
                install_with_retry(pkg, max_attempts = 2, from_source = TRUE)
            }
        } else {
            cat("✅ Already installed:", pkg, "\n")
        }
    }
    
    # Final verification
    cat("\n=== Installation Summary ===\n")
    essential_packages <- c("Seurat", "SingleR", "ggplot2", "dplyr", "BiocManager")
    all_good <- TRUE
    
    for (pkg in essential_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            cat("✅", pkg, ":", as.character(packageVersion(pkg)), "\n")
        } else {
            cat("❌", pkg, ": NOT INSTALLED\n")
            all_good <- FALSE
        }
    }
    
    if (all_good) {
        cat("\n🎉 R package installation completed successfully!\n")
    } else {
        cat("\n⚠️  Some essential packages failed to install. Check the log for details.\n")
    }
    
    return(TRUE)
    
}, error = function(e) {
    cat("\n❌ CRITICAL ERROR during R package installation:\n")
    cat("Error message:", e$message, "\n")
    cat("Please check your R and Rtools installation.\n")
    return(FALSE)
})

# Exit with appropriate code
if (!main_result) {
    quit(status = 1, save = "no")
}
EOF

echo "Starting R package installation..."
echo "This may take 30-60 minutes depending on your system and internet connection."
echo ""

# Run R installation with logging
Rscript "$INSTALL_R_SCRIPT" 2>&1 | tee "$LOG_FILE"
R_EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "==================================================================="
if [ $R_EXIT_CODE -eq 0 ]; then
    echo "✅ R package installation completed successfully!"
    echo "Check the summary in the log file: $LOG_FILE"
else
    echo "❌ R package installation encountered errors"
    echo "Check the log file for details: $LOG_FILE"
    echo ""
    echo "Common solutions:"
    echo "1. Ensure Rtools is properly installed and in PATH"
    echo "2. Run Git Bash as Administrator"
    echo "3. Check internet connection"
    echo "4. Update R to the latest version"
fi
echo "==================================================================="

# Clean up
if [ -f "$INSTALL_R_SCRIPT" ]; then
    rm "$INSTALL_R_SCRIPT"
fi

exit $R_EXIT_CODE 