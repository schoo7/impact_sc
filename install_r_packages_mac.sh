#!/bin/bash

echo "==================================================================="
echo "IMPACT-sc R Package Installation for macOS"
echo "==================================================================="
echo "This script installs all required R packages for IMPACT-sc on macOS."
echo "Prerequisites: R and build tools (Xcode CLI) must be installed"
echo ""

# Detect Mac architecture
ARCH=$(uname -m)
if [[ "$ARCH" == "arm64" ]]; then
    echo "🍎 Detected Apple Silicon Mac (M1/M2/M3)"
    MAC_TYPE="apple_silicon"
elif [[ "$ARCH" == "x86_64" ]]; then
    echo "💻 Detected Intel Mac"
    MAC_TYPE="intel"
else
    echo "❓ Unknown Mac architecture: $ARCH"
    MAC_TYPE="unknown"
fi

# Common Mac R installation paths
MAC_R_PATHS=(
    "/opt/homebrew/bin/Rscript"           # Homebrew on Apple Silicon
    "/usr/local/bin/Rscript"              # Homebrew on Intel
    "/Library/Frameworks/R.framework/Resources/bin/Rscript"  # CRAN installer
    "/Applications/R.app/Contents/MacOS/R"  # R.app
)

# Try to find R automatically
echo "Checking for R installation..."
R_FOUND=""
for r_path in "${MAC_R_PATHS[@]}"; do
    if [[ -f "$r_path" ]]; then
        R_FOUND="$r_path"
        echo "✅ Found R at: $R_FOUND"
        break
    fi
done

if [[ -z "$R_FOUND" ]]; then
    echo "❌ R not found in common locations."
    echo "Please install R first:"
    echo "  Option 1: brew install r"
    echo "  Option 2: Download from https://cran.r-project.org/bin/macosx/"
    exit 1
fi

# Check for Xcode Command Line Tools
if ! xcode-select -p &> /dev/null; then
    echo "❌ Xcode Command Line Tools not found."
    echo "Please install: xcode-select --install"
    exit 1
else
    echo "✅ Xcode Command Line Tools detected"
fi

# Set environment variables for Apple Silicon compilation
if [[ "$MAC_TYPE" == "apple_silicon" ]]; then
    echo "Setting up Apple Silicon compilation environment..."
    export LDFLAGS="-L/opt/homebrew/lib"
    export CPPFLAGS="-I/opt/homebrew/include"
    export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig"
fi

LOG_FILE="r_package_install_mac.log"
echo "Installation log will be saved to: $LOG_FILE"
echo ""

# Create R installation script with Mac-specific optimizations
INSTALL_R_SCRIPT="temp_install_r_packages_mac.R"
cat << 'EOF' > "$INSTALL_R_SCRIPT"
# Mac-optimized R package installation script
main_result <- tryCatch({

    cat("--- R Environment Details (macOS) ---\n")
    cat("R.version.string:", R.version.string, "\n")
    cat("Platform:", R.version$platform, "\n")
    cat("Architecture:", Sys.info()[["machine"]], "\n")
    cat("R_HOME:", R.home(), "\n")
    cat("Library paths:", paste(.libPaths(), collapse = "; "), "\n")
    
    # Check build tools
    cat("Build tools availability:\n")
    cat("make:", Sys.which('make'), "\n")
    cat("gcc:", Sys.which('gcc'), "\n")
    cat("g++:", Sys.which('g++'), "\n")
    cat("gfortran:", Sys.which('gfortran'), "\n")
    cat("cmake:", Sys.which('cmake'), "\n")
    
    # Set up library path
    user_lib_path <- .libPaths()[1]
    cat("Installing packages into:", user_lib_path, "\n")
    
    if (!dir.exists(user_lib_path)) {
        dir.create(user_lib_path, recursive = TRUE, showWarnings = TRUE)
    }
    
    # Configure CRAN mirror
    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    options(timeout = 600)
    
    # Mac-specific compilation flags for Apple Silicon
    if (Sys.info()[["machine"]] == "arm64") {
        cat("Configuring for Apple Silicon...\n")
        Sys.setenv(LDFLAGS = "-L/opt/homebrew/lib")
        Sys.setenv(CPPFLAGS = "-I/opt/homebrew/include")
        Sys.setenv(PKG_CONFIG_PATH = "/opt/homebrew/lib/pkgconfig")
    }
    
    # Use parallel installation
    ncpus <- max(1, parallel::detectCores() - 1)
    cat("Using", ncpus, "CPU cores for installation\n")
    
    # Function to install Bioconductor packages with Mac optimizations
    install_bioc_with_retry <- function(pkg_name, max_attempts = 3) {
        cat("Installing Bioconductor package:", pkg_name, "\n")
        
        for (attempt in 1:max_attempts) {
            success <- FALSE
            
            # Try binary first (faster)
            tryCatch({
                BiocManager::install(pkg_name, ask = FALSE, update = FALSE, 
                                   Ncpus = ncpus, force = TRUE)
                if (requireNamespace(pkg_name, quietly = TRUE)) {
                    cat("✅ Successfully installed:", pkg_name, "\n")
                    success <- TRUE
                }
            }, error = function(e) {
                cat("❌ Attempt", attempt, "failed for", pkg_name, ":", e$message, "\n")
            })
            
            if (success) break
            
            # If binary failed and last attempt, try source
            if (attempt == max_attempts) {
                cat("Trying source compilation for", pkg_name, "...\n")
                tryCatch({
                    BiocManager::install(pkg_name, type = "source", ask = FALSE, 
                                       update = FALSE, Ncpus = ncpus, force = TRUE)
                    if (requireNamespace(pkg_name, quietly = TRUE)) {
                        cat("✅ Successfully installed from source:", pkg_name, "\n")
                        success <- TRUE
                    }
                }, error = function(e) {
                    cat("❌ Source compilation also failed for", pkg_name, "\n")
                })
            }
        }
        return(success)
    }
    
    # Function to install CRAN packages with retry
    install_cran_with_retry <- function(pkg_name, max_attempts = 3) {
        for (attempt in 1:max_attempts) {
            success <- FALSE
            tryCatch({
                install.packages(pkg_name, Ncpus = ncpus)
                if (requireNamespace(pkg_name, quietly = TRUE)) {
                    cat("✅ Successfully installed:", pkg_name, "\n")
                    success <- TRUE
                }
            }, error = function(e) {
                cat("❌ Attempt", attempt, "failed for", pkg_name, ":", e$message, "\n")
            })
            if (success) break
        }
        return(success)
    }
    
    # Install BiocManager first
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        cat("Installing BiocManager...\n")
        install_cran_with_retry("BiocManager")
    }
    library(BiocManager)
    
    # Set appropriate Bioconductor version
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
    
    # Core CRAN packages installation
    cat("\n=== Installing Core CRAN Packages ===\n")
    core_cran_packages <- c(
        "devtools", "remotes", "BiocManager",
        "dplyr", "ggplot2", "Matrix", "tibble", "tidyr", 
        "viridis", "reshape2", "pheatmap", "cowplot", "ggpubr", "patchwork",
        "stringr", "reticulate", "future", "future.apply",
        "Rcpp", "RcppAnnoy", "igraph", "jsonlite",
        "openxlsx", "readxl", "data.table"
    )
    
    for (pkg in core_cran_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install_cran_with_retry(pkg)
        } else {
            cat("✅ Already installed:", pkg, "\n")
        }
    }
    
    # Install Seurat
    cat("\n=== Installing Seurat ===\n")
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        install_cran_with_retry("Seurat")
    }
    
    # Core Bioconductor packages
    cat("\n=== Installing Bioconductor Packages ===\n")
    core_bioc_packages <- c(
        "GenomeInfoDb", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
        "SingleCellExperiment", "SingleR", "scater", "scran",
        "celldex", "decoupleR", "UCell", "scDblFinder",
        "clusterProfiler", "enrichplot", "GSVA"
    )
    
    for (pkg in core_bioc_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install_bioc_with_retry(pkg)
        } else {
            cat("✅ Already installed:", pkg, "\n")
        }
    }
    
    # Install GitHub packages
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
    
    # Special handling for harmony and monocle3
    cat("\n=== Installing Special Packages ===\n")
    special_packages <- c("harmony", "monocle3")
    
    for (pkg in special_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat("Installing special package:", pkg, "\n")
            install_cran_with_retry(pkg)
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
    cat("Please check your R installation and build tools.\n")
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
"$R_FOUND" "$INSTALL_R_SCRIPT" 2>&1 | tee "$LOG_FILE"
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
    echo "1. Install missing build tools: xcode-select --install"
    echo "2. Update Homebrew and packages: brew update && brew upgrade"
    echo "3. For Apple Silicon: brew install llvm libomp"
    echo "4. Check internet connection"
fi
echo "==================================================================="

# Apple Silicon specific notes
if [[ "$MAC_TYPE" == "apple_silicon" ]]; then
    echo ""
    echo "🍎 Apple Silicon Notes:"
    echo "- Some packages may take longer to compile"
    echo "- If issues persist, try installing via Rosetta 2:"
    echo "  arch -x86_64 R # Run R in Rosetta mode"
fi

# Clean up
if [ -f "$INSTALL_R_SCRIPT" ]; then
    rm "$INSTALL_R_SCRIPT"
fi

exit $R_EXIT_CODE 