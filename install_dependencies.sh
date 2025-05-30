#!/bin/bash

# Unified dependency installer for IMPACT-sc
# Handles cross-platform compatibility and version management

set -euo pipefail

# -------------------------
# Configuration
# -------------------------
LOG_FILE="install_dependencies.log"
R_VERSION=$(R --version | head -1 | awk '{print $3}')
OS_NAME=$(uname -s)
ARCH=$(uname -m)

# Bioconductor versions mapped to R versions
BIOC_VERSIONS="4.4:3.19 4.3:3.18 4.2:3.16 4.1:3.14 4.0:3.12"

# -------------------------
# Functions
# -------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

fail() {
    log "ERROR: $*"
    exit 1
}

detect_bioc_version() {
    local r_major_minor=$(echo "$R_VERSION" | awk -F. '{print $1"."$2}')
    
    # Find matching Bioconductor version
    for pair in $BIOC_VERSIONS; do
        r_ver=${pair%%:*}
        bioc_ver=${pair#*:}
        if [[ "$r_major_minor" == "$r_ver" ]]; then
            echo "$bioc_ver"
            return 0
        fi
    done
    
    fail "Unsupported R version: $R_VERSION. Supported versions: $(echo "$BIOC_VERSIONS" | tr ' ' ',')"
}

install_r_packages() {
    local bioc_version=$(detect_bioc_version)
    log "Installing R packages for Bioconductor $bioc_version (R $R_VERSION)"
    
    # Install BiocManager for the detected version
    Rscript -e "
    if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
    BiocManager::install(version = '$bioc_version', ask = FALSE)
    "
    
    # Core packages
    Rscript -e "
    packages <- c('Seurat', 'SingleR', 'ggplot2', 'BiocParallel')
    BiocManager::install(packages, ask = FALSE, update = TRUE)
    "
}

install_python_deps() {
    log "Setting up Python environment"
    
    # Check for conda/mamba
    if ! command -v mamba &> /dev/null; then
        if ! command -v conda &> /dev/null; then
            fail "Neither conda nor mamba found. Please install Miniconda first."
        fi
        CONDA_CMD="conda"
    else
        CONDA_CMD="mamba"
    fi
    
    # Remove existing environment if present
    if $CONDA_CMD env list | grep -q "^impact_sc"; then
        log "Removing existing impact_sc environment"
        $CONDA_CMD env remove -n impact_sc
    fi
    
    # Create fresh environment
    log "Creating new impact_sc environment"
    $CONDA_CMD env create -f environment.yml
    
    # Verify installation
    log "Verifying Python packages"
    $CONDA_CMD run -n impact_sc python -c "
import pandas, scanpy, torch
print(f'Pandas {pandas.__version__}, Scanpy {scanpy.__version__}, PyTorch {torch.__version__}')
print('✅ All Python packages installed successfully!')
    "
}

# -------------------------
# Main Execution
# -------------------------

log "Starting IMPACT-sc dependency installation"
log "Detected system: $OS_NAME $ARCH"
log "Detected R version: $R_VERSION"

case "$OS_NAME" in
    Darwin)
        # macOS specific setup
        if [[ "$ARCH" == "arm64" ]]; then
            log "Apple Silicon detected - ensuring Rosetta compatibility"
            export LDFLAGS="-L/opt/homebrew/lib"
            export CPPFLAGS="-I/opt/homebrew/include"
        fi
        ;;
    Linux)
        # Linux specific setup
        ;;
    MINGW*|CYGWIN*|MSYS*)
        # Windows specific setup
        ;;
    *)
        fail "Unsupported operating system: $OS_NAME"
        ;;
esac

# Install dependencies
install_r_packages
install_python_deps

log "Installation completed successfully!"
log "Please run 'conda activate impact_sc' before using IMPACT-sc" 