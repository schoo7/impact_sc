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
    packages <- c('Seurat', 'SingleR', 'ggplot2', 'BiocParallel', 'celldex')
    BiocManager::install(packages, ask = FALSE, update = TRUE)
    "
    
    # Install CARD from GitHub
    Rscript -e "
    if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools')
    devtools::install_github('YingMa0107/CARD')
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
        log "macOS detected"
        if [[ "$ARCH" == "arm64" ]]; then
            log "Apple Silicon detected - ensuring Rosetta compatibility"
            export LDFLAGS="-L/opt/homebrew/lib"
            export CPPFLAGS="-I/opt/homebrew/include"
        fi
        ;;
    Linux)
        # Linux specific setup
        log "Linux detected"
        
        # Detect Linux distribution
        if [[ -f /etc/os-release ]]; then
            . /etc/os-release
            DISTRO=$ID
            log "Distribution: $DISTRO"
        else
            DISTRO="unknown"
            log "WARNING: Could not detect Linux distribution"
        fi
        
        # Install system dependencies based on distribution
        case "$DISTRO" in
            ubuntu|debian)
                log "Installing system dependencies for Debian/Ubuntu"
                if command -v apt-get &> /dev/null; then
                    sudo apt-get update
                    sudo apt-get install -y build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
                        libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
                        libtiff5-dev libjpeg-dev cmake pkg-config gfortran
                fi
                ;;
            centos|rhel|fedora)
                log "Installing system dependencies for RHEL/CentOS/Fedora"
                if command -v yum &> /dev/null; then
                    sudo yum groupinstall -y "Development Tools"
                    sudo yum install -y libcurl-devel openssl-devel libxml2-devel \
                        fontconfig-devel harfbuzz-devel fribidi-devel freetype-devel \
                        libpng-devel libtiff-devel libjpeg-turbo-devel cmake pkgconfig gcc-gfortran
                elif command -v dnf &> /dev/null; then
                    sudo dnf groupinstall -y "Development Tools"
                    sudo dnf install -y libcurl-devel openssl-devel libxml2-devel \
                        fontconfig-devel harfbuzz-devel fribidi-devel freetype-devel \
                        libpng-devel libtiff-devel libjpeg-turbo-devel cmake pkgconfig gcc-gfortran
                fi
                ;;
            arch)
                log "Installing system dependencies for Arch Linux"
                if command -v pacman &> /dev/null; then
                    sudo pacman -S --needed base-devel curl openssl libxml2 fontconfig \
                        harfbuzz fribidi freetype2 libpng libtiff libjpeg-turbo cmake pkgconfig gcc-fortran
                fi
                ;;
            *)
                log "WARNING: Unknown Linux distribution. You may need to install build tools manually."
                ;;
        esac
        ;;
    MINGW*|CYGWIN*|MSYS*)
        # Windows specific setup
        log "Windows (Git Bash/MSYS) detected"

        # Detect R installation
        if ! command -v R &> /dev/null; then
            log "R not found in PATH."
            log "Please provide the full path to the R executable (e.g., C:/Program Files/R/R-4.2.3/bin/R.exe)."
            log "To find the R path:"
            log "1. Open File Explorer and navigate to 'C:\\Program Files\\R'."
            log "2. Look for a folder named 'R-<version>' (e.g., 'R-4.2.3')."
            log "3. Inside, go to 'bin' and note the path to 'R.exe'."
            log "4. Enter the path below (use forward slashes or double backslashes)."
            
            read -p "Enter the path to R.exe: " R_PATH
            if [[ -f "$R_PATH" ]]; then
                R_DIR=$(dirname "$R_PATH")
                export PATH="$R_DIR:$PATH"
                log "Added $R_DIR to PATH."
            else
                fail "Invalid R path. Please ensure the path points to R.exe and try again."
            fi
        fi

        # Verify Rtools
        if ! command -v make &> /dev/null; then
            log "Rtools not found in PATH. Checking common locations..."
            if [[ -d "/c/Rtools" ]]; then
                export PATH="/c/Rtools/bin:/c/Rtools/mingw_64/bin:$PATH"
            else
                log "Rtools not found. Please ensure Rtools is installed and added to PATH."
                log "Download Rtools from: https://cran.r-project.org/bin/windows/Rtools/"
                fail "Rtools is required for compiling R packages on Windows."
            fi
        fi

        # Set R library path to avoid permission issues
        export R_LIBS_USER="$HOME/R/win-library/$(echo "$R_VERSION" | awk -F. '{print $1"."$2}')"
        mkdir -p "$R_LIBS_USER"
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