#!/bin/bash
# Data download script for IMPACT-sc
# Downloads demo data, models, and reference data

set -euo pipefail

# -------------------------
# Configuration
# -------------------------
LOG_FILE="download_data.log"
DATA_DIR="data"
DEMO_DIR="$DATA_DIR/demo"
MODELS_DIR="$DATA_DIR/models"
REFERENCE_DIR="$DATA_DIR/reference"

# URLs and file information
PBMC3K_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
PBMC3K_FILE="pbmc3k_filtered_gene_bc_matrices.tar.gz"

# URL for the default Azimuth reference dataset
AZIMUTH_REF_URL="https://seurat.nygenome.org/azimuth/demo_datasets/bmcite_demo.rds"
AZIMUTH_REF_FILE="bmcite_demo.rds"


# --- Attempt to add user-specified R path to the script's PATH ---
# This is useful if R is not in the system PATH and was not set by install_dependencies.sh
# The path F:\R-4.2.3\bin becomes /f/R-4.2.3/bin in Git Bash on Windows.
# !!!! CRITICAL: CHANGE THIS PATH TO YOUR ACTUAL R INSTALLATION BIN DIRECTORY !!!!
# Example: If R is installed in C:\Program Files\R\R-4.3.0, then use "/c/Program Files/R/R-4.3.0/bin"
USER_R_BIN_PATH_DOWNLOAD="/f/R-4.2.3/bin/x64" 


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

# Function to try and ensure Rscript is findable
ensure_rscript_path() {
    if ! command -v Rscript &> /dev/null; then
        log "Rscript not found in current PATH. Attempting to use USER_R_BIN_PATH_DOWNLOAD: $USER_R_BIN_PATH_DOWNLOAD"
        if [ -d "$USER_R_BIN_PATH_DOWNLOAD" ]; then
            export PATH="$USER_R_BIN_PATH_DOWNLOAD:$PATH"
            log "Added $USER_R_BIN_PATH_DOWNLOAD to PATH for this script session."
            if ! command -v Rscript &> /dev/null; then
                log "Rscript still not found after adding USER_R_BIN_PATH_DOWNLOAD."
                return 1 # Indicate failure
            else
                log "Rscript found after adding USER_R_BIN_PATH_DOWNLOAD: $(command -v Rscript)"
                return 0 # Indicate success
            fi
        else
            log "USER_R_BIN_PATH_DOWNLOAD ('$USER_R_BIN_PATH_DOWNLOAD') does not exist or is not a directory."
            return 1 # Indicate failure
        fi
    else
        log "Rscript found in PATH: $(command -v Rscript)"
        return 0 # Indicate success
    fi
}


create_directories() {
    log "Creating data directories"
    mkdir -p "$DEMO_DIR" "$MODELS_DIR" "$REFERENCE_DIR"
}

download_demo_data() {
    log "Downloading PBMC3k demo data (10x Genomics)"
    
    # Ensure we are in the correct directory for download, then return
    current_dir=$(pwd)
    mkdir -p "$DEMO_DIR" # Ensure demo directory exists
    cd "$DEMO_DIR" || fail "Could not change to directory $DEMO_DIR"
    
    # Download if not exists
    if [[ ! -f "$PBMC3K_FILE" ]]; then
        log "Downloading $PBMC3K_FILE..."
        if command -v wget &> /dev/null; then
            wget --no-check-certificate "$PBMC3K_URL" -O "$PBMC3K_FILE" # Added --no-check-certificate for potential SSL issues
        elif command -v curl &> /dev/null; then
            curl -kL "$PBMC3K_URL" -o "$PBMC3K_FILE" # Added -k for potential SSL issues
        else
            fail "Neither wget nor curl found. Please install one of them."
        fi
    else
        log "$PBMC3K_FILE already exists, skipping download."
    fi
    
    # Extract if not already extracted
    if [[ ! -d "filtered_gene_bc_matrices" ]]; then
        log "Extracting PBMC3k data..."
        tar -xzf "$PBMC3K_FILE" || fail "Failed to extract $PBMC3K_FILE"
        log "Demo data extracted to $DEMO_DIR/filtered_gene_bc_matrices"
    else
        log "Demo data already extracted in $DEMO_DIR/filtered_gene_bc_matrices"
    fi
    
    cd "$current_dir" || fail "Could not return to original directory $current_dir"
}

download_models() {
    log "Downloading Cell2Sentence model"

    # Check if conda command is available
    if ! command -v conda &> /dev/null; then
        fail "Conda command not found. Please ensure conda is in your PATH and install_dependencies.sh was run successfully."
    fi

    # Get Conda base path
    CONDA_BASE_DIR=$(conda info --base 2>/dev/null)
    if [[ $? -ne 0 || -z "$CONDA_BASE_DIR" ]]; then
        fail "Could not determine Conda base directory. Is conda initialized correctly? Try running 'conda init bash' and restarting your shell."
    fi
    log "Conda base directory found: $CONDA_BASE_DIR"

    log "Attempting to initialize Conda environment for script session..."
    # Try eval hook first, as it's often more robust in Git Bash
    OLD_SHELL_OPTS=$(set +o) 
    set +e 
    eval "$(conda shell.bash hook)"
    EVAL_STATUS=$?
    
    # Restore shell options carefully
    set +eu; # Temporarily disable errexit and nounset for safety during restore
    if [[ "$OLD_SHELL_OPTS" == *pipefail* ]]; then set -o pipefail; else set +o pipefail; fi
    if [[ "$OLD_SHELL_OPTS" == *errexit* ]]; then set -o errexit; else set +o errexit; fi
    if [[ "$OLD_SHELL_OPTS" == *nounset* ]]; then set -o nounset; else set +o nounset; fi
    # Re-enable standard strict mode if it was on before
    if [[ "$OLD_SHELL_OPTS" == *errexit* && "$OLD_SHELL_OPTS" == *nounset* && "$OLD_SHELL_OPTS" == *pipefail* ]]; then
        set -euo pipefail
    elif [[ "$OLD_SHELL_OPTS" == *errexit* && "$OLD_SHELL_OPTS" == *nounset* ]]; then
        set -eu
    fi


    if [[ $EVAL_STATUS -ne 0 ]]; then
        log "WARNING: 'eval \"\$(conda shell.bash hook)\"' failed with status $EVAL_STATUS. Attempting to source conda.sh directly."
        CONDA_SH_PATH="$CONDA_BASE_DIR/etc/profile.d/conda.sh"
        if [[ -f "$CONDA_SH_PATH" ]]; then
            # shellcheck source=/dev/null
            source "$CONDA_SH_PATH"
            log "Sourced conda.sh from $CONDA_SH_PATH"
        else
            fail "Failed to initialize Conda using both 'eval \$(conda shell.bash hook)' and sourcing conda.sh. Ensure Conda is correctly installed and initialized for bash (e.g., run 'conda init bash' and restart terminal)."
        fi
    else
        log "Successfully initialized Conda using 'eval \$(conda shell.bash hook)'."
    fi
    
    # Activate the target environment
    log "Attempting to activate impact_sc environment..."
    conda activate impact_sc || fail "Failed to activate impact_sc environment. Ensure it was created by install_dependencies.sh and Conda is initialized correctly for this shell session."
    log "Successfully activated impact_sc conda environment."
    
    # Download model using Python
    mkdir -p "$MODELS_DIR"
    log "Attempting to download model to $MODELS_DIR"

    # Using python from the activated conda environment
    python -c "
import os
from transformers import pipeline, AutoTokenizer, AutoModelForCausalLM
import torch # Keep torch import for consistency, though not directly used for download path

# Make models_dir path absolute for robustness inside Python
models_dir_abs = os.path.abspath('$MODELS_DIR')
os.makedirs(models_dir_abs, exist_ok=True) # Ensure directory exists

model_name = 'vandijklab/C2S-Pythia-410m-cell-type-prediction'

print(f'Downloading {model_name} to cache directory: {models_dir_abs}')
print('This may take several minutes...')

try:
    # Download tokenizer and model, explicitly using cache_dir
    tokenizer = AutoTokenizer.from_pretrained(model_name, cache_dir=models_dir_abs, local_files_only=False)
    model = AutoModelForCausalLM.from_pretrained(model_name, cache_dir=models_dir_abs, local_files_only=False)
    
    print('✅ Cell2Sentence model components downloaded successfully!')
    print(f'Model components cached in: {models_dir_abs}')
    
except Exception as e:
    print(f'❌ Error downloading/loading model: {e}')
    raise # Re-raise the exception to make it clearer if Python script fails
    " || fail "Python script for model download failed."
    
    # Deactivate conda environment
    conda deactivate || log "Warning: Failed to deactivate impact_sc environment."
    log "Deactivated impact_sc conda environment."
}

download_reference_data() {
    log "Downloading reference data (Azimuth PBMC CITE-seq Demo)"
    
    mkdir -p "$REFERENCE_DIR"
    local ref_file_path="$REFERENCE_DIR/$AZIMUTH_REF_FILE"

    if [[ ! -f "$ref_file_path" ]]; then
        log "Downloading $AZIMUTH_REF_FILE..."
        if command -v wget &> /dev/null; then
            wget --no-check-certificate "$AZIMUTH_REF_URL" -O "$ref_file_path"
        elif command -v curl &> /dev/null; then
            curl -kL "$AZIMUTH_REF_URL" -o "$ref_file_path"
        else
            fail "Neither wget nor curl found. Please install one of them to download the reference data."
        fi
        log "✅ Reference data downloaded to: $ref_file_path"
    else
        log "$AZIMUTH_REF_FILE already exists, skipping download."
    fi
}

check_disk_space() {
    log "Checking available disk space"
    local required_gb=5 # GB
    
    if command -v df &> /dev/null; then
        available_kb=$(df -k . | awk 'NR==2 {print $4}') 
        if [[ -n "$available_kb" && "$available_kb" -gt 0 ]]; then
            available_gb=$((available_kb / 1024 / 1024))
            log "Available disk space: ${available_gb}GB"
            if [[ $available_gb -lt $required_gb ]]; then
                log "WARNING: Low disk space. Available: ${available_gb}GB. Recommended: ${required_gb}GB+"
                read -p "Continue anyway? (y/N): " -n 1 -r
                echo # Move to a new line
                if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                    fail "Download cancelled due to insufficient disk space."
                fi
            else
                log "Sufficient disk space available."
            fi
        else
            log "Warning: Could not accurately determine available disk space. Proceeding with caution."
        fi
    else
        log "Warning: 'df' command not found. Unable to check disk space. Proceeding with caution."
    fi
}

verify_downloads() {
    log "Verifying downloaded data"
    local all_ok=true
    
    if [[ -d "$DEMO_DIR/filtered_gene_bc_matrices/hg19" ]]; then 
        log "✅ Demo data: PBMC3k data available ($DEMO_DIR/filtered_gene_bc_matrices/hg19)"
    else
        log "❌ Demo data: Missing or not properly extracted ($DEMO_DIR/filtered_gene_bc_matrices/hg19)"
        all_ok=false
    fi
    
    # More robust check for Hugging Face model cache
    # Look for a 'snapshots' directory and then a subdirectory with a commit hash
    # and then common model files within that hash directory.
    local model_base_name="models--vandijklab--C2S-Pythia-410m-cell-type-prediction"
    local model_path_check="$MODELS_DIR/$model_base_name/snapshots"
    if [[ -d "$model_path_check" ]]; then
        # Find any subdirectory in snapshots (this would be a commit hash)
        local commit_hash_dir
        commit_hash_dir=$(find "$model_path_check" -mindepth 1 -maxdepth 1 -type d -print -quit)
        if [[ -n "$commit_hash_dir" && \
              ( -f "$commit_hash_dir/pytorch_model.bin" || -f "$commit_hash_dir/model.safetensors" ) && \
              -f "$commit_hash_dir/config.json" && \
              -f "$commit_hash_dir/tokenizer.json" ]]; then
            log "✅ Model data: Cell2Sentence model components appear to be available in $MODELS_DIR"
        else
            log "❌ Model data: Cell2Sentence model components (bin/safetensors, config, tokenizer) seem to be missing or incomplete within $commit_hash_dir (under $MODELS_DIR)."
            all_ok=false
        fi
    else
        log "❌ Model data: Expected Hugging Face cache structure for Cell2Sentence model not found in $MODELS_DIR (missing $model_base_name/snapshots)."
        all_ok=false
    fi
    
    local ref_file_path_verify="$REFERENCE_DIR/$AZIMUTH_REF_FILE"
    if [[ -f "$ref_file_path_verify" && $(stat -c%s "$ref_file_path_verify" 2>/dev/null || stat -f%z "$ref_file_path_verify" 2>/dev/null || echo 0) -gt 1000000 ]]; then 
        log "✅ Reference data: $AZIMUTH_REF_FILE available ($REFERENCE_DIR)"
    else
        log "❌ Reference data: $AZIMUTH_REF_FILE missing or too small ($REFERENCE_DIR)"
        all_ok=false
    fi

    if [[ "$all_ok" = true ]]; then
        log "All main data components appear to be successfully downloaded/extracted."
    else
        log "Some data components might be missing or incomplete. Please check logs."
    fi
}

# -------------------------
# Main Execution
# -------------------------

# Clear log file at start
>"$LOG_FILE"

log "Starting IMPACT-sc data download"
log "This will download demo data, models, and reference data."
log "Approximate total download size (can vary): ~3-5GB (some models/references can be large)."
log "Ensure you have a stable internet connection."

check_disk_space
create_directories

download_demo_data
download_models
download_reference_data

verify_downloads

log "Data download script completed!"
log "Data should be structured as follows:"
log "├── $DATA_DIR/"
log "│   ├── demo/             - Demo datasets (e.g., PBMC3k)"
log "│   ├── models/           - Pre-trained models (e.g., Cell2Sentence)"
log "│   └── reference/        - Reference datasets (e.g., $AZIMUTH_REF_FILE)"
log ""
log "You should now be able to run the IMPACT-sc pipeline with the downloaded data."
