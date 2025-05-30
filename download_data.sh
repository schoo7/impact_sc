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

create_directories() {
    log "Creating data directories"
    mkdir -p "$DEMO_DIR" "$MODELS_DIR" "$REFERENCE_DIR"
}

download_demo_data() {
    log "Downloading PBMC3k demo data (10x Genomics)"
    
    cd "$DEMO_DIR"
    
    # Download if not exists
    if [[ ! -f "$PBMC3K_FILE" ]]; then
        log "Downloading $PBMC3K_FILE..."
        if command -v wget &> /dev/null; then
            wget "$PBMC3K_URL" -O "$PBMC3K_FILE"
        elif command -v curl &> /dev/null; then
            curl -L "$PBMC3K_URL" -o "$PBMC3K_FILE"
        else
            fail "Neither wget nor curl found. Please install one of them."
        fi
    else
        log "Demo data already exists, skipping download"
    fi
    
    # Extract if not already extracted
    if [[ ! -d "filtered_gene_bc_matrices" ]]; then
        log "Extracting PBMC3k data..."
        tar -xzf "$PBMC3K_FILE"
        log "Demo data extracted to $DEMO_DIR/filtered_gene_bc_matrices"
    else
        log "Demo data already extracted"
    fi
    
    cd - > /dev/null
}

download_models() {
    log "Downloading Cell2Sentence model"
    
    # Activate conda environment
    if command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate impact_sc || fail "Failed to activate impact_sc environment. Run install_dependencies.sh first."
    else
        fail "Conda not found. Please run install_dependencies.sh first."
    fi
    
    # Download model using Python
    python -c "
import os
from transformers import pipeline, AutoTokenizer, AutoModelForCausalLM
import torch

models_dir = '$MODELS_DIR'
model_name = 'vandijklab/C2S-Pythia-410m-cell-type-prediction'

print(f'Downloading {model_name}...')
print('This may take several minutes...')

try:
    # Download tokenizer and model
    tokenizer = AutoTokenizer.from_pretrained(model_name, cache_dir=models_dir)
    model = AutoModelForCausalLM.from_pretrained(model_name, cache_dir=models_dir)
    
    # Test the pipeline
    pipe = pipeline('text-generation', model=model_name, cache_dir=models_dir)
    print('✅ Cell2Sentence model downloaded successfully!')
    print(f'Model cached in: {models_dir}')
    
except Exception as e:
    print(f'❌ Error downloading model: {e}')
    exit(1)
    "
}

download_reference_data() {
    log "Downloading reference data (HumanPrimaryCellAtlasData)"
    
    # Download reference data using R
    Rscript -e "
    library(celldex)
    
    # Set reference data directory
    ref_dir <- '$REFERENCE_DIR'
    dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
    
    cat('Downloading HumanPrimaryCellAtlasData...\\n')
    cat('This may take several minutes...\\n')
    
    tryCatch({
        # Download reference data
        ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)
        
        # Save reference data
        ref_file <- file.path(ref_dir, 'HumanPrimaryCellAtlasData.rds')
        saveRDS(ref.data, ref_file)
        
        cat('✅ Reference data downloaded successfully!\\n')
        cat('Reference data saved to:', ref_file, '\\n')
        cat('Data dimensions:', dim(ref.data), '\\n')
        
    }, error = function(e) {
        cat('❌ Error downloading reference data:', e\$message, '\\n')
        quit(status = 1)
    })
    "
}

check_disk_space() {
    log "Checking available disk space"
    
    # Check available space (require at least 5GB)
    if command -v df &> /dev/null; then
        available_kb=$(df . | tail -1 | awk '{print $4}')
        available_gb=$((available_kb / 1024 / 1024))
        
        if [[ $available_gb -lt 5 ]]; then
            log "WARNING: Low disk space. Available: ${available_gb}GB. Recommended: 5GB+"
            read -p "Continue anyway? (y/N): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                fail "Download cancelled due to insufficient disk space"
            fi
        else
            log "Sufficient disk space available: ${available_gb}GB"
        fi
    fi
}

verify_downloads() {
    log "Verifying downloaded data"
    
    # Check demo data
    if [[ -d "$DEMO_DIR/filtered_gene_bc_matrices" ]]; then
        log "✅ Demo data: PBMC3k data available"
    else
        log "❌ Demo data: Missing"
    fi
    
    # Check model cache
    if [[ -d "$MODELS_DIR" && $(find "$MODELS_DIR" -name "*C2S*" -o -name "*Pythia*" | wc -l) -gt 0 ]]; then
        log "✅ Model data: Cell2Sentence model available"
    else
        log "❌ Model data: Missing"
    fi
    
    # Check reference data
    if [[ -f "$REFERENCE_DIR/HumanPrimaryCellAtlasData.rds" ]]; then
        log "✅ Reference data: HumanPrimaryCellAtlasData available"
    else
        log "❌ Reference data: Missing"
    fi
}

# -------------------------
# Main Execution
# -------------------------

log "Starting IMPACT-sc data download"
log "This will download demo data, models, and reference data"
log "Total download size: ~3-5GB"

# Check prerequisites
check_disk_space

# Create directory structure
create_directories

# Download all data
download_demo_data
download_models
download_reference_data

# Verify downloads
verify_downloads

log "Data download completed!"
log "Data structure:"
log "├── $DEMO_DIR/ - Demo datasets (PBMC3k)"
log "├── $MODELS_DIR/ - Pre-trained models (Cell2Sentence)"
log "└── $REFERENCE_DIR/ - Reference data (HumanPrimaryCellAtlas)"
log ""
log "You can now run the IMPACT-sc pipeline with the downloaded data." 