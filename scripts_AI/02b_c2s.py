# IMPACT-sc Script: 02b_run_cell2sentence.py
# Purpose: Run Cell2Sentence AI model for cell embeddings and predictions.
# This script should be run in the 'impact_sc' conda environment.

import anndata
import pandas as pd
import scanpy as sc
import cell2sentence as cs
import numpy as np
import random
from collections import Counter
import os
import multiprocessing
import time
import sys

# --- Robust System Detection ---
def detect_pytorch_and_device():
    """Detect PyTorch availability and optimal device"""
    try:
        import torch
        torch_available = True
        torch_version = torch.__version__
        print(f"✅ PyTorch {torch_version} detected")
        
        # Device detection with detailed info
        if torch.cuda.is_available():
            device = "cuda"
            gpu_name = torch.cuda.get_device_name(0)
            gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1024**3
            print(f"🚀 CUDA GPU detected: {gpu_name} ({gpu_memory:.1f}GB)")
        elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
            device = "mps"
            print(f"🍎 Apple Silicon MPS acceleration available")
        else:
            device = "cpu"
            print(f"💻 Using CPU processing")
            
        return torch, device, torch_available
        
    except ImportError:
        print("⚠️ PyTorch not found - falling back to CPU-only mode")
        return None, "cpu", False

def setup_multiprocessing():
    """Setup optimal multiprocessing configuration"""
    try:
        n_cores = multiprocessing.cpu_count()
        optimal_cores = min(n_cores, 8)  # Cap at 8 cores
        
        # Set environment variables for parallel processing
        os.environ["OMP_NUM_THREADS"] = str(optimal_cores)
        os.environ["MKL_NUM_THREADS"] = str(optimal_cores)
        
        print(f"💾 System: {n_cores} cores detected, using {optimal_cores} cores")
        return optimal_cores, True
        
    except Exception as e:
        print(f"⚠️ Multiprocessing setup failed: {e}")
        print("📝 Falling back to single-core processing")
        return 1, False

def progress_callback(step, total, message="Processing"):
    """Simple progress indicator"""
    percent = (step / total) * 100
    bar_length = 30
    filled_length = int(bar_length * step // total)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    print(f'\r{message}: |{bar}| {percent:.1f}% ({step}/{total})', end='', flush=True)
    if step == total:
        print()  # New line when complete

# --- Initialize System ---
torch, DEVICE, torch_available = detect_pytorch_and_device()
N_CORES, multiprocessing_available = setup_multiprocessing()

# Configure PyTorch if available
if torch_available:
    torch.set_num_threads(N_CORES)
    
# --- Parameters & Setup ---
SEED = 1234
N_GENES = 200  # Increased from 5 to 200 for better representation
random.seed(SEED)
np.random.seed(SEED)
if torch_available:
    torch.manual_seed(SEED)

print(f"\n🔬 Analysis Parameters:")
print(f"   • Number of genes: {N_GENES}")
print(f"   • Random seed: {SEED}")
print(f"   • Multiprocessing: {'✅ Enabled' if multiprocessing_available else '❌ Disabled'}")

# Get file paths and species from environment variables set by the R script
h5ad_file_path = os.getenv("H5AD_FILE_PATH")
c2s_output_dir = os.getenv("C2S_OUTPUT_DIR")
c2s_embeddings_csv = os.getenv("C2S_EMBEDDINGS_CSV")
c2s_predicted_csv = os.getenv("C2S_PREDICTED_CSV")
current_species = os.getenv("SPECIES_FOR_C2S", "human") 

# Get Cell2Sentence model path or name from environment variable
# Fallback to a default Hugging Face model name if the environment variable is not set
default_hf_model_name = "james-y-u/C2S-Pythia-410m-cell-type-prediction"
c2s_model_identifier = os.getenv("C2S_MODEL_PATH_OR_NAME", default_hf_model_name)


if not all([h5ad_file_path, c2s_output_dir, c2s_embeddings_csv, c2s_predicted_csv]):
    print("Error: One or more environment variables for file paths are not set.")
    print("Please ensure H5AD_FILE_PATH, C2S_OUTPUT_DIR, C2S_EMBEDDINGS_CSV, C2S_PREDICTED_CSV are set.")
    exit(1)

if not c2s_model_identifier: 
    print(f"Error: Cell2Sentence model identifier (C2S_MODEL_PATH_OR_NAME) is not set and no default is available.")
    exit(1)


if not os.path.exists(h5ad_file_path):
    print(f"Error: H5AD file not found at {h5ad_file_path}")
    exit(1)

if not os.path.exists(c2s_output_dir):
    try:
        os.makedirs(c2s_output_dir, exist_ok=True)
        print(f"Created C2S output directory: {c2s_output_dir}")
    except Exception as e:
        print(f"Error creating C2S output directory {c2s_output_dir}: {e}")
        exit(1)

print(f"--- Running Cell2Sentence for species: {current_species} ---")
print(f"Loading h5ad from: {h5ad_file_path}")

adata = anndata.read_h5ad(h5ad_file_path)
adata.obs.head()
adata.obs["cell_type"] = adata.obs["seurat_clusters"].astype(str)
Counter(adata.obs["cell_type"])

adata.obs["tissue"] = adata.obs.get("tissue", "unknown_tissue")
adata.obs["batch_condition"] = adata.obs.get("orig.ident", "unknown_batch") 
adata.obs["organism"] = adata.obs.get("organism", current_species)
adata.obs["sex"] = adata.obs.get("sex", "unknown_sex") 
obs_cols = ["cell_type", "tissue", "batch_condition", "organism", "sex"]

print("Converting AnnData to Arrow format for Cell2Sentence...")
start_time = time.time()
try:
    # Try with parallel processing first
    if multiprocessing_available:
        print("🔄 Attempting parallel data conversion...")
        arrow_ds, vocab = cs.CSData.adata_to_arrow(
            adata, 
            random_state=SEED, 
            sentence_delimiter=' ', 
            label_col_names=obs_cols,
            num_workers=min(N_CORES, 4)  # Parallelize data conversion
        )
    else:
        raise Exception("Multiprocessing not available, using single-core")
    print("✅ Arrow dataset created (parallel mode).")
except Exception as e:
    print(f"⚠️ Parallel conversion failed: {e}")
    print("🔄 Falling back to single-core conversion...")
    try:
        arrow_ds, vocab = cs.CSData.adata_to_arrow(
            adata, 
            random_state=SEED, 
            sentence_delimiter=' ', 
            label_col_names=obs_cols
        )
        print("✅ Arrow dataset created (single-core mode).")
    except Exception as e2:
        print(f"❌ Data conversion failed: {e2}")
        exit(1)

conversion_time = time.time() - start_time
print(f"⏱️ Data conversion completed in {conversion_time:.2f} seconds")

csdata_save_dir = os.path.join(c2s_output_dir, "csdata_arrow")
if not os.path.exists(csdata_save_dir):
    os.makedirs(csdata_save_dir, exist_ok=True)

print("Creating CSData object...")
start_time = time.time()
try:
    # Try with parallel processing first
    if multiprocessing_available:
        print("🔄 Attempting parallel CSData creation...")
        csdata = cs.CSData.csdata_from_arrow(
            arrow_dataset=arrow_ds, 
            vocabulary=vocab, 
            save_dir=csdata_save_dir, 
            save_name="cell_embedding", 
            dataset_backend="arrow",
            num_workers=min(N_CORES, 4)  # Parallelize CSData creation
        )
    else:
        raise Exception("Multiprocessing not available, using single-core")
    print("✅ CSData object created (parallel mode).")
except Exception as e:
    print(f"⚠️ Parallel CSData creation failed: {e}")
    print("🔄 Falling back to single-core CSData creation...")
    try:
        csdata = cs.CSData.csdata_from_arrow(
            arrow_dataset=arrow_ds, 
            vocabulary=vocab, 
            save_dir=csdata_save_dir, 
            save_name="cell_embedding", 
            dataset_backend="arrow"
        )
        print("✅ CSData object created (single-core mode).")
    except Exception as e2:
        print(f"❌ CSData creation failed: {e2}")
        exit(1)

csdata_time = time.time() - start_time
print(f"⏱️ CSData creation completed in {csdata_time:.2f} seconds")

# Use the model identifier from environment variable (or default)
model_path_to_load = c2s_model_identifier
# Define a directory where models downloaded from Hugging Face might be cached or where local models are expected.
# The `save_dir` in CSModel might be used as `cache_dir` by the underlying `transformers` library.
# If `model_path_to_load` is an absolute path to a local model, `cache_dir` might be ignored by `transformers`.
# It's good practice to have a dedicated cache/save directory for models.
model_cache_or_save_dir = os.path.join(c2s_output_dir, "c2s_model_cache") 
if not os.path.exists(model_cache_or_save_dir):
    os.makedirs(model_cache_or_save_dir, exist_ok=True)

print(f"Attempting to load Cell2Sentence model using identifier: {model_path_to_load}")
print(f"Model cache/save directory set to: {model_cache_or_save_dir}")

try:
    # If model_path_to_load is a local path, cell2sentence should handle it.
    # If it's a Hugging Face name, it will try to download.
    # The `save_dir` parameter in `cs.CSModel` might be used for different purposes
    # depending on whether the model is loaded locally or from Hugging Face.
    # For Hugging Face models, it can act as a `cache_dir`.
    csmodel = cs.CSModel(
        model_name_or_path=model_path_to_load, 
        save_dir=model_cache_or_save_dir, 
        save_name="c2s_model_instance" # This name is more for if you were fine-tuning and saving a new model.
    )
    print("Cell2Sentence model loaded successfully.")
except Exception as e:
    print(f"Error loading Cell2Sentence model using identifier '{model_path_to_load}': {e}")
    print("Please ensure you have a stable internet connection if it's a Hugging Face model name.")
    print("If it's a local path, ensure the path is correct and the model files are intact (e.g., contains 'config.json', 'pytorch_model.bin').")
    print("If the model was recently made private or moved from Hugging Face, this could also cause issues.")
    exit(1)

print("\n🧬 Embedding cells...")
start_time = time.time()
try:
    # Try optimized embedding with device and parallel processing
    if torch_available and (DEVICE != "cpu" or multiprocessing_available):
        print(f"🔄 Attempting optimized embedding (device: {DEVICE}, cores: {N_CORES})...")
        embedded_cells = cs.tasks.embed_cells(
            csdata=csdata, 
            csmodel=csmodel, 
            n_genes=N_GENES,
            device=DEVICE,          # Use optimal device (GPU/MPS/CPU)
            num_workers=min(N_CORES, 4) if multiprocessing_available else 1
        )
        print("✅ Cell embeddings generated (optimized mode).")
    else:
        raise Exception("No optimization available, using basic mode")
        
except Exception as e:
    print(f"⚠️ Optimized embedding failed: {e}")
    print("🔄 Falling back to basic embedding...")
    try:
        embedded_cells = cs.tasks.embed_cells(
            csdata=csdata, 
            csmodel=csmodel, 
            n_genes=N_GENES
        )
        print("✅ Cell embeddings generated (basic mode).")
    except Exception as fallback_e:
        print(f"❌ All embedding methods failed: {fallback_e}")
        exit(1)

# Store embeddings in AnnData
adata.obsm["c2s_cell_embeddings"] = embedded_cells
embedding_time = time.time() - start_time
n_cells = adata.shape[0]
print(f"⏱️ Cell embedding completed in {embedding_time:.2f} seconds")
print(f"📊 Processed {n_cells} cells ({n_cells/embedding_time:.1f} cells/sec)")

print("\n🔬 Predicting cell types...")
start_time = time.time()
try:
    # Try optimized prediction with device and parallel processing
    if torch_available and (DEVICE != "cpu" or multiprocessing_available):
        print(f"🔄 Attempting optimized prediction (device: {DEVICE}, cores: {N_CORES})...")
        predicted_df = cs.tasks.predict_cell_types_of_data(
            csdata=csdata, 
            csmodel=csmodel, 
            n_genes=N_GENES,
            device=DEVICE,          # Use optimal device
            num_workers=min(N_CORES, 4) if multiprocessing_available else 1
        )
        print("✅ Cell type prediction complete (optimized mode).")
    else:
        raise Exception("No optimization available, using basic mode")
        
except Exception as e:
    print(f"⚠️ Optimized prediction failed: {e}")
    print("🔄 Falling back to basic prediction...")
    try:
        predicted_df = cs.tasks.predict_cell_types_of_data(
            csdata=csdata, 
            csmodel=csmodel, 
            n_genes=N_GENES
        )
        print("✅ Cell type prediction complete (basic mode).")
    except Exception as fallback_e:
        print(f"⚠️ Basic prediction failed: {fallback_e}")
        print("📝 Creating empty DataFrame - prediction results unavailable")
        predicted_df = pd.DataFrame()

prediction_time = time.time() - start_time
print(f"⏱️ Cell type prediction completed in {prediction_time:.2f} seconds")

# Summary statistics
total_time = conversion_time + csdata_time + embedding_time + prediction_time
print(f"\n📈 Performance Summary:")
print(f"   • Data conversion: {conversion_time:.2f}s")
print(f"   • CSData creation: {csdata_time:.2f}s") 
print(f"   • Cell embedding: {embedding_time:.2f}s")
print(f"   • Type prediction: {prediction_time:.2f}s")
print(f"   • Total time: {total_time:.2f}s")
print(f"   • Processing rate: {n_cells/total_time:.1f} cells/sec")

try:
    embedded_df = pd.DataFrame(embedded_cells, index=adata.obs_names)
    print(f"Saving C2S cell embeddings to: {c2s_embeddings_csv}")
    embedded_df.to_csv(c2s_embeddings_csv)
except Exception as e:
    print(f"Error saving cell embeddings: {e}")

if not predicted_df.empty:
    try:
        print(f"Saving C2S predicted cell types to: {c2s_predicted_csv}")
        predicted_df.to_csv(c2s_predicted_csv, index=False)
    except Exception as e:
        print(f"Error saving predicted cell types: {e}")
else:
    print("Skipping saving of predicted cell types as prediction failed or produced no results.")

print("--- Cell2Sentence Python script finished. ---")
