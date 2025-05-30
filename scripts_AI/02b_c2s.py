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
import torch
import multiprocessing

# --- Performance Optimization Setup ---
# Set number of CPU cores to use
N_CORES = min(multiprocessing.cpu_count(), 8)  # Use up to 8 cores
os.environ["OMP_NUM_THREADS"] = str(N_CORES)
os.environ["MKL_NUM_THREADS"] = str(N_CORES)
torch.set_num_threads(N_CORES)

# Check for GPU availability
DEVICE = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
print(f"🚀 Performance Settings:")
print(f"   • Using {N_CORES} CPU cores")
print(f"   • Device: {DEVICE}")
print(f"   • PyTorch threads: {torch.get_num_threads()}")

# --- Parameters & Setup ---
SEED = 1234
BATCH_SIZE = 32  # Process cells in batches for better memory management
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

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
try:
    arrow_ds, vocab = cs.CSData.adata_to_arrow(
        adata, 
        random_state=SEED, 
        sentence_delimiter=' ', 
        label_col_names=obs_cols,
        num_workers=min(N_CORES, 4)  # Parallelize data conversion
    )
    print("Arrow dataset created.")
except Exception as e:
    print(f"Error converting AnnData to Arrow format: {e}")
    exit(1)

csdata_save_dir = os.path.join(c2s_output_dir, "csdata_arrow")
if not os.path.exists(csdata_save_dir):
    os.makedirs(csdata_save_dir, exist_ok=True)

try:
    csdata = cs.CSData.csdata_from_arrow(
        arrow_dataset=arrow_ds, 
        vocabulary=vocab, 
        save_dir=csdata_save_dir, 
        save_name="cell_embedding", 
        dataset_backend="arrow",
        num_workers=min(N_CORES, 4)  # Parallelize CSData creation
    )
    print("CSData object created.")
except Exception as e:
    print(f"Error creating CSData object: {e}")
    exit(1)

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

print("Embedding cells...")
try:
    # Optimized embedding with batch processing and device specification
    embedded_cells = cs.tasks.embed_cells(
        csdata=csdata, 
        csmodel=csmodel, 
        n_genes=5,
        batch_size=BATCH_SIZE,  # Process in batches
        device=DEVICE,          # Use optimal device (GPU/MPS/CPU)
        num_workers=min(N_CORES, 4)  # Parallel data loading
    )
    adata.obsm["c2s_cell_embeddings"] = embedded_cells
    print("Cell embeddings generated and stored in AnnData.")
except Exception as e:
    print(f"Error embedding cells: {e}")
    # Fallback to basic embedding if optimized version fails
    try:
        print("Attempting fallback embedding without optimization...")
        embedded_cells = cs.tasks.embed_cells(csdata=csdata, csmodel=csmodel, n_genes=5)
        adata.obsm["c2s_cell_embeddings"] = embedded_cells
        print("Fallback embedding successful.")
    except Exception as fallback_e:
        print(f"Fallback embedding also failed: {fallback_e}")
        exit(1)

print("Predicting cell types...")
try:
    # Optimized prediction with batch processing
    predicted_df = cs.tasks.predict_cell_types_of_data(
        csdata=csdata, 
        csmodel=csmodel, 
        n_genes=5,
        batch_size=BATCH_SIZE,  # Process in batches
        device=DEVICE,          # Use optimal device
        num_workers=min(N_CORES, 4)  # Parallel processing
    )
    print("Cell type prediction complete.")
except Exception as e:
    print(f"Error predicting cell types: {e}")
    # Fallback to basic prediction
    try:
        print("Attempting fallback prediction without optimization...")
        predicted_df = cs.tasks.predict_cell_types_of_data(csdata=csdata, csmodel=csmodel, n_genes=5)
        print("Fallback prediction successful.")
    except Exception as fallback_e:
        print(f"Fallback prediction also failed: {fallback_e}")
        predicted_df = pd.DataFrame()

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
