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

# --- Parameters & Setup ---
SEED = 1234
random.seed(SEED)
np.random.seed(SEED)

# Get file paths and species from environment variables set by the R script
h5ad_file_path = os.getenv("H5AD_FILE_PATH")
c2s_output_dir = os.getenv("C2S_OUTPUT_DIR")
c2s_embeddings_csv = os.getenv("C2S_EMBEDDINGS_CSV")
c2s_predicted_csv = os.getenv("C2S_PREDICTED_CSV")
current_species = os.getenv("SPECIES_FOR_C2S", "human") # Default to human if not set

if not all([h5ad_file_path, c2s_output_dir, c2s_embeddings_csv, c2s_predicted_csv]):
    print("Error: One or more environment variables for file paths are not set.")
    print("Please ensure H5AD_FILE_PATH, C2S_OUTPUT_DIR, C2S_EMBEDDINGS_CSV, C2S_PREDICTED_CSV are set.")
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

try:
    adata = anndata.read_h5ad(h5ad_file_path)
    print("AnnData object loaded successfully.")
except Exception as e:
    print(f"Error loading AnnData object: {e}")
    exit(1)

# Ensure the obs column from R (cell_type_low_res) is present
label_col_name = 'cell_type_low_res' 
if label_col_name not in adata.obs.columns:
    if 'cell_type' in adata.obs.columns: # Fallback to 'cell_type' if 'cell_type_low_res' is missing
        label_col_name = 'cell_type'
        print(f"Warning: '{label_col_name}' not found, using 'cell_type' as label column.")
    else:
        print(f"Error: Column '{label_col_name}' (or 'cell_type') not found in adata.obs.")
        print(f"Available columns: {adata.obs.columns.tolist()}")
        exit(1)
print(f"Cell type distribution from '{label_col_name}': {Counter(adata.obs[label_col_name])}")

# Add dummy annotations if required by CSData or use actual metadata if available
adata.obs["tissue"] = adata.obs.get("tissue", "unknown_tissue") 
adata.obs["batch_condition"] = adata.obs.get("orig.ident", "unknown_batch") # Use orig.ident if present
adata.obs["organism"] = adata.obs.get("organism", current_species)
adata.obs["sex"] = adata.obs.get("sex", "unknown_sex") 
obs_cols = [label_col_name, "tissue", "batch_condition", "organism", "sex"]

print("Converting AnnData to Arrow format for Cell2Sentence...")
try:
    # Cell2Sentence might expect raw counts for its internal gene selection.
    # Ensure adata.X is appropriate (e.g., counts or log-normalized counts based on C2S docs)
    # If SaveH5Seurat saved logcounts to adata.X, this might be fine.
    arrow_ds, vocab = cs.CSData.adata_to_arrow(
        adata, 
        random_state=SEED, 
        sentence_delimiter=' ', 
        label_col_names=obs_cols
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
        dataset_backend="arrow"
    )
    print("CSData object created.")
except Exception as e:
    print(f"Error creating CSData object: {e}")
    exit(1)

# Load pretrained model
model_path = "C2S-Pythia-410m-cell-type-prediction" # Pretrained model name
model_save_dir = os.path.join(c2s_output_dir, "csmodel")
if not os.path.exists(model_save_dir):
    os.makedirs(model_save_dir, exist_ok=True)

try:
    csmodel = cs.CSModel(
        model_name_or_path=model_path, 
        save_dir=model_save_dir, 
        save_name="cell_embedding_prediction"
    )
    print("Cell2Sentence model loaded.")
except Exception as e:
    print(f"Error loading Cell2Sentence model: {e}")
    exit(1)

# Embed cells
print("Embedding cells...")
try:
    # n_genes=200 is a parameter for the embedding process
    embedded_cells = cs.tasks.embed_cells(csdata=csdata, csmodel=csmodel, n_genes=200)
    adata.obsm["c2s_cell_embeddings"] = embedded_cells
    print("Cell embeddings generated and stored in AnnData.")
except Exception as e:
    print(f"Error embedding cells: {e}")
    exit(1)

# Predict cell types (if the model is for prediction)
print("Predicting cell types...")
try:
    predicted_df = cs.tasks.predict_cell_types_of_data(csdata=csdata, csmodel=csmodel, n_genes=200)
    print("Cell type prediction complete.")
except Exception as e:
    print(f"Error predicting cell types: {e}")
    # Continue to save embeddings even if prediction fails
    predicted_df = pd.DataFrame() # Empty dataframe

# Save results
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
