#!/usr/bin/env python3

import json
import os
import sys # Import sys for sys.exit()
from typing import List, Dict, Any
import subprocess # Ensure subprocess is imported

# Configuration for Ollama (if used for more advanced interaction in future)
OLLAMA_MODEL_NAME_DEFAULT = "gemma3:12b-it-qat" # You can change the default Ollama model if needed
OLLAMA_BASE_URL_DEFAULT = "http://localhost:11434" # Base URL for Ollama API

def normalize_path(path_str: str) -> str:
    """Converts all backslashes to forward slashes in a path string."""
    if path_str is None:
        return None
    return path_str.replace("\\", "/")

def ask_question(prompt: str, default_value: str = None, choices: List[str] = None) -> str:
    """Asks a question and returns the user's input or default.
    If choices are provided, validates input against choices."""
    while True:
        full_prompt = prompt
        if default_value is not None: # Check if default_value is actually provided
            full_prompt += f" (default: {default_value})"
        
        response = input(f"{full_prompt}: ").strip()
        
        if not response and default_value is not None: # User pressed Enter and default exists
            response = default_value
        
        if choices:
            choice_map = {c.lower(): c for c in choices} # Case-insensitive matching for choices
            if response.lower() in choice_map:
                return choice_map[response.lower()] # Return the original casing of the choice
            else:
                print(f"Invalid choice. Please select from: {', '.join(choices)}")
        elif response: # If not using choices, any non-empty response is valid
            return response
        # If response is empty AND no default_value is provided, and it's not a choice-based question
        # (implicitly mandatory if no default and not choices)
        elif default_value is None and not response: 
             print("This field cannot be empty. Please provide a value.")
        else: # Covers case where response is empty but default_value was None (and not choices)
            # This case should ideally be caught by the above if it's truly mandatory.
            # If it reaches here, it means an empty response for a non-mandatory field without default.
            return response # Return empty string if allowed

def ask_for_paths(prompt: str, allow_multiple: bool = False, is_optional: bool = False, optional_default_skip: str = "skip", ensure_file: bool = False, ensure_dir: bool = False, default_path_suggestion: str = None) -> List[str]:
    """Asks for file/directory paths, validates them, and normalizes to forward slashes."""
    paths = []
    
    # Normalize default_path_suggestion if provided
    normalized_default_suggestion = normalize_path(default_path_suggestion)

    prompt_to_display = prompt
    if normalized_default_suggestion:
        prompt_to_display += f" (suggested: {normalized_default_suggestion})"

    actual_prompt_text = prompt_to_display 
    if is_optional:
        actual_prompt_text += f" (enter '{optional_default_skip}' to skip or if not applicable)"
    
    first_path = True
    while True:
        if not allow_multiple and not first_path: 
            break # Exit if only one path is allowed and we've got it

        current_prompt = actual_prompt_text if first_path else "Enter additional path"
        if first_path and is_optional and not allow_multiple : 
             pass # No change to prompt for single optional path
        elif not first_path and allow_multiple and is_optional: # For multiple optional paths, allow skipping further additions
            current_prompt += f" (or type '{optional_default_skip}' to finish adding paths)"

        path_input = input(f"{current_prompt}: ").strip()

        if is_optional and path_input.lower() == optional_default_skip:
            if allow_multiple and not first_path: # If adding multiple and user types skip for an additional path
                break 
            return [] # Skip this entire path question

        if not path_input and normalized_default_suggestion and first_path: # Use suggestion if input is empty
            path_input = normalized_default_suggestion
            print(f"Using suggested path: {path_input}")

        # Handle empty input for optional paths without a suggestion
        if is_optional and not path_input and not normalized_default_suggestion: 
            print(f"No path entered. Assuming skip for optional path: {prompt}")
            return []

        # Handle empty input for mandatory paths
        if not path_input and not normalized_default_suggestion and not is_optional:
            print("This path cannot be empty. Please provide a value.")
            continue # Re-ask for the path

        # Normalize the user's input path
        normalized_input_path = normalize_path(path_input)
        abs_path = normalize_path(os.path.abspath(normalized_input_path)) # Get absolute path and normalize it too
        
        path_exists = os.path.exists(abs_path)
        error_message = ""

        if not path_exists:
            if ensure_file or ensure_dir: # Error only if we need to ensure it's a file/dir
                 error_message = f"Error: Path '{normalized_input_path}' (resolved to '{abs_path}') does not exist."
        elif ensure_file and not os.path.isfile(abs_path):
            error_message = f"Error: Path '{normalized_input_path}' (resolved to '{abs_path}') is not a file."
        elif ensure_dir and not os.path.isdir(abs_path):
            error_message = f"Error: Path '{normalized_input_path}' (resolved to '{abs_path}') is not a directory."
        
        if not error_message:
            # If path validation is required (ensure_file/ensure_dir) or if path exists, store abs_path.
            # Otherwise (e.g. model name), store the (normalized) user input.
            if ensure_file or ensure_dir or path_exists:
                paths.append(abs_path)
            else:
                paths.append(normalized_input_path) 
            
            if not allow_multiple:
                break # Done if only one path needed
            else: # Ask if user wants to add more paths
                first_path = False 
                another = ask_question("Do you want to add another path? (yes/no)", "no", choices=["yes", "no"]).lower()
                if another != 'yes':
                    break
        else: # Path validation failed
            if is_optional:
                error_message += f" Please re-enter or type '{optional_default_skip}' to skip."
            else:
                error_message += " Please re-enter."
            print(error_message)
    return paths


def select_modules() -> List[str]:
    """Allows user to select which modules to run."""
    all_modules = {
        "1": "01_data_processing",
        "2a": "02a_harmony_c2s_prep",
        "2b": "02b_c2s",
        "2c": "02c_load_c2s_result",
        "3": "03_cell_type_annotation",
        "4a": "04a_basic_visualization",
        "4b": "04b_DE_gsea",
        "4c": "04c_decoupler (Requires R >= 4.3.0)",
        "4d": "04d_ucell_scores",
        "4e": "04e_pseudotime", 
        "4f": "04f_query_projection",
        "4g": "04g_card (Requires R >= 4.3.0)"
    }
    print("\nAvailable IMPACT-sc Modules:")
    for key, name in all_modules.items():
        print(f"  {key}: {name}")

    selected_keys_str = ask_question("Enter the keys of modules to run, separated by commas (e.g., 1,2a,2b,2c,3,4a)")
    selected_keys = [key.strip().lower() for key in selected_keys_str.split(',') if key.strip()]

    selected_modules_list = []
    valid_keys_found = False
    for key in selected_keys:
        # Extract the module name from the display text
        module_name = all_modules.get(key, "").split(" (")[0]
        if module_name:
            selected_modules_list.append(module_name)
            valid_keys_found = True
        else:
            print(f"Warning: Module key '{key}' is not valid and will be ignored.")

    if not valid_keys_found :
        print("No valid modules selected or input was empty. Defaulting to Module 1 only.")
        return ["01_data_processing"]
    return selected_modules_list

def ask_for_dotplot_genes() -> List[Dict[str, Any]]:
    """Asks the user for gene groups for DotPlot."""
    print("\n--- DotPlot Gene Configuration ---")
    dotplot_groups = []
    while True:
        add_group = ask_question("Do you want to add a gene group for DotPlot? (yes/no)", "yes", choices=["yes", "no"]).lower()
        if add_group != 'yes':
            break

        group_name = ask_question("Enter the name for this gene group (e.g., B_cell_markers)")
        genes_str = ask_question(f"Enter comma-separated genes for '{group_name}' (e.g., MS4A1,CD79A)")
        genes_list = [gene.strip() for gene in genes_str.split(',') if gene.strip()]

        if group_name and genes_list:
            dotplot_groups.append({"name": group_name, "genes": genes_list})
            print(f"Group '{group_name}' with genes {genes_list} added.")
        else:
            print("Group name or gene list was empty. Group not added.")

    if not dotplot_groups:
        print("No gene groups provided for DotPlot. The R script will skip DotPlot if no genes are configured via environment variables.")
    return dotplot_groups

def check_downloaded_data() -> Dict[str, str]:
    """Check for downloaded data and return normalized paths if found."""
    current_dir = normalize_path(os.path.dirname(os.path.abspath(__file__)))
    data_paths = {
        "demo_data": None,
        "c2s_model": None,
        "default_reference_data": None
    }
    
    demo_path_base = normalize_path(os.path.join(current_dir, "data", "demo", "filtered_gene_bc_matrices"))
    if os.path.exists(demo_path_base):
        hg19_path = normalize_path(os.path.join(demo_path_base, "hg19"))
        if os.path.exists(hg19_path):
            matrix_file = normalize_path(os.path.join(hg19_path, "matrix.mtx"))
            genes_file = normalize_path(os.path.join(hg19_path, "genes.tsv"))
            barcodes_file = normalize_path(os.path.join(hg19_path, "barcodes.tsv"))
            
            genes_file_gz = genes_file + ".gz"
            barcodes_file_gz = barcodes_file + ".gz"

            if os.path.exists(matrix_file) and \
               (os.path.exists(genes_file) or os.path.exists(genes_file_gz)) and \
               (os.path.exists(barcodes_file) or os.path.exists(barcodes_file_gz)):
                data_paths["demo_data"] = hg19_path
                print(f"Found demo data: {hg19_path}")
            else:
                print(f"Demo directory found but missing some 10x files (matrix.mtx, genes.tsv/genes.tsv.gz, barcodes.tsv/barcodes.tsv.gz): {hg19_path}")
        else:
            data_paths["demo_data"] = demo_path_base 
            print(f"Found demo data (base directory): {demo_path_base}")
    else:
        print(f"Demo data directory not found at: {demo_path_base}")
    
    models_dir = normalize_path(os.path.join(current_dir, "data", "models"))
    if os.path.exists(models_dir):
        import glob
        model_files = glob.glob(normalize_path(os.path.join(models_dir, "**", "*C2S*")), recursive=True) + \
                     glob.glob(normalize_path(os.path.join(models_dir, "**", "*Pythia*")), recursive=True)
        if model_files:
            data_paths["c2s_model"] = "vandijklab/C2S-Pythia-410m-cell-type-prediction"
            print(f"Found potential Cell2Sentence model cache in: {models_dir} (will use HuggingFace name for robustness)")
        else:
            print(f"Models directory exists but no obvious cached C2S/Pythia files found: {models_dir}")
    
    ref_path = normalize_path(os.path.join(current_dir, "data", "reference", "bmcite_demo.rds"))
    if os.path.exists(ref_path):
        data_paths["default_reference_data"] = ref_path
        print(f"Found default reference data: {ref_path}")
    else:
        print(f"Default reference data (bmcite_demo.rds) not found at: {ref_path}")
    
    return data_paths

def ask_demo_mode() -> bool:
    """Ask if user wants to use demo mode."""
    print("\n" + "="*60)
    print("IMPACT-sc Setup Mode Selection")
    print("="*60)
    
    downloaded_data = check_downloaded_data()
    has_demo_data = downloaded_data.get("demo_data") is not None
    
    if has_demo_data:
        print("Demo data detected! You can run in demo mode.")
        print("Demo data includes:")
        print("   - PBMC3k single-cell dataset (10x Genomics)")
        print("   - Pre-configured parameters for quick testing")
        print("")
        
        mode_choice = ask_question(
            "Choose setup mode:\n  [demo] - Use demo data and auto-configure\n  [custom] - Configure your own data",
            "demo",
            choices=["demo", "custom"]
        ).lower()
        
        return mode_choice == "demo"
    else:
        print("No demo data found. Run './download_data.sh' first to enable demo mode.")
        print("Proceeding with custom setup...")
        return False

def main():
    print("--- Welcome to IMPACT-sc Interactive Setup ---")
    
    downloaded_data = check_downloaded_data()
    use_demo_mode = ask_demo_mode()
    
    print("\n--- Rscript Executable Path ---")
    default_rscript_path_suggestion = auto_detect_rscript()
    
    rscript_exe_paths = ask_for_paths(
        "Enter the full path to your Rscript executable (recommended: .../R-4.2.3/bin/x64/Rscript.exe)",
        ensure_file=True,
        default_path_suggestion=default_rscript_path_suggestion
    )
    if not rscript_exe_paths: 
        print("CRITICAL ERROR: Rscript executable path not provided or invalid. Exiting.")
        sys.exit(1)
    rscript_executable_path = rscript_exe_paths[0]
    print(f"Rscript executable path set to: {rscript_executable_path}")

    if use_demo_mode:
        return setup_demo_mode(downloaded_data, rscript_executable_path)
    else:
        return setup_custom_mode(downloaded_data, rscript_executable_path)


def setup_demo_mode(downloaded_data: Dict[str, str], rscript_executable_path: str):
    """Setup demo mode using downloaded data and pre-configured parameters."""
    print("\nSetting up DEMO MODE...")
    print("Using pre-configured parameters for PBMC3k dataset.")
    
    include_c2s = ask_question(
        "Include Cell2Sentence (C2S) analysis? This provides AI-powered cell type prediction but takes longer",
        "no",
        choices=["yes", "no"]
    ).lower() == "yes"

    params: Dict[str, Any] = {}
    pipeline_base_dir = normalize_path(os.path.dirname(os.path.abspath(__file__)))
    
    params["input_r_scripts_dir"] = normalize_path(os.path.join(pipeline_base_dir, "scripts_AI"))
    params["input_python_scripts_dir"] = normalize_path(os.path.join(pipeline_base_dir, "scripts_AI"))
    params["rscript_executable_path"] = rscript_executable_path
    params["species"] = "human"
    params["output_directory"] = normalize_path(os.path.abspath("demo_output"))
    
    if downloaded_data.get("demo_data"):
        params["input_data_paths"] = [downloaded_data["demo_data"]]
        print(f"Using demo data: {downloaded_data['demo_data']}")
    else:
        print("Demo data not found! This should have been caught by ask_demo_mode.")
        print("Please run './download_data.sh' to download demo data or configure custom paths.")
        return False 
    
    if include_c2s:
        print("Including Cell2Sentence analysis (slower but more comprehensive)")
        params["selected_modules"] = [
            "01_data_processing", "02a_harmony_c2s_prep", "02b_c2s",
            "02c_load_c2s_result", "03_cell_type_annotation", "04a_basic_visualization"
        ]
        params["h5ad_path_for_c2s"] = normalize_path(os.path.join(params["output_directory"], "02_module2_for_c2s.h5ad"))
        params["c2s_model_path_or_name"] = downloaded_data.get("c2s_model", "vandijklab/C2S-Pythia-410m-cell-type-prediction")
        if downloaded_data.get("c2s_model"): print(f"Using potentially cached Cell2Sentence model (name: {params['c2s_model_path_or_name']})")
        else: print("No cached model hint found. Will download C2S model on first use.")
    else:
        print("Skipping Cell2Sentence for faster demo (using traditional annotation only)")
        params["selected_modules"] = ["01_data_processing", "02a_harmony_c2s_prep", "03_cell_type_annotation", "04a_basic_visualization"]
        params["h5ad_path_for_c2s"] = None
        params["c2s_model_path_or_name"] = None
    
    print("Setting SingleR reference to the downloaded local file for demo.")
    default_ref_path = normalize_path(os.path.join(pipeline_base_dir, "data", "reference", "bmcite_demo.rds"))
    
    if os.path.exists(default_ref_path):
        params["local_singler_ref_path"] = default_ref_path
        params["local_singler_ref_label_col"] = "celltype.l1"
        print(f"Using local reference for demo: {default_ref_path}")
    else:
        print(f"CRITICAL ERROR for Demo Mode: The default reference file was not found at the expected path: {default_ref_path}")
        print("Please run './download_data.sh' to download the required data.")
        return False

    params["final_cell_type_source"] = "auto"
    params["reduction_method"] = "umap"
    params["featureplot_genes"] = "CD3D,CD14,MS4A1,FCGR3A,LYZ,PPBP"
    params["dotplot_gene_groups"] = [
        {"name": "T_cell_markers", "genes": ["CD3D", "CD3E", "CD8A", "CD4"]},
        {"name": "B_cell_markers", "genes": ["MS4A1", "CD79A", "CD79B"]},
        {"name": "Myeloid_markers", "genes": ["CD14", "LYZ", "FCGR3A", "CST3"]}
    ]
    params["de_gsea_plot_gene"] = "CD3D"
    params["msigdb_category"] = "H"
    params["ucell_plot_pathway_name"] = ""
    params["conditional_paths"] = {"query_rds_path": None, "query_species": None, "palantir_start_cell": None}
    params["spatial_data_rds_path"] = None # Not used in demo
    params["ollama_model_name"] = OLLAMA_MODEL_NAME_DEFAULT
    params["ollama_base_url"] = OLLAMA_BASE_URL_DEFAULT
    
    save_params(params)
    print("\nDemo mode setup complete!")
    return True

def auto_detect_rscript() -> str:
    """Auto-detect Rscript executable and return normalized path."""
    import subprocess
    import platform
    
    rscript_path = None
    try:
        if platform.system() == "Windows":
            result = subprocess.run(['where', 'Rscript.exe'], capture_output=True, text=True, check=False, timeout=5)
        else:
            result = subprocess.run(['which', 'Rscript'], capture_output=True, text=True, check=False, timeout=5)
        
        if result.returncode == 0 and result.stdout.strip():
            rscript_path = result.stdout.strip().split('\n')[0]
    except (FileNotFoundError, subprocess.TimeoutExpired, subprocess.CalledProcessError):
        pass 
    
    if rscript_path and os.path.exists(normalize_path(rscript_path)) and os.path.isfile(normalize_path(rscript_path)):
        return normalize_path(rscript_path)

    fallback_paths = []
    if platform.system() == "Darwin":
        fallback_paths = ["/opt/homebrew/bin/Rscript", "/usr/local/bin/Rscript", "/Library/Frameworks/R.framework/Resources/bin/Rscript"]
    elif platform.system() == "Windows":
        program_files = os.environ.get("ProgramFiles", "C:\\Program Files")
        versions = ["R-4.2.3", "R-4.3.1", "R-4.3.0", "R-4.4.2", "R-4.1.3", "R-4.0.5"]
        for ver in versions:
            fallback_paths.append(os.path.join(program_files, "R", ver, "bin", "x64", "Rscript.exe"))
            fallback_paths.append(os.path.join(program_files, "R", ver, "bin", "Rscript.exe"))
    else: 
        fallback_paths = ["/usr/local/bin/Rscript", "/usr/bin/Rscript"]
    
    for path_attempt in fallback_paths:
        normalized_attempt = normalize_path(path_attempt)
        if os.path.exists(normalized_attempt) and os.path.isfile(normalized_attempt):
            return normalized_attempt
    
    return normalize_path("Rscript")


def setup_custom_mode(downloaded_data: Dict[str, str], rscript_executable_path: str):
    """Setup custom mode with user inputs, using downloaded data as defaults, and normalize paths."""
    print("\nSetting up CUSTOM MODE...")
    print("Please provide the following information:")

    params: Dict[str, Any] = {}
    params["rscript_executable_path"] = rscript_executable_path

    print("\n--- Script Locations ---")
    pipeline_base_dir = normalize_path(os.path.dirname(os.path.abspath(__file__)))
    default_scripts_dir_suggestion = normalize_path(os.path.join(pipeline_base_dir, "scripts_AI"))
    if not os.path.isdir(default_scripts_dir_suggestion):
        default_scripts_dir_suggestion_alt = normalize_path(os.path.join(pipeline_base_dir, "..", "scripts_AI"))
        if os.path.isdir(default_scripts_dir_suggestion_alt):
            default_scripts_dir_suggestion = default_scripts_dir_suggestion_alt
        else: 
            default_scripts_dir_suggestion = None

    scripts_dir_paths = ask_for_paths(
        "Enter the full path to the directory containing R and Python module scripts",
        ensure_dir=True,
        default_path_suggestion=default_scripts_dir_suggestion
    )
    if not scripts_dir_paths: 
        print("CRITICAL ERROR: Scripts directory not provided or invalid. Exiting.")
        sys.exit(1)
    params["input_r_scripts_dir"] = scripts_dir_paths[0] 
    params["input_python_scripts_dir"] = scripts_dir_paths[0]

    print("\n--- Input Data ---")
    demo_suggestion = downloaded_data.get("demo_data")
    if demo_suggestion: print(f"Tip: Demo data available at: {demo_suggestion}")
    
    input_data_paths_val = ask_for_paths(
        "Enter the full path to your primary input scRNA-seq data file(s) (e.g., feature-barcode matrix directory, or an .RDS file like 'ori.RDS')",
        allow_multiple=True, is_optional=False, default_path_suggestion=demo_suggestion
    )
    if not input_data_paths_val:
        print("CRITICAL ERROR: Input data path(s) not provided. Exiting.")
        sys.exit(1)
    params["input_data_paths"] = input_data_paths_val

    params["species"] = ask_question("Enter the species ('human' or 'mouse')", "human", choices=["human", "mouse"]).lower()

    print("\n--- Output Configuration ---")
    output_dir_input = ask_question("Enter the full path for your desired output/results folder", "demo_output")
    params["output_directory"] = normalize_path(os.path.abspath(output_dir_input))


    print("\n--- Module Selection ---")
    params["selected_modules"] = select_modules()

    params["h5ad_path_for_c2s"] = None 
    params["c2s_model_path_or_name"] = None 
    if "02b_c2s" in params["selected_modules"]:
        print("\n--- Cell2Sentence (Module 02b) Specific Input ---")
        
        h5ad_c2s_paths = ask_for_paths(
            "Enter the full path to the H5AD file for Cell2Sentence input (e.g., F:/R_PROJECT/impact/02_module2_for_c2s.h5ad)",
            is_optional=False, 
            ensure_file=True
        )

        if h5ad_c2s_paths: 
            params["h5ad_path_for_c2s"] = h5ad_c2s_paths[0]
        else:
            print("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no valid H5AD input path was provided.")
            sys.exit(1)

        default_c2s_model = downloaded_data.get("c2s_model", "vandijklab/C2S-Pythia-410m-cell-type-prediction")
        if downloaded_data.get("c2s_model"): print(f"Tip: Cached model name available: {default_c2s_model}")
        
        c2s_model_input_list = ask_for_paths(
            f"Enter the Cell2Sentence model path (if local, e.g., F:/path/to/model_folder) or Hugging Face model name",
            is_optional=False, ensure_dir=False, default_path_suggestion=default_c2s_model
        )
        if c2s_model_input_list:
            c2s_model_input = c2s_model_input_list[0]
            if (os.path.sep in c2s_model_input or "/" in c2s_model_input) and os.path.isdir(c2s_model_input):
                params["c2s_model_path_or_name"] = c2s_model_input 
                print(f"Using local Cell2Sentence model from directory: {params['c2s_model_path_or_name']}")
            else:
                params["c2s_model_path_or_name"] = c2s_model_input 
                print(f"Using Cell2Sentence model (name or non-validated path): {params['c2s_model_path_or_name']}")
        else:
            print("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no Cell2Sentence model path/name was provided.")
            sys.exit(1)

    # --- [MODIFICATION] ---
    # The following block that asks for CollecTRI and PROGENy CSV paths has been removed
    # because the R script will now download them automatically.
    if "04c_decoupler" in params["selected_modules"]:
        print("\n--- DecoupleR Network Information (Module 04c) ---")
        print("Note: This module requires R version 4.3.0 or higher to run.")
        print("CollecTRI and PROGENy networks will be downloaded automatically based on the selected species.")

    params["spatial_data_rds_path"] = None
    if "04g_card" in params["selected_modules"]:
        print("\n--- CARD (Module 04g) Specific Input ---")
        print("Note: This module requires R version 4.3.0 or higher to run.")
        spatial_rds_paths = ask_for_paths(
            "Enter the full path to the spatial data RDS file (for CARD)",
            allow_multiple=False, is_optional=False, ensure_file=True
        )
        if spatial_rds_paths:
            params["spatial_data_rds_path"] = spatial_rds_paths[0]
        else:
            print("CRITICAL ERROR: Module 04g selected, but no spatial data RDS path provided. Exiting.")
            sys.exit(1)

    if "03_cell_type_annotation" in params["selected_modules"]:
        print("\n--- Final Cell Type Source (Module 03) ---")
        params["final_cell_type_source"] = ask_question(
            "Which annotation source should be used for the final 'cell_type' column? (auto, seurat, c2s, singler)",
            "auto", choices=["auto", "seurat", "c2s", "singler"]
        ).lower()
    else: params["final_cell_type_source"] = "auto"

    if "04a_basic_visualization" in params["selected_modules"]:
        print("\n--- Basic Visualization (Module 04a) Specific Inputs ---")
        params["reduction_method"] = ask_question(
            "Enter the preferred reduction method for plotting",
            default_value="umap",
            choices=["umap_c2s", "umap", "harmony", "pca"]
        )
        params["featureplot_genes"] = ask_question("Enter comma-separated genes for FeaturePlot. Leave empty to skip.", "")
        params["dotplot_gene_groups"] = ask_for_dotplot_genes()
    else:
        params["reduction_method"] = "umap"
        params["featureplot_genes"] = ""
        params["dotplot_gene_groups"] = []

    if "04b_DE_gsea" in params["selected_modules"]:
        print("\n--- Differential Expression & GSEA (Module 04b) Specific Inputs ---")
        params["de_gsea_plot_gene"] = ask_question("Enter a single gene for violin plot in module 4b. Leave empty to skip.", "")
    else: params["de_gsea_plot_gene"] = ""

    if "04d_ucell_scores" in params["selected_modules"]:
        print("\n--- UCell Gene Scores (Module 04d) Specific Inputs ---")
        params["msigdb_category"] = ask_question(f"Enter MSigDB category for UCell (e.g., H, C2, C5).", "H").upper()
        params["ucell_plot_pathway_name"] = ask_question("Enter specific pathway name to plot for UCell. Leave empty for first.", "")
    else:
        params["msigdb_category"] = "H"
        params["ucell_plot_pathway_name"] = ""

    if "03_cell_type_annotation" in params["selected_modules"]:
        print("\n--- SingleR Reference Configuration (Module 03) ---")
        print("This module requires a local SingleR reference data file in .RDS format.")

        ref_suggestion = downloaded_data.get("default_reference_data")
        if ref_suggestion:
             print(f"Tip: A downloaded reference is available at: {ref_suggestion}")

        local_ref_paths = ask_for_paths(
            f"Enter the full path to your local SingleR reference RDS file for {params['species']}",
            allow_multiple=False,
            is_optional=False,
            ensure_file=True,
            default_path_suggestion=ref_suggestion
        )
        if local_ref_paths:
            params["local_singler_ref_path"] = local_ref_paths[0]
        else:
            print("CRITICAL ERROR: Module 03 was selected, but no valid reference file path was provided. Exiting.")
            sys.exit(1)

        default_label_suggestion = "label.main"
        if params["local_singler_ref_path"] and "bmcite_demo.rds" in params["local_singler_ref_path"]:
            default_label_suggestion = "celltype.l1"

        label_col = ask_question(
            "Enter the name of the metadata column in your reference file that contains cell type labels",
            default_value=default_label_suggestion
        )
        params["local_singler_ref_label_col"] = label_col
        print(f"Using local reference: {params['local_singler_ref_path']} with label column: '{label_col}'")
    else:
        params["local_singler_ref_path"] = None
        params["local_singler_ref_label_col"] = None

    params["conditional_paths"] = {"query_rds_path": None, "query_species": None, "palantir_start_cell": None}
    if "04e_pseudotime" in params["selected_modules"]:
        print("\n--- Pseudotime Analysis (Module 04e) Specific Input ---")
        params["conditional_paths"]["palantir_start_cell"] = ask_question("Enter start cell for Palantir (e.g., AAACCCAAGTCGAAGG-1)", "first_cell_barcode")

    if "04f_query_projection" in params["selected_modules"]:
        print("\n--- Query Dataset Projection (Module 04f) Specific Input ---")
        query_rds_paths = ask_for_paths("Enter full path to query.RDS file", False, False, ensure_file=True)
        if query_rds_paths: params["conditional_paths"]["query_rds_path"] = query_rds_paths[0]
        else:
            print(f"CRITICAL ERROR: Module 04f selected, but no query.RDS path provided.")
            sys.exit(1)
        params["conditional_paths"]["query_species"] = ask_question("Enter species of query dataset ('human' or 'mouse')", params["species"], ["human", "mouse"]).lower()

    params["ollama_model_name"] = OLLAMA_MODEL_NAME_DEFAULT
    params["ollama_base_url"] = OLLAMA_BASE_URL_DEFAULT
    
    save_params(params)
    print("\nCustom mode setup complete!")
    return True


def save_params(params: Dict[str, Any]) -> bool:
    """Save parameters to JSON file. Expects all paths in params to be pre-normalized."""
    try:
        output_dir = params["output_directory"]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")

        params_path = normalize_path(os.path.join(output_dir, "impact_sc_params.json"))
        with open(params_path, 'w', encoding='utf-8') as f:
            json.dump(params, f, indent=4)
        print(f"\nParameters saved to: {params_path}")

        print(f"\n--- Setup Complete ---")
        print(f"Next steps:")
        print(f"1. Ensure all R and Python dependencies have been correctly installed.")
        print(f"2. Activate the 'impact_sc' conda environment: conda activate impact_sc")
        print(f"3. Run the pipeline using: python run_impact_sc_pipeline.py {params_path}")
        
        return True
    except IOError as e:
        print(f"Error: Could not write parameters file. {e}")
        return False
    except Exception as e: 
        print(f"An unexpected error occurred while saving parameters: {e}")
        return False

if __name__ == "__main__":
    if main():
        sys.exit(0)
    else:
        print("Setup did not complete successfully.")
        sys.exit(1)
