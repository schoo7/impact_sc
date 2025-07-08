#!/usr/bin/env python3

import json
import os
import sys
from typing import List, Dict, Any
import subprocess
import urllib.request
import urllib.error

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
        if default_value is not None:
            full_prompt += f" (default: {default_value})"
        
        response = input(f"{full_prompt}: ").strip()
        
        if not response and default_value is not None:
            response = default_value
        
        if choices:
            choice_map = {c.lower(): c for c in choices}
            if "," in response:
                responses = [r.strip().lower() for r in response.split(',')]
                valid_responses = [choice_map[r] for r in responses if r in choice_map]
                if len(valid_responses) == len(responses):
                    return ",".join(valid_responses)
                else:
                    print(f"Invalid choice detected. Please select from: {', '.join(choices)}")
                    print("You can provide multiple values separated by commas.")
            elif response.lower() in choice_map:
                return choice_map[response.lower()]
            else:
                print(f"Invalid choice. Please select from: {', '.join(choices)}")
        elif response:
            return response
        elif default_value is None and not response: 
             print("This field cannot be empty. Please provide a value.")
        else:
            return response

def ask_for_paths(prompt: str, allow_multiple: bool = False, is_optional: bool = False, optional_default_skip: str = "skip", ensure_file: bool = False, ensure_dir: bool = False, default_path_suggestion: str = None) -> List[str]:
    """Asks for file/directory paths, validates them, and normalizes to forward slashes."""
    paths = []
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
            break

        current_prompt = actual_prompt_text if first_path else "Enter additional path"
        if not first_path and allow_multiple and is_optional:
            current_prompt += f" (or type '{optional_default_skip}' to finish adding paths)"

        path_input = input(f"{current_prompt}: ").strip()

        if is_optional and path_input.lower() == optional_default_skip:
            if allow_multiple and not first_path:
                break 
            return []

        if not path_input and normalized_default_suggestion and first_path:
            path_input = normalized_default_suggestion
            print(f"Using suggested path: {path_input}")

        if is_optional and not path_input and not normalized_default_suggestion: 
            print(f"No path entered. Assuming skip for optional path: {prompt}")
            return []

        if not path_input and not normalized_default_suggestion and not is_optional:
            print("This path cannot be empty. Please provide a value.")
            continue

        normalized_input_path = normalize_path(path_input)
        # For non-existence checks, use the user's input path directly for flexibility.
        # For existence checks, resolve to an absolute path.
        path_to_check = normalize_path(os.path.abspath(normalized_input_path)) if (ensure_dir or ensure_file) else normalized_input_path

        path_exists = os.path.exists(path_to_check)
        error_message = ""

        if not path_exists:
            # Only error out if we are explicitly told the path must exist.
            if ensure_file or ensure_dir:
                 error_message = f"Error: Path '{normalized_input_path}' (resolved to '{path_to_check}') does not exist."
        elif ensure_file and not os.path.isfile(path_to_check):
            error_message = f"Error: Path '{normalized_input_path}' (resolved to '{path_to_check}') is not a file."
        elif ensure_dir and not os.path.isdir(path_to_check):
            error_message = f"Error: Path '{normalized_input_path}' (resolved to '{path_to_check}') is not a directory."
        
        if not error_message:
            # If path validation is required and passes, use the absolute path.
            # Otherwise, use the user-provided normalized path (e.g., for output files).
            path_to_add = path_to_check if (ensure_dir or ensure_file) else normalized_input_path
            paths.append(path_to_add)
            
            if not allow_multiple:
                break
            else:
                first_path = False 
                another = ask_question("Do you want to add another path? (yes/no)", "no", choices=["yes", "no"]).lower()
                if another != 'yes':
                    break
        else:
            if is_optional:
                error_message += f" Please re-enter or type '{optional_default_skip}' to skip."
            else:
                error_message += " Please re-enter."
            print(error_message)
    return paths

def select_modules() -> List[str]:
    """Allows user to select which modules to run."""
    all_modules = {
        "1": "01_data_processing", "2a": "02a_harmony_c2s_prep", "2b": "02b_c2s",
        "2c": "02c_load_c2s_result", "3": "03_cell_type_annotation", "4a": "04a_basic_visualization",
        "4b": "04b_DE_gsea", "4c": "04c_decoupler (Requires R >= 4.3.0)", "4d": "04d_ucell_scores",
        "4e": "04e_pseudotime", "4f": "04f_query_projection", "4g": "04g_card (Requires R >= 4.3.0)",
        "4h": "04h_cell_chat"
    }
    print("\nAvailable IMPACT-sc Modules:")
    for key, name in all_modules.items():
        print(f"  {key}: {name}")
    selected_keys_str = ask_question("Enter the keys of modules to run, separated by commas (e.g., 1,2a,2b,2c,3,4a,4h)")
    selected_keys = [key.strip().lower() for key in selected_keys_str.split(',') if key.strip()]
    selected_modules_list = []
    valid_keys_found = False
    for key in selected_keys:
        module_name = all_modules.get(key, "").split(" (")[0]
        if module_name:
            selected_modules_list.append(module_name)
            valid_keys_found = True
        else:
            print(f"Warning: Module key '{key}' is not valid and will be ignored.")
    if not valid_keys_found:
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
        print("No gene groups provided for DotPlot.")
    return dotplot_groups

def check_downloaded_data() -> Dict[str, str]:
    """Check for downloaded data and return normalized paths if found."""
    current_dir = normalize_path(os.path.dirname(os.path.abspath(__file__)))
    data_paths = {"demo_data": None, "c2s_model": None, "default_reference_data": None}
    demo_path_base = normalize_path(os.path.join(current_dir, "data", "demo", "filtered_gene_bc_matrices"))
    if os.path.exists(demo_path_base):
        hg19_path = normalize_path(os.path.join(demo_path_base, "hg19"))
        if os.path.exists(hg19_path):
            matrix_file = normalize_path(os.path.join(hg19_path, "matrix.mtx"))
            genes_file = normalize_path(os.path.join(hg19_path, "genes.tsv"))
            barcodes_file = normalize_path(os.path.join(hg19_path, "barcodes.tsv"))
            genes_file_gz, barcodes_file_gz = genes_file + ".gz", barcodes_file + ".gz"
            if os.path.exists(matrix_file) and (os.path.exists(genes_file) or os.path.exists(genes_file_gz)) and (os.path.exists(barcodes_file) or os.path.exists(barcodes_file_gz)):
                data_paths["demo_data"] = hg19_path
                print(f"Found demo data: {hg19_path}")
            else:
                print(f"Demo directory found but missing some 10x files in: {hg19_path}")
        else:
            data_paths["demo_data"] = demo_path_base
            print(f"Found demo data (base directory): {demo_path_base}")
    else:
        print(f"Demo data directory not found at: {demo_path_base}")
    models_dir = normalize_path(os.path.join(current_dir, "data", "models"))
    if os.path.exists(models_dir):
        import glob
        model_files = glob.glob(normalize_path(os.path.join(models_dir, "**", "*C2S*")), recursive=True) + glob.glob(normalize_path(os.path.join(models_dir, "**", "*Pythia*")), recursive=True)
        if model_files:
            data_paths["c2s_model"] = "vandijklab/C2S-Pythia-410m-cell-type-prediction"
            print(f"Found potential Cell2Sentence model cache in: {models_dir}")
        else:
            print(f"Models directory exists but no obvious cached C2S/Pythia files found: {models_dir}")
    ref_path = normalize_path(os.path.join(current_dir, "data", "reference", "bmcite_demo.rds"))
    if os.path.exists(ref_path):
        data_paths["default_reference_data"] = ref_path
        print(f"Found default reference data: {ref_path}")
    else:
        print(f"Default reference data (bmcite_demo.rds) not found at: {ref_path}")
    return data_paths

def run_ai_assistant_flow() -> str:
    """Handles the AI assistant user interaction before custom setup."""
    print("\n--- AI Assistant Mode ---")
    print("This mode helps configure the pipeline by understanding your analysis goals.")
    user_description = ask_question(
        "\nPlease describe your dataset and the analysis you would like to perform.\nFor example: 'I have 10x Genomics human PBMC data and I want to find cell type markers and see how different clusters communicate.'"
    )
    return user_description

def get_ai_configuration(user_description: str, ollama_model_name: str, ollama_base_url: str) -> Dict[str, Any]:
    """Connects to Ollama to get a suggested pipeline configuration."""
    prompt = f"""
You are an expert bioinformatician configuring a single-cell analysis pipeline. Based on the user's request, you must decide which modules to run and what parameters to use.

The user's request is: "{user_description}"

Here are the available modules:
- "01_data_processing": Basic filtering and normalization.
- "02a_harmony_c2s_prep": Data integration using Harmony.
- "02b_c2s": Cell type prediction using Cell2Sentence (requires 02a).
- "02c_load_c2s_result": Load C2S results back into the Seurat object (requires 02b).
- "03_cell_type_annotation": Cell type annotation using SingleR or ceLLama.
- "04a_basic_visualization": Generate UMAPs, feature plots, and dot plots.
- "04b_DE_gsea": Differential expression and gene set enrichment analysis.
- "04c_decoupler": Transcription factor activity analysis.
- "04d_ucell_scores": Calculate gene signature scores.
- "04e_pseudotime": Pseudotime analysis with Palantir.
- "04f_query_projection": Project a query dataset onto a reference.
- "04g_card": Spatial deconvolution.
- "04h_cell_chat": Cell-cell communication analysis.

Here are the parameters you can set:
- "remove_doublets" (boolean): Recommended for most datasets.
- "regress_cell_cycle" (boolean): Use if cell cycle is a major source of variation.
- "qc_min_nfeature_rna" (integer): Typical values: 200-500.
- "qc_max_nfeature_rna" (integer): Typical values: 4000-8000.
- "qc_max_percent_mt" (integer): Typical values: 5-20.
- "pca_dims" (integer): Default is 50.
- "cluster_resolution" (float): Higher for more clusters. Typical: 0.1-1.2.
- "dims_for_clustering" (integer): Default is 50.
- "annotation_method" (string): "singler" or "cellama". Choose "singler" for standard reference-based annotation.
- "final_cell_type_source" (string): "auto", "seurat", "c2s", "singler", or "cellama". Usually same as annotation_method.
- "featureplot_genes" (string): Comma-separated list of genes. Suggest common markers if possible.
- "liana_method" (string): e.g., "logfc", "natmi". Default to "logfc" if unsure.
- "cellchat_source_groups" (string): Comma-separated cluster IDs for LIANA source.
- "cellchat_target_groups" (string): Comma-separated cluster IDs for LIANA target.

Based on the user's request, provide a JSON object with your recommended configuration.
Only output the raw JSON object. Do not include any other text, explanations, or markdown formatting like ```json.

Example:
{{
  "selected_modules": ["01_data_processing", "03_cell_type_annotation", "04a_basic_visualization"],
  "remove_doublets": true,
  "qc_min_nfeature_rna": 200,
  "qc_max_percent_mt": 15,
  "cluster_resolution": 0.5,
  "annotation_method": "singler"
}}
"""
    api_url = f"{ollama_base_url}/api/generate"
    payload = {
        "model": ollama_model_name,
        "prompt": prompt,
        "stream": False,
        "format": "json"
    }
    headers = {"Content-Type": "application/json"}
    try:
        req = urllib.request.Request(api_url, data=json.dumps(payload).encode('utf-8'), headers=headers, method='POST')
        with urllib.request.urlopen(req) as response:
            response_body = response.read().decode('utf-8')
            response_json = json.loads(response_body)
            # The actual generated JSON content is in the 'response' key
            ai_config_str = response_json.get("response", "{}")
            return json.loads(ai_config_str)
    except urllib.error.URLError as e:
        print(f"\nError connecting to Ollama at {api_url}: {e.reason}")
        print("Please ensure the Ollama server is running and accessible.")
        return None
    except json.JSONDecodeError:
        print("\nError: Failed to decode the JSON response from the AI model.")
        print(f"Received: {ai_config_str}")
        return None
    except Exception as e:
        print(f"\nAn unexpected error occurred while communicating with the AI: {e}")
        return None

def select_setup_mode(downloaded_data: Dict[str, str]) -> str:
    """Ask user to select the setup mode: demo, custom, or AI-assisted."""
    print("\n" + "="*60 + "\nIMPACT-sc Setup Mode Selection\n" + "="*60)
    has_demo_data = downloaded_data.get("demo_data") is not None
    mode_choices, prompt_lines = ["custom", "ai"], ["Choose setup mode:", "  [custom] - Manually configure all parameters.", "  [ai]     - (New!) Describe your project to get help with configuration."]
    default_choice = "custom"
    if has_demo_data:
        print("Demo data detected! You can run in demo mode.")
        mode_choices.insert(0, "demo")
        prompt_lines.insert(1, "  [demo]   - Use pre-downloaded demo data for a quick test.")
        default_choice = "demo"
    else:
        print("No demo data found. Run './download_data.sh' first to enable demo mode.")
    return ask_question("\n".join(prompt_lines), default_choice, choices=mode_choices).lower()

def main():
    """Main function to drive the interactive setup."""
    print("--- Welcome to IMPACT-sc Interactive Setup ---")
    downloaded_data = check_downloaded_data()
    mode = select_setup_mode(downloaded_data)
    
    print("\n--- Rscript Executable Path ---")
    default_rscript_path_suggestion = auto_detect_rscript()
    rscript_exe_paths = ask_for_paths(
        "Enter the full path to your Rscript executable (e.g., .../R-4.3.1/bin/Rscript)",
        ensure_file=True, default_path_suggestion=default_rscript_path_suggestion
    )
    if not rscript_exe_paths:
        print("CRITICAL ERROR: Rscript executable path not provided or invalid. Exiting.")
        sys.exit(1)
    rscript_executable_path = rscript_exe_paths[0]
    print(f"Rscript executable path set to: {rscript_executable_path}")

    if mode == 'demo':
        return setup_demo_mode(downloaded_data, rscript_executable_path)
    elif mode == 'ai':
        return setup_ai_mode(downloaded_data, rscript_executable_path)
    else:
        return setup_custom_mode(downloaded_data, rscript_executable_path)

def setup_ai_mode(downloaded_data: Dict[str, str], rscript_executable_path: str) -> bool:
    """Drives the AI-assisted setup flow."""
    user_description = run_ai_assistant_flow()
    
    ollama_model = ask_question("Enter the Ollama model to use", OLLAMA_MODEL_NAME_DEFAULT)
    ollama_url = ask_question("Enter the Ollama base URL", OLLAMA_BASE_URL_DEFAULT)

    print("\nConnecting to AI Assistant to generate pipeline configuration... Please wait.")
    ai_params = get_ai_configuration(user_description, ollama_model, ollama_url)

    if not ai_params:
        print("\nAI-assisted setup failed. Could not retrieve a valid configuration. Exiting.")
        return False

    print("\nAI has suggested the following configuration:")
    print(json.dumps(ai_params, indent=2))
    
    if ask_question("Do you want to proceed with this configuration?", "yes", ["yes", "no"]).lower() == 'no':
        print("Setup aborted by user.")
        return False

    params = ai_params
    params["rscript_executable_path"] = rscript_executable_path
    params["ollama_model_name"] = ollama_model
    params["ollama_base_url"] = ollama_url
    
    print("\n--- Please provide the necessary paths for the pipeline ---")
    pipeline_base_dir = normalize_path(os.path.dirname(os.path.abspath(__file__)))
    default_scripts_dir = normalize_path(os.path.join(pipeline_base_dir, "scripts_AI"))
    scripts_dir_paths = ask_for_paths("Enter path to the module scripts directory", ensure_dir=True, default_path_suggestion=default_scripts_dir)
    if not scripts_dir_paths: sys.exit("CRITICAL ERROR: Scripts directory not provided.")
    params["input_r_scripts_dir"] = params["input_python_scripts_dir"] = scripts_dir_paths[0]

    input_data_paths = ask_for_paths("Enter path to your input scRNA-seq data (matrix directory or .RDS)", allow_multiple=True)
    if not input_data_paths: sys.exit("CRITICAL ERROR: Input data path(s) not provided.")
    params["input_data_paths"] = input_data_paths

    params["species"] = ask_question("Enter species ('human' or 'mouse')", "human", ["human", "mouse"]).lower()
    output_dir = ask_question("Enter path for the output folder", "demo_output")
    params["output_directory"] = normalize_path(os.path.abspath(output_dir))

    selected_modules = params.get("selected_modules", [])
    if "02b_c2s" in selected_modules:
        print("\n--- AI-selected pipeline requires Cell2Sentence inputs ---")
        # The H5AD file is an output of module 02a. We ask the user to confirm its future path.
        default_h5ad_path_suggestion = normalize_path(os.path.join(params.get("output_directory", "."), "02_module2_for_c2s.h5ad"))
        h5ad_path = ask_for_paths(
            "Confirm or enter the path for the H5AD file to be generated for C2S",
            ensure_file=False,  # Corrected: File doesn't exist yet.
            default_path_suggestion=default_h5ad_path_suggestion,
            allow_multiple=False
        )
        if not h5ad_path: 
            sys.exit("CRITICAL ERROR: H5AD path for C2S not provided.")
        params["h5ad_path_for_c2s"] = h5ad_path[0]

        # Ask for the C2S model path or name.
        default_c2s_model = downloaded_data.get("c2s_model", "vandijklab/C2S-Pythia-410m-cell-type-prediction")
        c2s_model_input_list = ask_for_paths(
            "Confirm or enter the C2S model path or Hugging Face name",
            is_optional=False,
            ensure_dir=False,
            ensure_file=False,
            default_path_suggestion=default_c2s_model,
            allow_multiple=False
        )
        if c2s_model_input_list:
            params["c2s_model_path_or_name"] = c2s_model_input_list[0]
        else:
            sys.exit("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no C2S model path/name was provided.")

    if "03_cell_type_annotation" in selected_modules and params.get("annotation_method") == "singler":
        print("\n--- AI requires SingleR reference file ---")
        ref_suggestion = downloaded_data.get("default_reference_data")
        ref_path = ask_for_paths("Enter path to your local SingleR reference RDS file", ensure_file=True, default_path_suggestion=ref_suggestion)
        if not ref_path: sys.exit("CRITICAL ERROR: SingleR reference not provided.")
        params["local_singler_ref_path"] = ref_path[0]
        params["local_singler_ref_label_col"] = ask_question("Enter the label column name in the reference", "celltype.l1")
    
    # Set other conditional params to None if not set by AI, to avoid key errors
    params.setdefault("conditional_paths", {})
    params.setdefault("spatial_data_rds_path", None)

    save_params(params)
    return True


def setup_demo_mode(downloaded_data: Dict[str, str], rscript_executable_path: str):
    """Setup demo mode using downloaded data and pre-configured parameters."""
    print("\nSetting up DEMO MODE...")
    print("Using pre-configured parameters for PBMC3k dataset.")

    # Per user request, C2S is disabled and Seurat method is used.
    # The 'include_c2s' question is removed.
    print("Skipping Cell2Sentence analysis and using Seurat method for annotation.")
    
    include_cellchat = ask_question(
        "Include Cell-Cell Communication (CellChat/LIANA) analysis (Module 04h)?",
        "yes",
        choices=["yes", "no"]
    ).lower() == "yes"


    params: Dict[str, Any] = {}
    pipeline_base_dir = normalize_path(os.path.dirname(os.path.abspath(__file__)))
    
    params["input_r_scripts_dir"] = normalize_path(os.path.join(pipeline_base_dir, "scripts_AI"))
    params["input_python_scripts_dir"] = normalize_path(os.path.join(pipeline_base_dir, "scripts_AI"))
    params["rscript_executable_path"] = rscript_executable_path
    params["species"] = "human"
    params["output_directory"] = normalize_path(os.path.abspath("demo_output"))
    
    # --- Data Processing Options for DEMO ---
    params["remove_doublets"] = False
    params["regress_cell_cycle"] = False
    params["qc_min_nfeature_rna"] = 200
    params["qc_max_nfeature_rna"] = 6000
    params["qc_max_percent_mt"] = 10
    params["pca_dims"] = 50
    params["cluster_resolution"] = 0.1
    params["dims_for_clustering"] = 50
    print("Demo mode will use default QC, PCA, and Clustering parameters for speed.")

    if downloaded_data.get("demo_data"):
        params["input_data_paths"] = [downloaded_data["demo_data"]]
        print(f"Using demo data: {downloaded_data['demo_data']}")
    else:
        print("Demo data not found! Please run './download_data.sh' to download demo data.")
        return False 
    
    # Per user request, modules are fixed to 01, 02a, 03, 04a. CellChat is optional.
    base_modules = ["01_data_processing", "02a_harmony_c2s_prep", "03_cell_type_annotation", "04a_basic_visualization"]
    
    # C2S is always skipped in this modified demo mode.
    params["h5ad_path_for_c2s"] = None
    params["c2s_model_path_or_name"] = None
    
    if include_cellchat:
        print("Including LIANA analysis.")
        base_modules.append("04h_cell_chat")
        params["cellchat_source_groups"] = "0,1,2,3"
        params["cellchat_target_groups"] = "4,5,6,7,8"
        params["liana_method"] = "logfc" # Set default method for demo mode
        print(f"Using demo source groups: {params['cellchat_source_groups']}")
        print(f"Using demo target groups: {params['cellchat_target_groups']}")
        print(f"Using LIANA method: {params['liana_method']}")

    params["selected_modules"] = base_modules

    # Per user request, using Seurat method, so SingleR reference is not needed.
    print("Using 'seurat' for annotation method. No external reference file is needed.")
    params["local_singler_ref_path"] = None
    params["local_singler_ref_label_col"] = None

    params["annotation_method"] = "seurat"
    params["final_cell_type_source"] = "seurat"
    params["cellama_temperature"] = 0.0
    params["ollama_model_name"] = OLLAMA_MODEL_NAME_DEFAULT

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
    scripts_dir_paths = ask_for_paths("Enter the full path to the directory containing R and Python module scripts", ensure_dir=True, default_path_suggestion=default_scripts_dir_suggestion)
    if not scripts_dir_paths: sys.exit("CRITICAL ERROR: Scripts directory not provided or invalid. Exiting.")
    params["input_r_scripts_dir"] = scripts_dir_paths[0] 
    params["input_python_scripts_dir"] = scripts_dir_paths[0]
    print("\n--- Input Data ---")
    demo_suggestion = downloaded_data.get("demo_data")
    if demo_suggestion: print(f"Tip: Demo data available at: {demo_suggestion}")
    input_data_paths_val = ask_for_paths("Enter the full path to your primary input scRNA-seq data file(s)", allow_multiple=True, is_optional=False, default_path_suggestion=demo_suggestion)
    if not input_data_paths_val: sys.exit("CRITICAL ERROR: Input data path(s) not provided. Exiting.")
    params["input_data_paths"] = input_data_paths_val
    params["species"] = ask_question("Enter the species ('human' or 'mouse')", "human", choices=["human", "mouse"]).lower()
    print("\n--- Output Configuration ---")
    output_dir_input = ask_question("Enter the full path for your desired output/results folder", "demo_output")
    params["output_directory"] = normalize_path(os.path.abspath(output_dir_input))
    print("\n--- Data Processing & Analysis Parameters (Modules 01 & 02a) ---")
    params["remove_doublets"] = ask_question("Remove potential doublets using scDblFinder?", "no", choices=["yes", "no"]).lower() == "yes"
    params["regress_cell_cycle"] = ask_question("Regress out cell cycle effects (S/G2M scores)?", "no", choices=["yes", "no"]).lower() == "yes"
    print("\n--- Quality Control (QC) Parameters ---")
    params["qc_min_nfeature_rna"] = int(ask_question("Enter minimum nFeature_RNA (genes per cell)", "200"))
    params["qc_max_nfeature_rna"] = int(ask_question("Enter maximum nFeature_RNA (genes per cell)", "6000"))
    params["qc_max_percent_mt"] = int(ask_question("Enter maximum mitochondrial gene percentage", "10"))
    print("\n--- PCA & Clustering Parameters ---")
    params["pca_dims"] = int(ask_question("Enter number of principal components (PCs) for PCA", "50"))
    params["cluster_resolution"] = float(ask_question("Enter clustering resolution for FindClusters", "0.1"))
    params["dims_for_clustering"] = int(ask_question("Enter number of dimensions for clustering (e.g., 50 for UMAP, 1024 for C2S)", "50"))
    print("\n--- Module Selection ---")
    params["selected_modules"] = select_modules()
    if "03_cell_type_annotation" in params["selected_modules"]:
        print("\n--- Cell Type Annotation Method (Module 03) ---")
        annotation_choice = ask_question("Choose annotation source for the final 'cell_type' column(seurat, c2s, singler, cellama)", "auto", choices=["auto", "seurat", "c2s", "singler", "cellama"]).lower()
        params["annotation_method"] = params["final_cell_type_source"] = annotation_choice
        print(f"Annotation method set to: '{annotation_choice}'.")
        params["cellama_temperature"], params["ollama_model_name"] = 0.0, OLLAMA_MODEL_NAME_DEFAULT
        if params["annotation_method"] == "cellama":
            print("\nThe 'temperature' parameter controls the randomness of the ceLLama model's output.")
            try:
                params["cellama_temperature"] = float(ask_question("Enter ceLLama temperature value", "0.0"))
            except ValueError:
                params["cellama_temperature"] = 0.0
            params["ollama_model_name"] = ask_question("Enter the Ollama model name for ceLLama", OLLAMA_MODEL_NAME_DEFAULT)
    else:
        params["annotation_method"], params["final_cell_type_source"], params["cellama_temperature"], params["ollama_model_name"] = "singler", "auto", 0.0, OLLAMA_MODEL_NAME_DEFAULT
    
    params["h5ad_path_for_c2s"], params["c2s_model_path_or_name"] = None, None
    if "02b_c2s" in params["selected_modules"]:
        print("\n--- Cell2Sentence (Module 02b) Specific Input ---")
        # The H5AD file is an output of module 02a. We ask the user to define its future path.
        # This path will be used as the input for module 02b.
        default_h5ad_path_suggestion = normalize_path(os.path.join(params.get("output_directory", "."), "02_module2_for_c2s.h5ad"))
        h5ad_c2s_paths = ask_for_paths(
            "Enter the full path where the H5AD file from Module 02a should be saved",
            is_optional=False,
            ensure_file=False,  # Corrected: This file does not exist yet.
            default_path_suggestion=default_h5ad_path_suggestion
        )
        if h5ad_c2s_paths:
            params["h5ad_path_for_c2s"] = h5ad_c2s_paths[0]
        else:
            sys.exit("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no H5AD file path was provided.")

        default_c2s_model = downloaded_data.get("c2s_model", "vandijklab/C2S-Pythia-410m-cell-type-prediction")
        if downloaded_data.get("c2s_model"):
            print(f"Tip: Cached model name available: {default_c2s_model}")
        
        # We don't check for file/dir existence as it could be a Hugging Face model name.
        c2s_model_input_list = ask_for_paths(
            "Enter the local C2S model path or a Hugging Face model name",
            is_optional=False,
            ensure_dir=False,
            ensure_file=False,
            default_path_suggestion=default_c2s_model
        )
        if c2s_model_input_list:
            params["c2s_model_path_or_name"] = c2s_model_input_list[0]
        else:
            sys.exit("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no C2S model path/name was provided.")

    if "04g_card" in params["selected_modules"]:
        print("\n--- CARD (Module 04g) Specific Input ---")
        spatial_rds_paths = ask_for_paths("Enter path to the spatial data RDS file for CARD", allow_multiple=False, is_optional=False, ensure_file=True)
        if spatial_rds_paths: params["spatial_data_rds_path"] = spatial_rds_paths[0]
        else: sys.exit("CRITICAL ERROR: Module 04g selected, but no spatial data RDS path provided.")
    else: params["spatial_data_rds_path"] = None
    if "04a_basic_visualization" in params["selected_modules"]:
        print("\n--- Basic Visualization (Module 04a) Specific Inputs ---")
        params["reduction_method"] = ask_question("Enter preferred reduction method for plotting", "umap", ["umap_c2s", "umap", "harmony", "pca"])
        params["featureplot_genes"] = ask_question("Enter comma-separated genes for FeaturePlot. Leave empty to skip.", "")
        params["dotplot_gene_groups"] = ask_for_dotplot_genes()
    else:
        params["reduction_method"], params["featureplot_genes"], params["dotplot_gene_groups"] = "umap", "", []
    params["de_gsea_plot_gene"] = ask_question("Enter a single gene for violin plot in module 4b. Leave empty to skip.", "") if "04b_DE_gsea" in params["selected_modules"] else ""
    if "04d_ucell_scores" in params["selected_modules"]:
        print("\n--- UCell Gene Scores (Module 04d) Specific Inputs ---")
        params["msigdb_category"] = ask_question(f"Enter MSigDB category for UCell (e.g., H, C2, C5).", "H").upper()
        params["ucell_plot_pathway_name"] = ask_question("Enter specific pathway name to plot for UCell. Leave empty for first.", "")
    else:
        params["msigdb_category"], params["ucell_plot_pathway_name"] = "H", ""
    if "04h_cell_chat" in params["selected_modules"]:
        print("\n--- Cell-Cell Communication (Module 04h) Specific Inputs ---")
        params["cellchat_source_groups"] = ask_question("Enter the source cell groups/clusters", "")
        params["cellchat_target_groups"] = ask_question("Enter the target cell groups/clusters", "")
        params["liana_method"] = ask_question(f"Enter LIANA method(s)", "logfc", ["natmi", "connectome", "logfc", "sca", "cellphonedb"])
    if "03_cell_type_annotation" in params["selected_modules"] and params.get("annotation_method") == "singler":
        print("\n--- SingleR Reference Configuration (Module 03) ---")
        ref_suggestion = downloaded_data.get("default_reference_data")
        if ref_suggestion: print(f"Tip: A downloaded reference is available at: {ref_suggestion}")
        local_ref_paths = ask_for_paths(f"Enter path to your local SingleR reference RDS file", is_optional=False, ensure_file=True, default_path_suggestion=ref_suggestion)
        if local_ref_paths: params["local_singler_ref_path"] = local_ref_paths[0]
        else: sys.exit("CRITICAL ERROR: SingleR selected, but no reference file provided.")
        default_label = "label.main"
        if params["local_singler_ref_path"] and "bmcite_demo.rds" in params["local_singler_ref_path"]: default_label = "celltype.l1"
        params["local_singler_ref_label_col"] = ask_question("Enter name of the metadata column with cell type labels", default_label)
    else:
        params["local_singler_ref_path"], params["local_singler_ref_label_col"] = None, None
    params["conditional_paths"] = {"query_rds_path": None, "query_species": None, "palantir_start_cell": None}
    if "04e_pseudotime" in params["selected_modules"]:
        params["conditional_paths"]["palantir_start_cell"] = ask_question("Enter start cell for Palantir", "first_cell_barcode")
    if "04f_query_projection" in params["selected_modules"]:
        query_rds_paths = ask_for_paths("Enter full path to query.RDS file", ensure_file=True)
        if query_rds_paths: params["conditional_paths"]["query_rds_path"] = query_rds_paths[0]
        else: sys.exit(f"CRITICAL ERROR: Module 04f selected, but no query.RDS path provided.")
        params["conditional_paths"]["query_species"] = ask_question("Enter species of query dataset", params["species"], ["human", "mouse"]).lower()
    params["ollama_base_url"] = OLLAMA_BASE_URL_DEFAULT
    save_params(params)
    print("\nCustom mode setup complete!")
    return True

def save_params(params: Dict[str, Any]) -> bool:
    """Save parameters to JSON file."""
    try:
        output_dir = params["output_directory"]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        params_path = normalize_path(os.path.join(output_dir, "impact_sc_params.json"))
        with open(params_path, 'w', encoding='utf-8') as f:
            json.dump(params, f, indent=4)
        print(f"\nParameters saved to: {params_path}")
        print(f"\n--- Setup Complete ---\nNext steps:")
        print(f"1. Ensure all dependencies are installed.")
        print(f"2. Activate the conda environment: conda activate impact_sc")
        print(f"3. Run the pipeline: python run_impact_sc_pipeline.py {params_path}")
        return True
    except (IOError, Exception) as e:
        print(f"Error saving parameters: {e}")
        return False

if __name__ == "__main__":
    if not main():
        print("Setup did not complete successfully.")
        sys.exit(1)
