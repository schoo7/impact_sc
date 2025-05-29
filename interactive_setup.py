#!/usr/bin/env python3

import json
import os
import sys # Import sys for sys.exit()
from typing import List, Dict, Any
import subprocess # Ensure subprocess is imported

# Configuration for Ollama (if used for more advanced interaction in future)
OLLAMA_MODEL_NAME_DEFAULT = "gemma3:12b-it-qat" # You can change the default Ollama model if needed
OLLAMA_BASE_URL_DEFAULT = "http://localhost:11434" # Base URL for Ollama API

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
            choice_map = {c.lower(): c for c in choices}
            if response.lower() in choice_map:
                return choice_map[response.lower()] 
            else:
                print(f"Invalid choice. Please select from: {', '.join(choices)}")
        elif response: 
            return response
        elif default_value is None and not response: # No default and no response, and field is mandatory implicitly
             print("This field cannot be empty. Please provide a value.")
        else: 
            return response


def ask_for_paths(prompt: str, allow_multiple: bool = False, is_optional: bool = False, optional_default_skip: str = "skip", ensure_file: bool = False, ensure_dir: bool = False, default_path_suggestion: str = None) -> List[str]:
    """Asks for file/directory paths and validates their existence and type if specified."""
    paths = []
    
    prompt_to_display = prompt
    if default_path_suggestion:
        prompt_to_display += f" (suggested: {default_path_suggestion})"

    actual_prompt_text = prompt_to_display 
    if is_optional:
        actual_prompt_text += f" (enter '{optional_default_skip}' to skip or if not applicable)"
    
    first_path = True
    while True:
        if not allow_multiple and not first_path: 
            break

        current_prompt = actual_prompt_text if first_path else "Enter additional path"
        if first_path and is_optional and not allow_multiple : 
             pass 
        elif not first_path and allow_multiple and is_optional: 
            current_prompt += f" (or type '{optional_default_skip}' to finish adding paths)"

        path_input = input(f"{current_prompt}: ").strip()

        if is_optional and path_input.lower() == optional_default_skip:
            if allow_multiple and not first_path: 
                break 
            return [] 

        if not path_input and default_path_suggestion and first_path:
            path_input = default_path_suggestion
            print(f"Using suggested path: {path_input}")

        if is_optional and not path_input and not default_path_suggestion: 
            print(f"No path entered. Assuming skip for optional path: {prompt}")
            return []

        if not path_input and not default_path_suggestion and not is_optional:
            print("This path cannot be empty. Please provide a value.")
            continue

        abs_path = os.path.abspath(path_input)
        path_exists = os.path.exists(abs_path)
        error_message = ""

        if not path_exists:
            # If ensure_file or ensure_dir is not True, we might accept a non-existing path
            # (e.g., for a Hugging Face model name that isn't a local path)
            if ensure_file or ensure_dir:
                 error_message = f"Error: Path '{path_input}' (resolved to '{abs_path}') does not exist."
            # If neither ensure_file nor ensure_dir, assume it might be a model name, so no error yet
        elif ensure_file and not os.path.isfile(abs_path):
            error_message = f"Error: Path '{path_input}' (resolved to '{abs_path}') is not a file."
        elif ensure_dir and not os.path.isdir(abs_path):
            error_message = f"Error: Path '{path_input}' (resolved to '{abs_path}') is not a directory."
        
        if not error_message:
            paths.append(abs_path if (ensure_file or ensure_dir or path_exists) else path_input) # Store original if not a validated local path
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
        "1": "01_data_processing",
        "2a": "02a_harmony_c2s_prep",
        "2b": "02b_c2s",
        "2c": "02c_load_c2s_result",
        "3": "03_cell_type_annotation",
        "4a": "04a_basic_visulization",
        "4b": "04b_DE_gsea",
        "4c": "04c_decoupler",
        "4d": "04d_ucell_scores",
        "4e": "04e_pseudotime", 
        "4f": "04f_query_projection"
    }
    print("\nAvailable IMPACT-sc Modules:")
    for key, name in all_modules.items():
        print(f"  {key}: {name}")

    selected_keys_str = ask_question("Enter the keys of modules to run, separated by commas (e.g., 1,2a,2b,2c,3,4a)")
    selected_keys = [key.strip().lower() for key in selected_keys_str.split(',')]

    selected_modules_list = []
    for key in selected_keys:
        if key in all_modules:
            selected_modules_list.append(all_modules[key])
        else:
            print(f"Warning: Module key '{key}' is not valid and will be ignored.")

    if not selected_modules_list:
        print("No valid modules selected. Defaulting to Module 1 only.")
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

def main():
    print("--- Welcome to IMPACT-sc Interactive Setup ---")
    print("Please provide the following information:")

    params: Dict[str, Any] = {}

    print("\n--- Script Locations ---")
    pipeline_base_dir = os.path.dirname(os.path.abspath(__file__))
    default_scripts_dir_suggestion = os.path.abspath(os.path.join(pipeline_base_dir, "..", "scripts.AI"))
    if not os.path.isdir(default_scripts_dir_suggestion):
        default_scripts_dir_suggestion = os.path.abspath(os.path.join(pipeline_base_dir, "scripts.AI"))

    scripts_dir_paths = ask_for_paths(
        "Enter the full path to the directory containing R and Python module scripts",
        ensure_dir=True,
        default_path_suggestion=default_scripts_dir_suggestion if os.path.isdir(default_scripts_dir_suggestion) else None
    )
    if not scripts_dir_paths: 
        print("CRITICAL ERROR: Scripts directory not provided or invalid. Exiting.")
        sys.exit(1)
    params["input_r_scripts_dir"] = scripts_dir_paths[0]
    params["input_python_scripts_dir"] = scripts_dir_paths[0]


    print("\n--- Rscript Executable Path ---")
    default_rscript_path_suggestion = None
    potential_r_paths = []
    if os.name == 'nt':
        try:
            result = subprocess.run(['where', 'Rscript.exe'], capture_output=True, text=True, check=False, timeout=5)
            if result.returncode == 0 and result.stdout.strip():
                default_rscript_path_suggestion = result.stdout.strip().splitlines()[0]
        except (FileNotFoundError, subprocess.TimeoutExpired):
             pass
        if not default_rscript_path_suggestion:
            potential_r_paths.extend([
                "C:\\Program Files\\R\\R-4.3.0\\bin\\x64\\Rscript.exe",
                "C:\\Program Files\\R\\R-4.2.3\\bin\\x64\\Rscript.exe",
                "F:\\R-4.2.3\\bin\\x64\\Rscript.exe"
            ])
    else:
        try:
            result = subprocess.run(['which', 'Rscript'], capture_output=True, text=True, check=False, timeout=5)
            if result.returncode == 0 and result.stdout.strip():
                default_rscript_path_suggestion = result.stdout.strip()
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
        if not default_rscript_path_suggestion:
            potential_r_paths.extend([
                "/usr/local/bin/Rscript",
                "/usr/bin/Rscript"
            ])
    if not default_rscript_path_suggestion:
        for r_path in potential_r_paths:
            if os.path.exists(r_path):
                default_rscript_path_suggestion = r_path
                break
    
    rscript_exe_paths = ask_for_paths(
        "Enter the full path to your Rscript executable",
        ensure_file=True,
        default_path_suggestion=default_rscript_path_suggestion
    )
    if not rscript_exe_paths: 
        print("CRITICAL ERROR: Rscript executable path not provided or invalid. Exiting.")
        sys.exit(1)
    params["rscript_executable_path"] = rscript_exe_paths[0]


    print("\n--- Input Data ---")
    params["input_data_paths"] = ask_for_paths(
        "Enter the full path to your primary input scRNA-seq data file(s) (e.g., feature-barcode matrix directory, or an .RDS file like 'ori.RDS')",
        allow_multiple=True,
        is_optional=False
    )
    params["species"] = ask_question("Enter the species ('human' or 'mouse')", "human", choices=["human", "mouse"]).lower()

    print("\n--- Output Configuration ---")
    default_output_dir_suggestion = "impact_output"
    if params["input_data_paths"]: 
        first_input_path = params["input_data_paths"][0]
        if os.path.isdir(first_input_path):
            first_input_parent_dir = first_input_path
        else: 
            first_input_parent_dir = os.path.dirname(first_input_path)
        default_output_dir_suggestion = os.path.join(first_input_parent_dir, "impact_output")

    output_dir_input = ask_question("Enter the full path for your desired output/results folder", default_output_dir_suggestion)
    params["output_directory"] = os.path.abspath(output_dir_input)


    print("\n--- Module Selection ---")
    params["selected_modules"] = select_modules()

    params["h5ad_path_for_c2s"] = None 
    params["c2s_model_path_or_name"] = None # Initialize
    if "02b_c2s" in params["selected_modules"]:
        print("\n--- Cell2Sentence (Module 02b) Specific Input ---")
        h5ad_c2s_paths = ask_for_paths(
            "Enter the full path to the H5AD file for Cell2Sentence input (e.g., F:\\R_PROJECT\\impact/02_module2_for_c2s.h5ad)",
            is_optional=False, 
            ensure_file=True
        )
        if h5ad_c2s_paths:
            params["h5ad_path_for_c2s"] = h5ad_c2s_paths[0]
            print(f"Using H5AD file for Cell2Sentence: {params['h5ad_path_for_c2s']}")
        else:
            print("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no valid H5AD input path was provided. The pipeline will likely fail.")
            sys.exit(1)

        # Ask for Cell2Sentence model path or name
        default_c2s_model_name = "james-y-u/C2S-Pythia-410m-cell-type-prediction"
        c2s_model_input_paths = ask_for_paths(
            f"Enter the Cell2Sentence model path (if local, e.g., F:\\path\\to\\model_folder) or Hugging Face model name",
            is_optional=False, # Assuming model is required for this module
            ensure_dir=False, # Set to True if you want to enforce it's a local DIR, False allows HF names
            default_path_suggestion=default_c2s_model_name
        )
        if c2s_model_input_paths:
            # If it's a local path, ask_for_paths would have made it absolute if it exists.
            # If it doesn't exist locally, it's treated as a name (e.g., from Hugging Face).
            # We need to check if the user-provided path is an existing directory.
            # If it is, we use the absolute path. Otherwise, we use the input as is (could be HF name).
            potential_local_path = os.path.abspath(c2s_model_input_paths[0])
            if os.path.isdir(potential_local_path):
                params["c2s_model_path_or_name"] = potential_local_path
                print(f"Using local Cell2Sentence model from directory: {params['c2s_model_path_or_name']}")
            else:
                params["c2s_model_path_or_name"] = c2s_model_input_paths[0] # Use as is (could be HF name)
                print(f"Using Cell2Sentence model (name or non-validated path): {params['c2s_model_path_or_name']}")
        else:
            print("CRITICAL ERROR IN SETUP: Module 02b_c2s selected, but no Cell2Sentence model path/name was provided.")
            sys.exit(1)


    params["collectri_csv_path"] = None
    params["progeny_csv_path"] = None

    if "04c_decoupler" in params["selected_modules"]:
        print("\n--- DecoupleR Network CSV Paths (Module 04c) ---")
        collectri_prompt = f"Enter the full path to the CollecTRI CSV file for {params['species']} (e.g., collectri_{params['species']}.csv)"
        collectri_paths = ask_for_paths(collectri_prompt, allow_multiple=False, is_optional=False, ensure_file=True)
        if collectri_paths:
            params["collectri_csv_path"] = collectri_paths[0]
            print(f"Using CollecTRI CSV: {params['collectri_csv_path']}")
        else:
            print(f"CRITICAL ERROR IN SETUP: Module 04c (DecoupleR) selected, but no CollecTRI CSV path was obtained. The R script will fail.")
            sys.exit(1)

        if params["species"] == "human":
            progeny_prompt = "Enter the full path to the PROGENy CSV file for human (e.g., progeny_human.csv)"
            progeny_paths = ask_for_paths(progeny_prompt, allow_multiple=False, is_optional=True, optional_default_skip="skip", ensure_file=True)
            if progeny_paths:
                params["progeny_csv_path"] = progeny_paths[0]
                print(f"Using PROGENy CSV: {params['progeny_csv_path']}")
            else:
                print("No PROGENy CSV path provided for human (or skipped). R script will skip PROGENy analysis if path is not set/valid.")
        else:
            print(f"PROGENy pathway analysis is typically for human. Skipping PROGENy CSV path question for {params['species']}.")
            params["progeny_csv_path"] = None
    else:
        params["collectri_csv_path"] = None
        params["progeny_csv_path"] = None


    if "03_cell_type_annotation" in params["selected_modules"]:
        print("\n--- Final Cell Type Source (Module 03) ---")
        params["final_cell_type_source"] = ask_question(
            "Which annotation source should be used for the final 'cell_type' column? (auto, seurat, c2s, singler)",
            "auto",
            choices=["auto", "seurat", "c2s", "singler"]
        ).lower()
    else:
        params["final_cell_type_source"] = "auto"

    if "04a_basic_visulization" in params["selected_modules"]:
        print("\n--- Basic Visualization (Module 04a) Specific Inputs ---")
        params["featureplot_genes"] = ask_question(
            "Enter comma-separated genes for FeaturePlot (e.g., CD3D,CD14,MS4A1). Leave empty to skip.",
            ""
        )
        params["dotplot_gene_groups"] = ask_for_dotplot_genes()
    else:
        params["featureplot_genes"] = ""
        params["dotplot_gene_groups"] = []

    if "04b_DE_gsea" in params["selected_modules"]:
        print("\n--- Differential Expression & GSEA (Module 04b) Specific Inputs ---")
        params["de_gsea_plot_gene"] = ask_question(
            "Enter a single gene name for the expression violin plot in module 4b (e.g., CD3D). Leave empty to skip.",
            ""
        )
    else:
        params["de_gsea_plot_gene"] = ""

    if "04d_ucell_scores" in params["selected_modules"]:
        print("\n--- UCell Gene Scores (Module 04d) Specific Inputs ---")
        msigdb_categories_common = ["H", "C2", "C5", "C7", "C8"]
        params["msigdb_category"] = ask_question(
            f"Enter the MSigDB category for UCell (e.g., H, C2, C5). Common choices: {', '.join(msigdb_categories_common)}. See MSigDB website for all.",
            "H"
        ).upper()
        params["ucell_plot_pathway_name"] = ask_question(
            "Enter the specific pathway/gene set name to plot for UCell (e.g., HALLMARK_APOPTOSIS). Leave empty to plot the first available.",
            ""
        )
    else:
        params["msigdb_category"] = "H"
        params["ucell_plot_pathway_name"] = ""


    params["local_singler_ref_path"] = None
    module_3_is_selected = "03_cell_type_annotation" in params["selected_modules"]

    if params["species"] in ["human", "mouse"]:
        prompt_text = f"Enter the full path to a local SingleR reference RDS file for {params['species']} (e.g., F:/R_PROJECT/impact/ref/{params['species']}ref.RDS)."
        is_ref_optional_for_module3 = not module_3_is_selected

        if module_3_is_selected:
            print(f"\n--- Required Local SingleR Reference for {params['species'].capitalize()} (Module 3 is selected) ---")
        else:
            print(f"\n--- Optional Local SingleR Reference for {params['species'].capitalize()} ---")
            prompt_text += " This is optional as Module 3 is not selected (but will be saved if provided)."

        local_ref_paths = ask_for_paths(
            prompt_text,
            allow_multiple=False,
            is_optional=is_ref_optional_for_module3,
            optional_default_skip="skip",
            ensure_file=True
        )

        if local_ref_paths:
            params["local_singler_ref_path"] = local_ref_paths[0]
            print(f"Using local SingleR reference: {params['local_singler_ref_path']}")
        elif module_3_is_selected:
            print(f"CRITICAL ERROR IN SETUP: Module 3 selected, but no local SingleR reference path was obtained. The R script will fail. Please ensure a valid path is entered.")
            sys.exit(1)
        else:
            print("No local SingleR reference provided (optional and skipped, or Module 3 not selected).")


    params["conditional_paths"] = {
        "query_rds_path": None,
        "query_species": None,
        "palantir_start_cell": None
    }


    if "04f_pseudotime" in params["selected_modules"]:
        print("\n--- Pseudotime Analysis (Module 04f) Specific Input ---")
        params["conditional_paths"]["palantir_start_cell"] = ask_question(
            "Enter the barcode/name of a known progenitor/start cell for Palantir pseudotime analysis (e.g., AAACCCAAGTCGAAGG-1)",
            "first_cell_barcode" 
        )

    if "04g_query_projection" in params["selected_modules"]:
        print("\n--- Query Dataset Projection (Module 04g) Specific Input ---")
        query_rds_paths = ask_for_paths(
            "Enter the full path to the query.RDS file (e.g., F:/R_PROJECT/impact/ref/query.RDS)",
            allow_multiple=False,
            is_optional=False,
            ensure_file=True
        )
        if query_rds_paths:
            params["conditional_paths"]["query_rds_path"] = query_rds_paths[0]
            print(f"Query RDS path set to: {params['conditional_paths']['query_rds_path']}")
        else:
            print(f"CRITICAL ERROR IN SETUP: Module 04g selected, but no query.RDS path was provided. R script will fail.")
            sys.exit(1)

        params["conditional_paths"]["query_species"] = ask_question(
             "Enter the species of the query dataset ('human' or 'mouse')",
             params["species"],
             choices=["human", "mouse"]
        ).lower()


    params["ollama_model_name"] = OLLAMA_MODEL_NAME_DEFAULT
    params["ollama_base_url"] = OLLAMA_BASE_URL_DEFAULT

    try:
        if not os.path.exists(params["output_directory"]):
            os.makedirs(params["output_directory"], exist_ok=True)
            print(f"Created output directory: {params['output_directory']}")

        params_path = os.path.join(params["output_directory"], "impact_sc_params.json")
        with open(params_path, 'w', encoding='utf-8') as f:
            json.dump(params, f, indent=4)
        print(f"\nParameters saved to: {params_path}")

        print(f"\n--- Setup Complete ---")
        print(f"Next steps:")
        print(f"1. Ensure all R and Python dependencies have been correctly installed (run install_dependencies.sh if needed).")
        print(f"2. Activate the 'impact_sc' conda environment: conda activate impact_sc")
        print(f"3. Run the pipeline using: python run_impact_sc_pipeline.py {params_path}")
        print(f"   IMPORTANT: The 'run_impact_sc_pipeline.py' script MUST read all relevant parameters from '{params_path}' and set them as appropriate environment variables.")
        print(f"   This includes 'h5ad_path_for_c2s' and 'c2s_model_path_or_name' for Module 02b, 'collectri_csv_path' and 'progeny_csv_path' for Module 04c, and the RDS path for Module 04g.")

    except IOError as e:
        print(f"Error: Could not write parameters file to {params_path}. {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
