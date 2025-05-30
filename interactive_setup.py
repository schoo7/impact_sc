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
        "4a": "04a_basic_visualization",
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

def check_downloaded_data() -> Dict[str, str]:
    """Check for downloaded data and return paths if found."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_paths = {
        "demo_data": None,
        "c2s_model": None,
        "reference_data": None
    }
    
    # Check for demo data
    demo_path = os.path.join(current_dir, "data", "demo", "filtered_gene_bc_matrices")
    if os.path.exists(demo_path):
        # Check if this is a 10x format directory structure
        hg19_path = os.path.join(demo_path, "hg19")
        if os.path.exists(hg19_path):
            # Check for required 10x files
            matrix_file = os.path.join(hg19_path, "matrix.mtx")
            genes_file = os.path.join(hg19_path, "genes.tsv")
            barcodes_file = os.path.join(hg19_path, "barcodes.tsv")
            
            if all(os.path.exists(f) for f in [matrix_file, genes_file, barcodes_file]):
                data_paths["demo_data"] = hg19_path  # Point to hg19 subdirectory
                print(f"✅ Found demo data: {hg19_path}")
            else:
                print(f"⚠️ Demo directory found but missing 10x files: {demo_path}")
        else:
            # Fallback to original path if no hg19 subdirectory
            data_paths["demo_data"] = demo_path
            print(f"✅ Found demo data: {demo_path}")
    else:
        print(f"⚠️ Demo data not found at: {demo_path}")
    
    # Check for Cell2Sentence model cache
    models_dir = os.path.join(current_dir, "data", "models")
    if os.path.exists(models_dir):
        # Look for cached model files
        import glob
        model_files = glob.glob(os.path.join(models_dir, "**/*C2S*"), recursive=True) + \
                     glob.glob(os.path.join(models_dir, "**/*Pythia*"), recursive=True)
        if model_files:
            data_paths["c2s_model"] = "vandijklab/C2S-Pythia-410m-cell-type-prediction"  # Use HF name with cache
            print(f"✅ Found Cell2Sentence model cache: {models_dir}")
        else:
            print(f"⚠️ Models directory exists but no cached files found: {models_dir}")
    
    # Check for reference data
    ref_path = os.path.join(current_dir, "data", "reference", "HumanPrimaryCellAtlasData.rds")
    if os.path.exists(ref_path):
        data_paths["reference_data"] = ref_path
        print(f"✅ Found reference data: {ref_path}")
    
    return data_paths

def ask_demo_mode() -> bool:
    """Ask if user wants to use demo mode."""
    print("\n" + "="*60)
    print("🚀 IMPACT-sc Setup Mode Selection")
    print("="*60)
    
    downloaded_data = check_downloaded_data()
    has_demo_data = downloaded_data["demo_data"] is not None
    
    if has_demo_data:
        print("✅ Demo data detected! You can run in demo mode.")
        print("📂 Demo data includes:")
        print("   • PBMC3k single-cell dataset (10x Genomics)")
        print("   • Pre-configured parameters for quick testing")
        print("")
        
        mode_choice = ask_question(
            "Choose setup mode:\n  [demo] - Use demo data and auto-configure\n  [custom] - Configure your own data",
            "demo",
            choices=["demo", "custom"]
        ).lower()
        
        return mode_choice == "demo"
    else:
        print("⚠️ No demo data found. Run './download_data.sh' first to enable demo mode.")
        print("Proceeding with custom setup...")
        return False

def main():
    print("--- Welcome to IMPACT-sc Interactive Setup ---")
    
    # Check for downloaded data
    downloaded_data = check_downloaded_data()
    use_demo_mode = ask_demo_mode()
    
    if use_demo_mode:
        return setup_demo_mode(downloaded_data)
    else:
        return setup_custom_mode(downloaded_data)


def setup_demo_mode(downloaded_data: Dict[str, str]):
    """Setup demo mode using downloaded data and pre-configured parameters."""
    print("\n🎯 Setting up DEMO MODE...")
    print("Using pre-configured parameters for PBMC3k dataset.")
    
    # Ask if user wants to include C2S (which can be slow)
    include_c2s = ask_question(
        "Include Cell2Sentence (C2S) analysis? This provides AI-powered cell type prediction but takes longer",
        "no",
        choices=["yes", "no"]
    ).lower() == "yes"

    params: Dict[str, Any] = {}
    pipeline_base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Basic configuration
    params["input_r_scripts_dir"] = os.path.abspath(os.path.join(pipeline_base_dir, "scripts_AI"))
    params["input_python_scripts_dir"] = os.path.abspath(os.path.join(pipeline_base_dir, "scripts_AI"))
    params["rscript_executable_path"] = auto_detect_rscript()
    params["species"] = "human"
    params["output_directory"] = os.path.abspath("demo_output")
    
    # Demo data configuration
    if downloaded_data["demo_data"]:
        params["input_data_paths"] = [downloaded_data["demo_data"]]
        print(f"📂 Using demo data: {downloaded_data['demo_data']}")
    else:
        print("❌ Demo data not found!")
        return False
    
    # Module selection based on C2S choice
    if include_c2s:
        print("🤖 Including Cell2Sentence analysis (slower but more comprehensive)")
        params["selected_modules"] = [
            "01_data_processing",
            "02a_harmony_c2s_prep", 
            "02b_c2s",
            "02c_load_c2s_result",
            "03_cell_type_annotation",
            "04a_basic_visualization"
        ]
        
        # Cell2Sentence configuration
        params["h5ad_path_for_c2s"] = os.path.join(params["output_directory"], "02_module2_for_c2s.h5ad")
        if downloaded_data["c2s_model"]:
            params["c2s_model_path_or_name"] = downloaded_data["c2s_model"]
            print(f"🤖 Using cached Cell2Sentence model")
        else:
            params["c2s_model_path_or_name"] = "vandijklab/C2S-Pythia-410m-cell-type-prediction"
            print("⚠️ No cached model found. Will download on first use.")
    else:
        print("⚡ Skipping Cell2Sentence for faster demo (using traditional annotation only)")
        params["selected_modules"] = [
            "01_data_processing",
            "03_cell_type_annotation",
            "04a_basic_visualization"
        ]
        
        # No C2S configuration needed
        params["h5ad_path_for_c2s"] = None
        params["c2s_model_path_or_name"] = None
    
    # Reference data
    if downloaded_data["reference_data"]:
        params["local_singler_ref_path"] = downloaded_data["reference_data"]
        print(f"📚 Using downloaded reference data")
    else:
        params["local_singler_ref_path"] = None
        print("⚠️ No reference data found. Will use online celldex.")
    
    # Default visualization genes for PBMC
    params["featureplot_genes"] = "CD3D,CD14,MS4A1,FCGR3A,LYZ,PPBP"
    params["dotplot_gene_groups"] = [
        {"name": "T_cell_markers", "genes": ["CD3D", "CD3E", "CD8A", "CD4"]},
        {"name": "B_cell_markers", "genes": ["MS4A1", "CD79A", "CD79B"]},
        {"name": "Myeloid_markers", "genes": ["CD14", "LYZ", "FCGR3A", "CST3"]}
    ]
    
    # Other defaults
    params["final_cell_type_source"] = "auto"
    params["de_gsea_plot_gene"] = "CD3D"
    params["collectri_csv_path"] = None
    params["progeny_csv_path"] = None
    params["msigdb_category"] = "H"
    params["ucell_plot_pathway_name"] = ""
    params["conditional_paths"] = {
        "query_rds_path": None,
        "query_species": None,
        "palantir_start_cell": None
    }
    params["ollama_model_name"] = "gemma3:12b-it-qat"
    params["ollama_base_url"] = "http://localhost:11434"
    
    # Save parameters
    save_params(params)
    print("\n🎉 Demo mode setup complete!")
    print("Run: conda activate impact_sc && python run_impact_sc_pipeline.py demo_output/impact_sc_params.json")
    return True


def auto_detect_rscript() -> str:
    """Auto-detect Rscript executable."""
    import subprocess
    import platform
    
    try:
        if platform.system() == "Windows":
            result = subprocess.run(['where', 'Rscript.exe'], capture_output=True, text=True, timeout=5)
        else:
            result = subprocess.run(['which', 'Rscript'], capture_output=True, text=True, timeout=5)
        
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip().split('\n')[0]
    except:
        pass
    
    # Fallback paths
    fallback_paths = []
    if platform.system() == "Darwin":  # macOS
        fallback_paths = [
            "/opt/homebrew/bin/Rscript",
            "/usr/local/bin/Rscript", 
            "/Library/Frameworks/R.framework/Resources/bin/Rscript"
        ]
    elif platform.system() == "Windows":
        fallback_paths = [
            "C:\\Program Files\\R\\R-4.3.0\\bin\\x64\\Rscript.exe",
            "C:\\Program Files\\R\\R-4.2.3\\bin\\x64\\Rscript.exe"
        ]
    else:  # Linux
        fallback_paths = ["/usr/local/bin/Rscript", "/usr/bin/Rscript"]
    
    for path in fallback_paths:
        if os.path.exists(path):
            return path
    
    return "Rscript"  # Default fallback


def setup_custom_mode(downloaded_data: Dict[str, str]):
    """Setup custom mode with user inputs, using downloaded data as defaults."""
    print("\n⚙️ Setting up CUSTOM MODE...")
    print("Please provide the following information:")

    params: Dict[str, Any] = {}

    print("\n--- Script Locations ---")
    pipeline_base_dir = os.path.dirname(os.path.abspath(__file__))
    default_scripts_dir_suggestion = os.path.abspath(os.path.join(pipeline_base_dir, "..", "scripts_AI"))
    if not os.path.isdir(default_scripts_dir_suggestion):
        default_scripts_dir_suggestion = os.path.abspath(os.path.join(pipeline_base_dir, "scripts_AI"))

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
    default_rscript_path_suggestion = auto_detect_rscript()
    
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
    # Use demo data as default if available
    demo_suggestion = downloaded_data["demo_data"] if downloaded_data["demo_data"] else None
    if demo_suggestion:
        print(f"💡 Tip: Demo data available at: {demo_suggestion}")
    
    params["input_data_paths"] = ask_for_paths(
        "Enter the full path to your primary input scRNA-seq data file(s) (e.g., feature-barcode matrix directory, or an .RDS file like 'ori.RDS')",
        allow_multiple=True,
        is_optional=False,
        default_path_suggestion=demo_suggestion
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

        # Use downloaded model as default
        default_c2s_model = downloaded_data["c2s_model"] if downloaded_data["c2s_model"] else "vandijklab/C2S-Pythia-410m-cell-type-prediction"
        if downloaded_data["c2s_model"]:
            print(f"💡 Tip: Cached model available: {default_c2s_model}")
        
        c2s_model_input_paths = ask_for_paths(
            f"Enter the Cell2Sentence model path (if local, e.g., F:\\path\\to\\model_folder) or Hugging Face model name",
            is_optional=False,
            ensure_dir=False,
            default_path_suggestion=default_c2s_model
        )
        if c2s_model_input_paths:
            potential_local_path = os.path.abspath(c2s_model_input_paths[0])
            if os.path.isdir(potential_local_path):
                params["c2s_model_path_or_name"] = potential_local_path
                print(f"Using local Cell2Sentence model from directory: {params['c2s_model_path_or_name']}")
            else:
                params["c2s_model_path_or_name"] = c2s_model_input_paths[0]
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

    if "04a_basic_visualization" in params["selected_modules"]:
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

        # Use downloaded reference data as default
        ref_suggestion = downloaded_data["reference_data"] if downloaded_data["reference_data"] else None
        if ref_suggestion:
            print(f"💡 Tip: Downloaded reference data available: {ref_suggestion}")
            prompt_text = f"Enter the full path to a local SingleR reference RDS file for {params['species']}"

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
            ensure_file=True,
            default_path_suggestion=ref_suggestion
        )

        if local_ref_paths:
            params["local_singler_ref_path"] = local_ref_paths[0]
            print(f"Using local SingleR reference: {params['local_singler_ref_path']}")
        elif module_3_is_selected:
            if downloaded_data["reference_data"]:
                print("⚠️ Using downloaded reference data since no custom path provided.")
                params["local_singler_ref_path"] = downloaded_data["reference_data"]
            else:
                print(f"CRITICAL ERROR IN SETUP: Module 3 selected, but no local SingleR reference path was obtained. The R script will fail. Please ensure a valid path is entered.")
                sys.exit(1)
        else:
            print("No local SingleR reference provided (optional and skipped, or Module 3 not selected).")


    params["conditional_paths"] = {
        "query_rds_path": None,
        "query_species": None,
        "palantir_start_cell": None
    }


    if "04e_pseudotime" in params["selected_modules"]:
        print("\n--- Pseudotime Analysis (Module 04e) Specific Input ---")
        params["conditional_paths"]["palantir_start_cell"] = ask_question(
            "Enter the barcode/name of a known progenitor/start cell for Palantir pseudotime analysis (e.g., AAACCCAAGTCGAAGG-1)",
            "first_cell_barcode" 
        )

    if "04f_query_projection" in params["selected_modules"]:
        print("\n--- Query Dataset Projection (Module 04f) Specific Input ---")
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
            print(f"CRITICAL ERROR IN SETUP: Module 04f selected, but no query.RDS path was provided. R script will fail.")
            sys.exit(1)

        params["conditional_paths"]["query_species"] = ask_question(
             "Enter the species of the query dataset ('human' or 'mouse')",
             params["species"],
             choices=["human", "mouse"]
        ).lower()


    params["ollama_model_name"] = "gemma3:12b-it-qat"
    params["ollama_base_url"] = "http://localhost:11434"

    # Save parameters
    save_params(params)
    print("\n✅ Custom mode setup complete!")
    return True


def save_params(params: Dict[str, Any]) -> bool:
    """Save parameters to JSON file."""
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
        print(f"1. Ensure all R and Python dependencies have been correctly installed.")
        print(f"2. Activate the 'impact_sc' conda environment: conda activate impact_sc")
        print(f"3. Run the pipeline using: python run_impact_sc_pipeline.py {params_path}")
        
        return True
    except IOError as e:
        print(f"Error: Could not write parameters file. {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

if __name__ == "__main__":
    if main():
        # If setup is successful, exit
        sys.exit(0)
    else:
        # If setup fails, exit with a non-zero code
        sys.exit(1)
