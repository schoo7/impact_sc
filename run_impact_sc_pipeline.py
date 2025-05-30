#!/usr/bin/env python3

import json
import subprocess
import os
import sys
import shutil # For shutil.which

def run_script(script_path: str, script_type: str, params: dict, module_name: str) -> bool:
    """
    Runs an R or Python script with necessary environment variables and parameters.
    """
    output_dir = params["output_directory"]
    base_script_name = os.path.basename(script_path)
    log_file_path = os.path.join(output_dir, f"{module_name}_log.txt")

    print(f"--- Running {module_name} ({base_script_name}) ---")
    print(f"Output log: {log_file_path}")

    env = os.environ.copy()
    env["IMPACT_SC_SPECIES"] = params.get("species", "human")
    env["IMPACT_SC_OUTPUT_DIR"] = output_dir 

    if params.get("input_data_paths") and len(params["input_data_paths"]) > 0:
        env["IMPACT_SC_INPUT_DATA_PATH"] = params["input_data_paths"][0]
    else:
        print(f"Warning: 'input_data_paths' is empty or not defined in params.json. Module {module_name} might fail if it requires IMPACT_SC_INPUT_DATA_PATH.")

    if params.get("input_data_paths"):
        env["IMPACT_SC_INPUT_DATA_PATHS"] = ";".join(params["input_data_paths"])

    # --- Environment variable settings for specific modules ---
    if module_name == "02b_c2s":
        h5ad_input_for_c2s = params.get("h5ad_path_for_c2s")
        c2s_model_path_or_name = params.get("c2s_model_path_or_name") # Get the model path/name
        
        if not (h5ad_input_for_c2s and os.path.exists(h5ad_input_for_c2s)):
            error_msg = f"ERROR: H5AD file path for Cell2Sentence (module 02b_c2s) is not set or invalid in params.json: '{h5ad_input_for_c2s}'."
            print(error_msg)
            with open(log_file_path, 'a', encoding='utf-8', errors='replace') as lf: lf.write(f"\nOrchestrator Error: {error_msg}\n")
            return False
        env["H5AD_FILE_PATH"] = h5ad_input_for_c2s
        
        if not c2s_model_path_or_name:
            error_msg = f"ERROR: Cell2Sentence model path/name (module 02b_c2s) is not set in params.json."
            print(error_msg)
            with open(log_file_path, 'a', encoding='utf-8', errors='replace') as lf: lf.write(f"\nOrchestrator Error: {error_msg}\n")
            return False
        env["C2S_MODEL_PATH_OR_NAME"] = c2s_model_path_or_name # Set the new environment variable

        c2s_module_specific_output_dir = os.path.join(output_dir, "cell2sentence_module_outputs")
        env["C2S_OUTPUT_DIR"] = c2s_module_specific_output_dir
        env["C2S_EMBEDDINGS_CSV"] = os.path.join(c2s_module_specific_output_dir, "c2s_embeddings.csv")
        env["C2S_PREDICTED_CSV"] = os.path.join(c2s_module_specific_output_dir, "c2s_predictions.csv")

        print(f"Setting environment variables for module {module_name}:")
        print(f"  H5AD_FILE_PATH = {env['H5AD_FILE_PATH']}")
        print(f"  C2S_MODEL_PATH_OR_NAME = {env['C2S_MODEL_PATH_OR_NAME']}") # Print the model path/name
        print(f"  C2S_OUTPUT_DIR = {env['C2S_OUTPUT_DIR']}")
        print(f"  C2S_EMBEDDINGS_CSV = {env['C2S_EMBEDDINGS_CSV']}")
        print(f"  C2S_PREDICTED_CSV = {env['C2S_PREDICTED_CSV']}")

    elif module_name == "03_cell_type_annotation":
        local_ref_path = params.get("local_singler_ref_path")
        if local_ref_path and isinstance(local_ref_path, str) and local_ref_path.strip():
            env["IMPACT_SC_LOCAL_SINGLER_REF_PATH"] = local_ref_path
            print(f"Setting IMPACT_SC_LOCAL_SINGLER_REF_PATH for Module 3 to: {local_ref_path}")
        else:
            env["IMPACT_SC_LOCAL_SINGLER_REF_PATH"] = ""
            print(f"Warning: 'local_singler_ref_path' not found or empty for Module 3. R script might fail if it requires this.")
            
        final_cell_type_source = params.get("final_cell_type_source")
        if final_cell_type_source:
            env["IMPACT_SC_FINAL_CELL_TYPE_SOURCE"] = final_cell_type_source
            print(f"Setting IMPACT_SC_FINAL_CELL_TYPE_SOURCE for Module 3 to: {final_cell_type_source}")
        else:
            env["IMPACT_SC_FINAL_CELL_TYPE_SOURCE"] = "auto"
            print(f"Warning: 'final_cell_type_source' not found in params. Module 3 will use its default ('auto').")

    elif module_name == "04a_basic_visualization":
        featureplot_genes = params.get("featureplot_genes", "")
        env["IMPACT_SC_FEATUREPLOT_GENES"] = featureplot_genes
        print(f"Setting IMPACT_SC_FEATUREPLOT_GENES for Module 4a to: '{featureplot_genes if featureplot_genes else 'empty (skip)'}'")

        dotplot_gene_groups = params.get("dotplot_gene_groups", [])
        try:
            dotplot_json_str = json.dumps(dotplot_gene_groups if dotplot_gene_groups else [])
            env["IMPACT_SC_DOTPLOT_GENES_JSON"] = dotplot_json_str
            print(f"Setting IMPACT_SC_DOTPLOT_GENES_JSON for Module 4a to: {dotplot_json_str}")
        except TypeError as e:
            print(f"Error: Could not serialize dotplot_gene_groups to JSON: {e}. Defaulting to empty JSON array.")
            env["IMPACT_SC_DOTPLOT_GENES_JSON"] = "[]"
        
    elif module_name == "04b_DE_gsea":
        de_gsea_plot_gene = params.get("de_gsea_plot_gene", "")
        env["IMPACT_SC_DE_GENE"] = de_gsea_plot_gene
        print(f"Setting IMPACT_SC_DE_GENE for Module 4b to: '{de_gsea_plot_gene if de_gsea_plot_gene else 'empty (skip)'}'")

    elif module_name == "04c_decoupler":
        collectri_csv = params.get("collectri_csv_path", "")
        env["IMPACT_SC_COLLECTRI_CSV_PATH"] = collectri_csv if collectri_csv else ""
        if collectri_csv:
            print(f"Setting IMPACT_SC_COLLECTRI_CSV_PATH for Module 4c to: {collectri_csv}")
        else:
            print("Warning: 'collectri_csv_path' not found or empty in params for Module 4c. R script will skip TF analysis.")

        progeny_csv = params.get("progeny_csv_path", "")
        env["IMPACT_SC_PROGENY_CSV_PATH"] = progeny_csv if progeny_csv else ""
        if progeny_csv:
            print(f"Setting IMPACT_SC_PROGENY_CSV_PATH for Module 4c to: {progeny_csv}")
        else:
            print("Info: 'progeny_csv_path' not found or empty in params for Module 4c. R script will skip PROGENy analysis if applicable.")

    elif module_name == "04d_ucell_scores":
        msigdb_category = params.get("msigdb_category", "H")
        env["IMPACT_SC_MSIGDB_CATEGORY"] = msigdb_category
        print(f"Setting IMPACT_SC_MSIGDB_CATEGORY for Module 4d to: '{msigdb_category}'")

        ucell_plot_pathway_name = params.get("ucell_plot_pathway_name", "")
        env["IMPACT_SC_UCELL_PLOT_PATHWAY_NAME"] = ucell_plot_pathway_name
        print(f"Setting IMPACT_SC_UCELL_PLOT_PATHWAY_NAME for Module 4d to: '{ucell_plot_pathway_name if ucell_plot_pathway_name else 'empty (plot first)'}'")

    elif module_name == "04e_pseudotime": 
        palantir_start_cell = params.get("conditional_paths", {}).get("palantir_start_cell")
        env["IMPACT_SC_PALANTIR_START_CELL"] = palantir_start_cell if palantir_start_cell else ""
        if palantir_start_cell: print(f"Setting IMPACT_SC_PALANTIR_START_CELL: {palantir_start_cell}")
        else: print("Warning: palantir_start_cell not set for Pseudotime.")

    # elif module_name == "04f_query_projection":  # TEMPORARILY DISABLED
    #     query_rds_path = params.get("conditional_paths", {}).get("query_rds_path")
    #     if query_rds_path and os.path.exists(query_rds_path):
    #         env["IMPACT_SC_QUERY_RDS_PATH"] = query_rds_path
    #         print(f"Setting IMPACT_SC_QUERY_RDS_PATH for Module 04f to: {query_rds_path}")
    #     else:
    #         env["IMPACT_SC_QUERY_RDS_PATH"] = ""
    #         print(f"Warning: query_rds_path ('{query_rds_path}') for Module 04f not set or invalid. R script might skip or fail.")
    #     
    #     query_species = params.get("conditional_paths", {}).get("query_species")
    #     env["IMPACT_SC_QUERY_SPECIES"] = query_species if query_species else params.get("species", "human")
    #     if query_species: print(f"Setting IMPACT_SC_QUERY_SPECIES: {query_species}")
    #     else: print(f"Warning: query_species not set for Query Projection. Defaulting to main dataset species: {env['IMPACT_SC_QUERY_SPECIES']}")

    # --- Script execution logic ---
    stdout_data = None
    stderr_data = None
    process = None 

    try:
        with open(log_file_path, 'w', encoding='utf-8', errors='replace') as log_file:
            cmd = []
            cwd_to_use = os.getcwd() 

            if script_type == "R":
                rscript_exe_from_params = params.get("rscript_executable_path")
                env_name_for_r = None 

                if rscript_exe_from_params:
                    norm_path = os.path.normpath(rscript_exe_from_params)
                    path_parts = norm_path.split(os.sep)
                    if "envs" in path_parts:
                        try:
                            env_idx = path_parts.index("envs")
                            if env_idx + 1 < len(path_parts):
                                potential_env_name = path_parts[env_idx + 1]
                                if (env_idx + 2 < len(path_parts) and path_parts[env_idx + 2].lower() in ("scripts", "bin", "r", "lib")):
                                     env_name_for_r = potential_env_name
                            if env_name_for_r:
                                print(f"Derived Conda env name for Rscript: {env_name_for_r} from path: {rscript_exe_from_params}")
                            else:
                                print(f"Warning: Could not reliably derive Conda env name for Rscript from path: {rscript_exe_from_params}")
                        except ValueError:
                            print(f"Warning: 'envs' directory not found in Rscript path, cannot derive Conda env name: {rscript_exe_from_params}")
                    else:
                        print(f"Info: Rscript path '{rscript_exe_from_params}' does not appear to be a standard Conda env path. Will attempt direct execution.")
                
                conda_exe_path = shutil.which("conda")
                if env_name_for_r and conda_exe_path:
                    cmd = [conda_exe_path, "run", "-n", env_name_for_r, "Rscript", script_path]
                    log_file.write(f"Executing R script using Conda: {' '.join(cmd)}\n")
                elif rscript_exe_from_params and os.path.isfile(rscript_exe_from_params):
                    cmd = [rscript_exe_from_params, script_path]
                    log_file.write(f"Executing R script directly using path from params: {' '.join(cmd)}\n")
                elif shutil.which("Rscript"):
                    cmd = ["Rscript", script_path]
                    log_file.write(f"Warning: 'rscript_executable_path' not used or invalid, and no Conda env derived. Executing Rscript using system PATH: {' '.join(cmd)}\n")
                else:
                    error_msg = f"FATAL ERROR: Rscript executable not found via provided path, Conda, or system PATH. Cannot run R script: {script_path}"
                    print(error_msg)
                    log_file.write(f"{error_msg}\n")
                    return False
                cwd_to_use = params.get("input_r_scripts_dir", os.getcwd())

            elif script_type == "python":
                target_conda_env_python = "impact_sc" 
                conda_exe = shutil.which("conda")

                if conda_exe:
                    cmd = [conda_exe, "run", "-n", target_conda_env_python, "python", script_path]
                    log_file.write(f"Attempting to execute Python script using 'conda run': {' '.join(cmd)}\n")
                else:
                    cmd = [sys.executable, script_path] 
                    log_file.write(f"Warning: Conda executable not found. Executing Python script using sys.executable: {sys.executable}\n")
                    log_file.write(f"Command: {' '.join(cmd)}\n")
                    log_file.write(f"Ensure that the '{target_conda_env_python}' environment is active, or that sys.executable points to the correct Python interpreter with all necessary packages.\n")
                cwd_to_use = params.get("input_python_scripts_dir", os.getcwd())
            
            else:
                error_msg = f"Error: Unknown script type '{script_type}' for {base_script_name}"
                print(error_msg)
                log_file.write(f"{error_msg}\n")
                return False

            log_file.write(f"Working directory: {cwd_to_use}\n\n")
            log_file.flush() 

            # Special handling for C2S script to show real-time progress
            if module_name == "02b_c2s":
                print(f"🚀 Running {module_name} with live progress monitoring...")
                # For C2S, show real-time output while also logging
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                                         env=env, cwd=cwd_to_use, text=True, encoding='utf-8', 
                                         errors='replace', bufsize=1, universal_newlines=True)
                
                output_lines = []
                for line in iter(process.stdout.readline, ''):
                    line = line.rstrip()
                    if line:
                        print(line)  # Show real-time progress
                        output_lines.append(line)
                        log_file.write(line + '\n')
                        log_file.flush()
                
                process.wait()
                stdout_data = '\n'.join(output_lines)
                stderr_data = ""
            else:
                # Normal execution for other scripts
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                         env=env, cwd=cwd_to_use, text=True, encoding='utf-8', errors='replace')
                stdout_data, stderr_data = process.communicate()

            if stdout_data and module_name != "02b_c2s":  # Already logged for C2S above
                log_file.write("\n--- stdout ---\n")
                log_file.write(stdout_data)
            if stderr_data:
                log_file.write("\n--- stderr ---\n")
                log_file.write(stderr_data)
            log_file.flush()

            if process.returncode == 0:
                print(f"Successfully ran {module_name}.")
                log_file.write(f"\nSuccessfully ran {module_name}.\n")
                return True
            else:
                print(f"Error: {module_name} ({base_script_name}) failed. Return code: {process.returncode}. Check log: {log_file_path}")
                log_file.write(f"\nError: {module_name} ({base_script_name}) failed. Return code: {process.returncode}.\n")
                if script_type == "python" and conda_exe:
                    full_output = (stdout_data or "") + (stderr_data or "")
                    if "No such environment" in full_output:
                         log_file.write(f"\nPotential issue: The Conda environment '{target_conda_env_python}' might not exist or is not accessible.\n")
                return False

    except FileNotFoundError as e:
        error_msg = f"Fatal error running {module_name} ({base_script_name}): Command not found ('{e.filename}' -> {e.strerror}). Ensure Conda/Rscript/Python is in PATH or paths in params.json are correct."
        print(error_msg)
        try:
            with open(log_file_path, 'a', encoding='utf-8', errors='replace') as log_file_append:
                log_file_append.write(f"\nPython orchestrator fatal error: Command not found - {e.filename}: {e.strerror}\n")
        except Exception: 
            pass
        return False
    except Exception as e:
        error_msg = f"Fatal error running {module_name} ({base_script_name}): {e}"
        print(error_msg)
        try:
            with open(log_file_path, 'a', encoding='utf-8', errors='replace') as log_file_append:
                log_file_append.write(f"\nPython orchestrator fatal error: {e}\n")
        except Exception:
            pass
        return False

def main():
    if len(sys.argv) != 2:
        print("Usage: python run_impact_sc_pipeline.py <path_to_params.json>")
        sys.exit(1)

    params_path = sys.argv[1]
    if not os.path.exists(params_path):
        print(f"Error: Parameters file not found at {params_path}")
        sys.exit(1)
    
    _params_temp = {}
    try:
        with open(params_path, 'r', encoding='utf-8') as f:
            _params_temp = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: Could not parse parameters file {params_path}. Invalid JSON: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading parameters file {params_path}: {e}")
        sys.exit(1)
    
    _params_temp['__params_file_path__'] = params_path 
    params = _params_temp

    print("--- Starting IMPACT-sc Pipeline Execution ---")
    print(f"Parameters loaded from: {params_path}")

    current_conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    expected_conda_env = "impact_sc"
    if current_conda_env != expected_conda_env:
        print("-" * 70)
        print(f"WARNING: You appear to be running this script from an unexpected Conda environment.")
        print(f"Current environment: '{current_conda_env if current_conda_env else 'Not in a Conda environment or CONDA_DEFAULT_ENV not set.'}'")
        print(f"Expected environment: '{expected_conda_env}'")
        print(f"To ensure all Python modules run correctly, please activate the '{expected_conda_env}' environment first:")
        print(f"  conda activate {expected_conda_env}")
        print(f"Then re-run the pipeline: ")
        print(f"  python {' '.join(sys.argv)}")
        print("The pipeline will attempt to use 'conda run' for Python scripts, which might mitigate this,")
        print("but it's best practice to run the main orchestrator script from the target environment.")
        print("-" * 70)

    output_directory = params.get("output_directory")
    if not output_directory:
        print("Error: 'output_directory' not specified in parameters file.")
        sys.exit(1)
    print(f"Output directory: {output_directory}")

    try:
        os.makedirs(output_directory, exist_ok=True)
    except OSError as e:
        print(f"Error: Could not create output directory {output_directory}: {e}")
        sys.exit(1)

    input_r_scripts_dir = params.get("input_r_scripts_dir")
    input_python_scripts_dir = params.get("input_python_scripts_dir")

    if not input_r_scripts_dir:
        print("Warning: 'input_r_scripts_dir' not specified. Assuming R scripts are relative to CWD or run_script's CWD.")
    if not input_python_scripts_dir:
        print("Warning: 'input_python_scripts_dir' not specified. Assuming Python scripts are relative to CWD or run_script's CWD.")

    script_map = {
        "01_data_processing": ("R", "01_data_processing.R"),
        "02a_harmony_c2s_prep": ("R", "02a_harmony_c2s_prep.R"),
        "02b_c2s": ("python", "02b_c2s.py"),
        "02c_load_c2s_result": ("R", "02c_load_c2s_result.R"),
        "03_cell_type_annotation": ("R", "03_cell_type_annotation.R"),
        "04a_basic_visualization": ("R", "04a_basic_visualization.R"),
        "04b_DE_gsea": ("R", "04b_DE_gsea.R"),
        "04c_decoupler": ("R", "04c_decoupler.R"),
        "04d_ucell_scores": ("R", "04d_ucell_scores.R"),
        "04e_pseudotime": ("R", "04e_pseudotime.R"),
        # "04f_query_projection": ("R", "04f_query_projection.R")  # TEMPORARILY DISABLED
    }

    selected_modules = params.get("selected_modules", [])
    if not selected_modules:
        print("Warning: No modules selected in 'selected_modules' parameter.")

    for module_name in selected_modules:
        # Temporarily skip 04f_query_projection module
        if module_name == "04f_query_projection":
            print(f"INFO: Module '{module_name}' is temporarily disabled. Skipping...")
            continue
            
        if module_name in script_map:
            script_type, script_file_name = script_map[module_name]

            script_base_dir = ""
            if script_type == "R":
                script_base_dir = input_r_scripts_dir if input_r_scripts_dir else os.getcwd()
            elif script_type == "python":
                script_base_dir = input_python_scripts_dir if input_python_scripts_dir else os.getcwd()
            
            script_full_path = os.path.join(script_base_dir, script_file_name)

            if not os.path.exists(script_full_path):
                pipeline_script_dir = os.path.dirname(os.path.abspath(__file__))
                potential_path_relative_to_pipeline = os.path.join(pipeline_script_dir, script_file_name)
                
                if os.path.exists(potential_path_relative_to_pipeline):
                    script_full_path = potential_path_relative_to_pipeline
                    print(f"Info: Found script {script_file_name} relative to pipeline script location: {script_full_path}")
                else:
                    potential_path_in_scripts_subdir = os.path.join(pipeline_script_dir, "scripts_AI", script_file_name)
                    if os.path.exists(potential_path_in_scripts_subdir):
                        script_full_path = potential_path_in_scripts_subdir
                        print(f"Info: Found script {script_file_name} in 'scripts_AI' subdirectory relative to pipeline script: {script_full_path}")
                    else:
                        print(f"Error: Script file {script_file_name} (expected at {script_full_path}, {potential_path_relative_to_pipeline}, or {potential_path_in_scripts_subdir}) for module {module_name} not found. Skipping.")
                        continue
            
            if not run_script(script_full_path, script_type, params, module_name): 
                print(f"Pipeline execution stopped due to error in module {module_name}.")
                sys.exit(1)
        else:
            print(f"Warning: Module '{module_name}' is selected but not defined in the script map. Skipping.")

    print("--- IMPACT-sc Pipeline Execution Finished ---")

if __name__ == "__main__":
    main()
