#!/usr/bin/env python3
"""
Test script to verify data download components.
This version checks for core system tools and Python model dependencies.
"""

import os
import sys
import shutil

def test_system_tools():
    """Checks for essential command-line tools that should be in the system's PATH."""
    print("Checking for Git...")
    if shutil.which("git"):
        print("✅ Git is available.")
        return True
    else:
        print("❌ Git is not available in the system's PATH.")
        print("   Please install Git and ensure it can be called from your terminal.")
        return False

def test_model_dependencies():
    """Test if we can load the transformers library and access the model."""
    try:
        # Check if the 'transformers' library is installed
        from transformers import pipeline, AutoTokenizer, AutoModelForCausalLM
        print("✅ Transformers library is available.")

        model_name = "vandijklab/C2S-Pythia-410m-cell-type-prediction"
        print(f"Verifying access to model: {model_name}")

        try:
            # Check only if the model is cached locally to avoid network downloads
            tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=True)
            print("✅ Model appears to be cached locally.")
        except EnvironmentError: # Catches OSError, etc., when the model is not found locally
            print("ℹ️  Model is not cached. It will be downloaded by the main pipeline on first use.")

        return True
    except ImportError as e:
        print(f"❌ Missing critical Python dependency for model testing: {e}")
        print("   Please ensure the 'impact_sc' conda environment is activated and its dependencies are installed.")
        return False
    except Exception as e:
        # Catch other potential exceptions during the test
        print(f"⚠️  An unexpected warning occurred during model dependency test: {e}")
        return True # Return True as this is a warning, not a critical failure

if __name__ == "__main__":
    print("Testing IMPACT-sc Prerequisites...")
    print("=" * 50)

    # Updated test list to include a cross-platform tool check
    tests = [
        ("System Tools", test_system_tools),
        ("Model Dependencies", test_model_dependencies)
    ]

    all_passed = True
    for name, test_func in tests:
        print(f"\n--- Testing {name} ---")
        passed = test_func()
        if not passed:
            all_passed = False

    print("\n" + "=" * 50)
    if all_passed:
        print("✅ All prerequisite checks passed! You should be ready to proceed.")
    else:
        print("❌ Some prerequisite checks failed. Please review the messages above.")
        print("   Ensure you have run 'install_dependencies.sh' or equivalent and have activated the 'impact_sc' environment.")
        sys.exit(1)
