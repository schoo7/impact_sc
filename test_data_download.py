#!/usr/bin/env python3
"""
Test script to verify data download components
"""

import os
import sys

def test_model_download():
    """Test if we can load the transformers library and access the model"""
    try:
        from transformers import pipeline, AutoTokenizer, AutoModelForCausalLM
        print("✅ Transformers library available")
        
        model_name = "vandijklab/C2S-Pythia-410m-cell-type-prediction"
        print(f"Testing model access: {model_name}")
        
        # Just test if we can create the pipeline (without downloading)
        # This will work if model is cached, otherwise will show what would happen
        try:
            tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=True)
            print("✅ Model already cached locally")
        except:
            print("ℹ️ Model not cached - would download on first use")
        
        return True
    except ImportError as e:
        print(f"❌ Missing dependency: {e}")
        return False
    except Exception as e:
        print(f"⚠️ Model test warning: {e}")
        return True

def test_r_packages():
    """Test if R packages are available"""
    import subprocess
    
    try:
        # Test celldex availability
        result = subprocess.run([
            'Rscript', '-e', 
            'if(!require(celldex, quietly=TRUE)) quit(status=1); cat("celldex available\\n")'
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ celldex package available")
            return True
        else:
            print("❌ celldex package not available")
            print("Install with: ./install_dependencies.sh")
            return False
            
    except Exception as e:
        print(f"❌ Error testing R packages: {e}")
        return False

def test_download_tools():
    """Test if download tools are available"""
    tools = ['wget', 'curl', 'tar']
    available = []
    
    for tool in tools:
        if os.system(f"command -v {tool} >/dev/null 2>&1") == 0:
            available.append(tool)
            print(f"✅ {tool} available")
        else:
            print(f"❌ {tool} not available")
    
    return len(available) >= 2  # Need at least wget/curl and tar

if __name__ == "__main__":
    print("Testing IMPACT-sc data download prerequisites...")
    print("=" * 50)
    
    tests = [
        ("Download tools", test_download_tools),
        ("Model download", test_model_download), 
        ("R packages", test_r_packages)
    ]
    
    all_passed = True
    for name, test_func in tests:
        print(f"\nTesting {name}:")
        passed = test_func()
        all_passed = all_passed and passed
    
    print("\n" + "=" * 50)
    if all_passed:
        print("✅ All tests passed! Ready to download data.")
    else:
        print("❌ Some tests failed. Run install_dependencies.sh first.")
        sys.exit(1) 