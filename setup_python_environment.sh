#!/bin/bash

echo "==================================================================="
echo "IMPACT-sc Python Environment Setup"
echo "==================================================================="
echo "This script creates a conda environment and installs Python packages"
echo "Prerequisites: Conda/Miniconda must be installed"
echo ""

ENV_NAME="impact_sc"
PYTHON_VERSION="3.9"
LOG_FILE="python_env_setup.log"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: conda not found in PATH"
    echo "Please install Miniconda or Anaconda:"
    echo "- Windows: https://docs.conda.io/en/latest/miniconda.html#windows-installers"
    echo "- macOS Intel: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg"
    echo "- macOS Apple Silicon: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.pkg"
    exit 1
fi

echo "✅ Found conda: $(conda --version)"

# Detect platform
if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="macOS"
    ARCH=$(uname -m)
    if [[ "$ARCH" == "arm64" ]]; then
        echo "🍎 Detected: macOS Apple Silicon"
    else
        echo "💻 Detected: macOS Intel"
    fi
elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    PLATFORM="Windows"
    echo "🪟 Detected: Windows"
else
    PLATFORM="Linux"
    echo "🐧 Detected: Linux"
fi

echo "Log will be saved to: $LOG_FILE"
echo ""

# Function to log and execute commands
log_and_run() {
    echo "Running: $*" | tee -a "$LOG_FILE"
    "$@" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# Step 1: Create/Update conda environment
echo "=== Step 1: Setting up conda environment '$ENV_NAME' ==="

if conda env list | grep -q "^$ENV_NAME "; then
    echo "⚠️  Environment '$ENV_NAME' already exists."
    read -p "Do you want to remove and recreate it? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        log_and_run conda env remove -n "$ENV_NAME" -y
    else
        echo "Updating existing environment..."
    fi
fi

if ! conda env list | grep -q "^$ENV_NAME "; then
    echo "Creating new conda environment: $ENV_NAME"
    log_and_run conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    if [ $? -ne 0 ]; then
        echo "❌ Failed to create conda environment"
        exit 1
    fi
fi

echo "✅ Environment '$ENV_NAME' ready"

# Step 2: Install core scientific packages via conda
echo ""
echo "=== Step 2: Installing core scientific packages via conda ==="

echo "Activating environment and installing packages..."
log_and_run conda install -n "$ENV_NAME" -c conda-forge -c bioconda \
    pandas numpy scipy matplotlib seaborn \
    scikit-learn jupyter jupyterlab \
    scanpy anndata openpyxl h5py \
    notebook ipywidgets -y

if [ $? -ne 0 ]; then
    echo "❌ Failed to install core packages via conda"
    exit 1
fi

echo "✅ Core scientific packages installed"

# Step 3: Install specialized packages via pip
echo ""
echo "=== Step 3: Installing specialized packages via pip ==="

# Activate environment and install pip packages
echo "Installing specialized packages via pip..."

# Create a temporary script for pip installations within the conda environment
PIP_INSTALL_SCRIPT="temp_pip_install.sh"
cat << EOF > "$PIP_INSTALL_SCRIPT"
#!/bin/bash
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo "Python version: \$(python --version)"
echo "Pip version: \$(pip --version)"

# Install cell2sentence and related packages
echo "Installing cell2sentence..."
pip install cell2sentence

# Install transformer-related packages
echo "Installing transformers and torch..."
pip install transformers torch

# Install additional utility packages
echo "Installing additional packages..."
pip install umap-learn leidenalg python-igraph

# Verify critical installations
echo ""
echo "=== Verification ==="
python -c "
try:
    import pandas as pd
    print('✅ pandas:', pd.__version__)
except ImportError as e:
    print('❌ pandas:', e)

try:
    import numpy as np
    print('✅ numpy:', np.__version__)
except ImportError as e:
    print('❌ numpy:', e)

try:
    import scanpy as sc
    print('✅ scanpy:', sc.__version__)
except ImportError as e:
    print('❌ scanpy:', e)

try:
    import cell2sentence
    print('✅ cell2sentence: available')
except ImportError as e:
    print('❌ cell2sentence:', e)

try:
    import transformers
    print('✅ transformers:', transformers.__version__)
except ImportError as e:
    print('❌ transformers:', e)

print('🎉 Python environment verification complete!')
"
EOF

chmod +x "$PIP_INSTALL_SCRIPT"
log_and_run bash "$PIP_INSTALL_SCRIPT"
PIP_EXIT_CODE=$?

# Clean up
rm -f "$PIP_INSTALL_SCRIPT"

if [ $PIP_EXIT_CODE -ne 0 ]; then
    echo "❌ Failed to install pip packages"
    exit 1
fi

# Step 4: Final verification and summary
echo ""
echo "=== Step 4: Final Environment Test ==="

# Test the complete environment
TEST_SCRIPT="temp_test_env.sh"
cat << EOF > "$TEST_SCRIPT"
#!/bin/bash
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo "Testing complete IMPACT-sc Python environment..."
echo ""

echo "Environment information:"
echo "Python: \$(python --version)"
echo "Location: \$(which python)"
echo "Conda env: \$CONDA_DEFAULT_ENV"
echo ""

python -c "
import sys
print('Python executable:', sys.executable)
print('Python path:', sys.path[0])
print('')

# Test all critical packages
packages_to_test = [
    'pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn',
    'sklearn', 'scanpy', 'anndata', 'cell2sentence', 
    'transformers', 'jupyter', 'jupyterlab'
]

failed = []
for pkg in packages_to_test:
    try:
        __import__(pkg)
        print(f'✅ {pkg}: OK')
    except ImportError:
        print(f'❌ {pkg}: FAILED')
        failed.append(pkg)

if failed:
    print(f'\\n⚠️  Failed packages: {failed}')
    sys.exit(1)
else:
    print('\\n🎉 All packages imported successfully!')
    print('\\nYour IMPACT-sc Python environment is ready!')
"
EOF

chmod +x "$TEST_SCRIPT"
log_and_run bash "$TEST_SCRIPT"
TEST_EXIT_CODE=$?

# Clean up
rm -f "$TEST_SCRIPT"

echo ""
echo "==================================================================="
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "✅ IMPACT-sc Python environment setup completed successfully!"
    echo ""
    echo "🚀 Next steps:"
    echo "1. Activate the environment: conda activate $ENV_NAME"
    echo "2. Run the interactive setup: python interactive_setup.py"
    echo ""
    echo "📝 Environment summary:"
    echo "- Environment name: $ENV_NAME"
    echo "- Python version: $PYTHON_VERSION"
    echo "- Platform: $PLATFORM"
    echo "- Log file: $LOG_FILE"
else
    echo "❌ Python environment setup encountered errors"
    echo "Check the log file for details: $LOG_FILE"
    echo ""
    echo "Common solutions:"
    echo "1. Update conda: conda update conda"
    echo "2. Clear conda cache: conda clean --all"
    echo "3. Check internet connection"
    echo "4. Try creating environment manually:"
    echo "   conda create -n $ENV_NAME python=$PYTHON_VERSION -y"
    echo "   conda activate $ENV_NAME"
    echo "   pip install cell2sentence scanpy transformers"
fi
echo "==================================================================="

exit $TEST_EXIT_CODE 