# IMPACT-sc: Integrated Multi-Pipeline Analysis and Characterization of Single-Cell Data

<p align="center">
  <img src="impact_sc_logo.png" alt="IMPACT-sc Logo" width="800"/>
</p>

**IMPACT-sc** is a modular pipeline for comprehensive single-cell RNA sequencing (scRNA-seq) data analysis, integrating R and Python scripts for everything from data processing to advanced downstream analyses.

---

## 🚀 Quick Start

**Choose your platform:**

<details>
<summary><strong>🪟 Windows Users - Click to expand</strong></summary>

### Windows Quick Setup
```bash
# 1. Install prerequisites: Git Bash, R + Rtools, Conda
# 2. Open Git Bash as Administrator
cd /path/to/impact_sc
bash install_r_packages_windows.sh
bash setup_python_environment.sh
conda activate impact_sc
python interactive_setup.py
```
</details>

<details>
<summary><strong>🍎 macOS Users - Click to expand</strong></summary>

### macOS Quick Setup  
```bash
# 1. Install prerequisites: Xcode CLI, Homebrew, R, Conda
xcode-select --install
brew install r cmake pkg-config gfortran
cd /path/to/impact_sc
./install_r_packages_mac.sh
bash setup_python_environment.sh
conda activate impact_sc
python interactive_setup.py
```
</details>

---

## 📋 Overview

IMPACT-sc consists of three main components:
1. **Dependency Installation** - Set up R and Python environments
2. **Interactive Configuration** - Generate analysis parameters  
3. **Pipeline Execution** - Run selected analysis modules

### Key Features:
- **Data Processing**: QC, filtering, normalization
- **Batch Correction**: Harmony integration  
- **Cell Type Annotation**: Seurat, SingleR, Cell2Sentence
- **Visualization**: UMAP, tSNE, feature plots
- **Differential Expression**: DGE and GSEA analysis
- **Pathway Analysis**: DecoupleR, PROGENy, UCell
- **Advanced Analysis**: Pseudotime, query projection

<p align="center">
  <img src="overview_image.png" alt="IMPACT-sc Pipeline Overview" width="800"/>
  <br/>
  <em>Figure: Overview of the IMPACT-sc pipeline workflow</em>
</p>

**Supported:** Human & Mouse | Windows & macOS (including Apple Silicon) | R + Python integration

---

## 🛠️ Installation Guide

### Step 1: Prerequisites

<details>
<summary><strong>🪟 Windows Prerequisites</strong></summary>

#### **1.1 Install Git Bash**
- Download from: https://git-scm.com/download/win
- **IMPORTANT**: Run Git Bash as Administrator for installations

#### **1.2 Install R and Rtools**
1. **Install R** from: https://cran.r-project.org/bin/windows/base/
2. **Install Rtools** (CRITICAL): https://cran.r-project.org/bin/windows/Rtools/
   - Choose version matching your R (e.g., Rtools43 for R 4.3.x)
   - **Check "Add to PATH" during installation**

#### **1.3 Install Conda**
- **Miniconda** (recommended): https://docs.conda.io/en/latest/miniconda.html#windows-installers
- Choose "Miniconda3 Windows 64-bit"

#### **1.4 Windows Locale Settings (IMPORTANT)**
1. Control Panel > Region > Administrative tab
2. Under 'Language for non-Unicode programs', click 'Change system locale...'
3. Ensure **'Beta: Use Unicode UTF-8' is UNCHECKED**
4. Set locale to **'English (United States)'**
5. **RESTART your computer** (required for R packages)

**Verify Installation:**
```bash
# In Git Bash
Rscript --version
make --version  # Should work if Rtools installed correctly
conda --version
```
</details>

<details>
<summary><strong>🍎 macOS Prerequisites</strong></summary>

#### **1.1 Install Xcode Command Line Tools**
```bash
xcode-select --install
```

#### **1.2 Install Homebrew**
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### **1.3 Install R**
```bash
# Option A: Homebrew (recommended)
brew install r

# Option B: Download from https://cran.r-project.org/bin/macosx/
```

#### **1.4 Install Build Tools**
```bash
brew install cmake pkg-config gfortran

# For Apple Silicon only
brew install llvm libomp
```

#### **1.5 Install Conda**
- **Miniconda**: https://docs.conda.io/en/latest/miniconda.html#macos-installers
- **Apple Silicon**: Choose `Miniconda3 macOS Apple M1 64-bit pkg`
- **Intel**: Choose `Miniconda3 macOS Intel x86 64-bit pkg`

**Verify Installation:**
```bash
Rscript --version
gcc --version
conda --version
```
</details>

#### **Optional: Ollama AI Enhancement**
```bash
# Download from: https://ollama.ai/download
ollama --version
ollama pull gemma2:9b  # Install recommended model
```

---

### Step 2: Platform-Specific Package Installation

<details>
<summary><strong>🪟 Windows R Package Installation</strong></summary>

1. **Open Git Bash as Administrator** (Right-click > "Run as administrator")
2. **Navigate to project directory:**
   ```bash
   cd /path/to/impact_sc
   ```
3. **Run Windows R installer:**
   ```bash
   bash install_r_packages_windows.sh
   ```

**What this installs:**
- Bioconductor packages (SingleR, Seurat, etc.)
- CRAN packages (ggplot2, dplyr, etc.)  
- GitHub packages (SeuratExtend, CARD, etc.)

**Time:** 30-60 minutes | **Log:** `r_package_install_windows.log`
</details>

<details>
<summary><strong>🍎 macOS R Package Installation</strong></summary>

```bash
# Navigate to project directory
cd /path/to/impact_sc

# Run macOS R installer
chmod +x install_r_packages_mac.sh
./install_r_packages_mac.sh
```

**What this installs:**
- Platform-optimized compilation for Apple Silicon
- Bioconductor packages (SingleR, Seurat, etc.)
- CRAN packages (ggplot2, dplyr, etc.)  
- GitHub packages (SeuratExtend, CARD, etc.)

**Time:** 30-60 minutes (longer on Apple Silicon) | **Log:** `r_package_install_mac.log`
</details>

---

### Step 3: Python Environment Setup

**All Platforms:**
```bash
# Create and set up Python environment
bash setup_python_environment.sh
```

**What this does:**
- Creates `impact_sc` conda environment with Python 3.9
- Installs scientific packages (pandas, numpy, scanpy, etc.)
- Installs specialized packages (cell2sentence, transformers)

**Time:** 15-30 minutes | **Log:** `python_env_setup.log`

---

### Step 4: Verification

**Test R Installation:**
```bash
Rscript -e "
library(Seurat)
library(SingleR) 
library(ggplot2)
cat('✅ R packages loaded successfully!\n')
cat('Seurat version:', as.character(packageVersion('Seurat')), '\n')
"
```

**Test Python Environment:**
```bash
conda activate impact_sc
python -c "
import pandas as pd
import scanpy as sc
import cell2sentence
print('✅ All Python packages installed successfully!')
print(f'Scanpy version: {sc.__version__}')
"
```

---

## 🔧 Configuration and Usage

### Step 1: Interactive Setup

```bash
# Activate environment
conda activate impact_sc

# Run interactive setup
python interactive_setup.py
```

**Configuration includes:**

<details>
<summary><strong>Platform-Specific Path Examples</strong></summary>

**🪟 Windows Paths:**
```bash
R executable: "C:\Program Files\R\R-4.3.0\bin\x64\Rscript.exe"
Scripts: "C:\path\to\impact_sc\scripts_AI"
Output: "C:\Users\YourName\Documents\impact_output"
```

**🍎 macOS Paths:**
```bash
# Apple Silicon
R executable: "/opt/homebrew/bin/Rscript"
# Intel Mac  
R executable: "/usr/local/bin/Rscript"
# CRAN R (both)
R executable: "/Library/Frameworks/R.framework/Resources/bin/Rscript"

Scripts: "/Users/yourusername/impact_sc/scripts_AI"
Output: "/Users/yourusername/impact_sc/output"
```
</details>

### Step 2: Run Pipeline

```bash
python run_impact_sc_pipeline.py /path/to/impact_sc_params.json
```

---

## 📊 Available Analysis Modules

| Module | Description | Requirements |
|--------|-------------|--------------|
| **01_data_processing** | QC, filtering, normalization | Raw scRNA-seq data |
| **02a_harmony_c2s_prep** | Batch correction prep | Processed data |
| **02b_c2s** | Cell2Sentence analysis | H5AD file, C2S model |
| **02c_load_c2s_result** | Load C2S results | C2S output |
| **03_cell_type_annotation** | Cell type annotation | SingleR reference |
| **04a_basic_visualization** | UMAP, tSNE, plots | Processed data |
| **04b_DE_gsea** | Differential expression | Annotated data |
| **04c_decoupler** | Pathway analysis | CollecTRI/PROGENy files |
| **04d_ucell_scores** | Gene signatures | MSigDB data |
| **04e_pseudotime** | Trajectory analysis | Start cell barcode |
| **04f_query_projection** | Query mapping | Query RDS file |

---

## 🐛 Troubleshooting

<details>
<summary><strong>🪟 Windows-Specific Issues</strong></summary>

| Problem | Solution |
|---------|----------|
| **Permission denied** | Run Git Bash as Administrator |
| **R packages won't compile** | Install/reinstall Rtools, check PATH |
| **Conda not found** | Add conda to PATH, restart terminal |
| **Locale errors** | Set locale to English (US), restart computer |
| **Make command not found** | Install Rtools with PATH option checked |

**Debugging:**
```bash
# Check R and Rtools
Rscript -e "Sys.which('make')"  # Should show path

# Check environment
conda info --envs  # Should show impact_sc
```

**Key Requirements:**
- ✅ Administrator privileges for Git Bash
- ✅ Rtools installed with PATH
- ✅ English (US) locale settings
- ✅ Computer restarted after locale change
</details>

<details>
<summary><strong>🍎 macOS-Specific Issues</strong></summary>

| Problem | Solution |
|---------|----------|
| **R packages won't compile** | Install Xcode CLI: `xcode-select --install` |
| **Permission errors** | Fix R library permissions |
| **Conda not found** | Add conda to PATH in `~/.zshrc` |
| **Build tool errors** | `brew install cmake gfortran` |
| **Apple Silicon issues** | Try Rosetta mode: `arch -x86_64 R` |

**Apple Silicon Environment Variables:**
```bash
# Add to ~/.zshrc if needed
export LDFLAGS="-L/opt/homebrew/lib"
export CPPFLAGS="-I/opt/homebrew/include"
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig"
```

**Debugging:**
```bash
# Check architecture
uname -m  # arm64 = Apple Silicon, x86_64 = Intel

# Check R paths
which Rscript
R --slave -e ".libPaths()"

# Check build tools
xcode-select -p  # Should show Xcode path
```

**Key Requirements:**
- ✅ Xcode Command Line Tools installed
- ✅ Homebrew installed and working  
- ✅ R accessible via command line
- ✅ Conda environment properly activated
</details>

### **Common Log Files:**
- `r_package_install_windows.log` / `r_package_install_mac.log` - R package installation
- `python_env_setup.log` - Python environment setup  
- `*_log.txt` in output directory - Individual module logs

---

## 📁 Project Structure

```
impact_sc/
├── README.md                      # This comprehensive guide
├── scripts_AI/                    # Analysis modules
├── install_r_packages_windows.sh  # Windows R installer
├── install_r_packages_mac.sh      # macOS R installer  
├── setup_python_environment.sh    # Python environment setup
├── interactive_setup.py           # Configuration script
├── run_impact_sc_pipeline.py     # Main pipeline
└── output/                        # Results (created during setup)
    ├── impact_sc_params.json     # Configuration file
    ├── *_log.txt                 # Module logs
    └── results/                   # Analysis outputs
```

---

## 🔬 Technical Details

### **System Requirements:**
- **R** (≥4.2) + Bioconductor packages
- **Python** (3.9) + scientific stack  
- **System tools** (build tools, conda)

### **Supported Features:**
- ✅ **Species**: Human and Mouse
- ✅ **Platforms**: Windows and macOS (including Apple Silicon)
- ✅ **Environments**: Conda virtual environments
- ✅ **Integration**: R + Python seamless workflow
- ✅ **Reproducibility**: JSON parameter configuration
- ✅ **Modularity**: Select only needed analysis steps

### **Performance Notes:**
- **Installation time**: 1-2 hours total
- **R packages**: 30-60 minutes
- **Python environment**: 15-30 minutes  
- **Pipeline execution**: Varies by data size and modules

---

## 🆘 Getting Help

1. **Check platform-specific troubleshooting sections above**
2. **Review log files** for detailed error messages
3. **Verify prerequisites** are properly installed
4. **Test each step individually** using verification commands

### **Cross-Platform Notes:**
- JSON parameter files work on both platforms
- Module functionality is identical across platforms
- File paths use platform-appropriate separators

---

## 📄 Citation

If you use IMPACT-sc in your research, please cite:

```
[Citation information to be added]
```

## 🤝 Contributing

Contributions welcome! Please ensure compatibility with both Windows and macOS platforms.

---

**⚡ Total setup time: ~1-2 hours | Get started by expanding your platform section above!** 