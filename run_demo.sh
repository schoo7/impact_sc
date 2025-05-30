#!/bin/bash

# IMPACT-sc Demo Launcher
# Automatically sets up and runs the pipeline in demo mode

set -euo pipefail

echo "🚀 IMPACT-sc Demo Launcher"
echo "=========================="

# Check if data is downloaded
if [[ ! -d "data/demo" ]]; then
    echo "❌ Demo data not found. Running data download first..."
    ./download_data.sh
fi

# Check if conda environment exists
if ! conda env list | grep -q "impact_sc"; then
    echo "❌ Conda environment 'impact_sc' not found. Please run './install_dependencies.sh' first."
    exit 1
fi

# Run interactive setup in demo mode
echo "🎯 Setting up demo configuration..."
conda run -n impact_sc python interactive_setup.py

# Check if demo parameters were created
if [[ -f "demo_output/impact_sc_params.json" ]]; then
    echo "✅ Demo configuration created successfully!"
    echo "🏃 Running IMPACT-sc pipeline in demo mode..."
    conda run -n impact_sc python run_impact_sc_pipeline.py demo_output/impact_sc_params.json
    echo "🎉 Demo completed! Check results in demo_output/"
else
    echo "❌ Demo configuration failed. Please check the setup."
    exit 1
fi 