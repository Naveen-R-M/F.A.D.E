#!/bin/bash

# Fix conda environment by installing missing dependencies

echo "🔧 Fixing conda environment dependencies..."

# Load modules
module purge
module load miniconda3/24.11.1

# Activate conda environment
CONDA_ENV_PATH="${SCRATCH}/conda-envs/fade"
source activate "$CONDA_ENV_PATH"

echo "Installing missing dependencies..."

# Install the missing google-generativeai package
pip install google-generativeai>=0.4.0

# Verify installation
echo "Verifying installation..."
python3 -c "
import google.generativeai as genai
print('✅ google-generativeai installed successfully')

# Test other imports
try:
    import sys
    sys.path.insert(0, '.')
    from utils.gemini_client import GeminiClient
    print('✅ GeminiClient imports successfully')
except Exception as e:
    print(f'❌ GeminiClient import failed: {e}')
"

echo "✅ Conda environment fixed!"
echo ""
echo "Now you can run:"
echo "  source setup_nextflow_env.sh"
echo "  ./run_nextflow_integrated.sh 'Your drug discovery query'"
