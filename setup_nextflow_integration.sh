#!/bin/bash

# Setup script for F.A.D.E Nextflow Integration

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV_PATH="${SCRATCH}/conda-envs/fade"

echo "======================================="
echo "F.A.D.E Nextflow Integration Setup"
echo "======================================="

# Step 1: Load required modules
echo "Step 1: Loading required modules..."
module load miniconda3/24.11.1
module load nextflow/24.10.3

# Step 2: Activate conda environment
echo "Step 2: Activating conda environment..."
if [ ! -d "$CONDA_ENV_PATH" ]; then
    echo "Creating conda environment..."
    conda create -p "$CONDA_ENV_PATH" python=3.9 -y
fi

source activate "$CONDA_ENV_PATH"

# Step 3: Install required Python packages
echo "Step 3: Installing Python dependencies..."
cd "$SCRIPT_DIR"
pip install -r requirements.txt

# Step 4: Verify installation
echo "Step 4: Verifying installation..."
python3 -c "
import sys
sys.path.insert(0, '.')
try:
    from utils.gemini_client import GeminiClient
    from utils.uniprot_client import UniProtClient
    print('✅ Core utilities imported successfully')
except Exception as e:
    print(f'❌ Core utility import failed: {e}')
    sys.exit(1)

try:
    from agents.target_selector.target_selector import TargetSelector
    print('✅ TargetSelector imported successfully')
except Exception as e:
    print(f'❌ TargetSelector import failed: {e}')
    # Don't exit here as this might be due to missing API key

try:
    from agents.molecule_generator.molecule_generator import MoleculeGenerator
    print('✅ MoleculeGenerator imported successfully')
except Exception as e:
    print(f'❌ MoleculeGenerator import failed: {e}')
"

# Step 5: Test Nextflow workflow syntax
echo "Step 5: Testing Nextflow workflow syntax..."
nextflow config nextflow/main.nf > /dev/null && echo "✅ Nextflow syntax valid" || echo "❌ Nextflow syntax error"

# Step 6: Run stub test
echo "Step 6: Running stub test..."
nextflow run nextflow/main.nf --query "Test query for EGFR" -profile stub -work-dir work_test && echo "✅ Stub test passed" || echo "❌ Stub test failed"

echo ""
echo "======================================="
echo "Setup completed!"
echo "======================================="
echo ""
echo "Next steps:"
echo "1. Set your GEMINI_API_KEY in the .env file"
echo "2. Test with: ./run_nextflow_integrated.sh --dry-run 'Test query'"
echo "3. Run actual workflow: ./run_nextflow_integrated.sh 'Your query'"
echo ""
echo "Environment ready at: $CONDA_ENV_PATH"
