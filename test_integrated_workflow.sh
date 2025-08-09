#!/bin/bash

# Test script for F.A.D.E Integrated Nextflow Workflow

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "======================================="
echo "F.A.D.E Integration Test Suite"
echo "======================================="

# Test 1: Dry run with stub processes
echo "Test 1: Running dry run with stub processes..."
./run_nextflow_integrated.sh --dry-run "Test molecules targeting EGFR"

if [ $? -eq 0 ]; then
    echo "✅ Dry run test passed"
else
    echo "❌ Dry run test failed"
    exit 1
fi

# Test 2: Check if wrapper scripts are executable
echo "Test 2: Checking wrapper scripts..."
for script in nextflow/bin/run_*.py; do
    if [ -x "$script" ]; then
        echo "✅ $script is executable"
    else
        echo "❌ $script is not executable"
        exit 1
    fi
done

# Test 3: Validate Python imports (if conda env exists)
echo "Test 3: Validating Python imports..."
if [ -d "${SCRATCH}/conda-envs/fade" ]; then
    source activate ${SCRATCH}/conda-envs/fade
    python3 -c "
import sys
sys.path.insert(0, '${SCRIPT_DIR}')
try:
    from agents.target_selector.target_selector import TargetSelector
    from agents.structure_predictor.structure_predictor import StructurePredictor
    from agents.molecule_generator.molecule_generator import MoleculeGenerator
    print('✅ All agent imports successful')
except Exception as e:
    print(f'❌ Import failed: {e}')
    sys.exit(1)
"
else
    echo "⚠️ Conda environment not found, skipping import test"
fi

echo "======================================="
echo "Integration tests completed!"
echo "======================================="

# Show available commands
echo ""
echo "Available commands:"
echo "1. Dry run test:     ./run_nextflow_integrated.sh --dry-run 'YOUR_QUERY'"
echo "2. Debug mode:       ./run_nextflow_integrated.sh --debug --trace 'YOUR_QUERY'"
echo "3. Full pipeline:    ./run_nextflow_integrated.sh 'YOUR_QUERY'"
echo ""
echo "Example:"
echo "./run_nextflow_integrated.sh 'Find molecules targeting KRAS G12D with good BBB permeability'"
