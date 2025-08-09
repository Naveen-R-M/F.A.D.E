#!/bin/bash

# Test Nextflow workflow with stub processes

set -e

echo "🧪 Testing F.A.D.E Nextflow Integration with Stubs"
echo "=================================================="

# Setup environment
echo "Setting up environment..."
unset JAVA_HOME JAVA_PATH JDK_HOME
module purge
module load miniconda3/24.11.1
module load nextflow/24.10.3

# Verify Nextflow
echo "Nextflow version:"
nextflow -version

# Change to nextflow directory
cd nextflow

# Test workflow syntax
echo ""
echo "Testing workflow syntax..."
nextflow config main.nf > /dev/null && echo "✅ Workflow syntax is valid" || { echo "❌ Workflow syntax error"; exit 1; }

# Run stub workflow
echo ""
echo "Running stub workflow..."
nextflow run main.nf \
    --query "Test molecules targeting EGFR for lung cancer" \
    --max_molecules 5 \
    --output_dir ../test_results_$(date +'%Y%m%d_%H%M%S') \
    -profile stub \
    -work-dir ../work_test

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Stub workflow completed successfully!"
    echo "✅ Nextflow integration is working"
    echo ""
    echo "Next steps:"
    echo "1. Setup conda environment: ./setup_nextflow_integration.sh"
    echo "2. Set API key: echo 'GEMINI_API_KEY=your_key' >> .env"
    echo "3. Run with agents: ./run_nextflow_integrated.sh 'Your query'"
else
    echo ""
    echo "❌ Stub workflow failed"
    echo "Check the error messages above"
fi
