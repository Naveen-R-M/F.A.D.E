#!/bin/bash

# F.A.D.E Nextflow Runner with Explicit API Key Handling

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load .env file
if [ -f "${SCRIPT_DIR}/.env" ]; then
    source "${SCRIPT_DIR}/.env"
    echo "✅ Loaded .env file"
    echo "   API key length: ${#GEMINI_API_KEY}"
else
    echo "❌ .env file not found"
    exit 1
fi

# Validate API key
if [ -z "$GEMINI_API_KEY" ]; then
    echo "❌ GEMINI_API_KEY not found in .env file"
    exit 1
fi

# Get query from command line
QUERY="$1"
if [ -z "$QUERY" ]; then
    echo "Usage: $0 'Your drug discovery query'"
    echo "Example: $0 'Find molecules targeting EGFR for lung cancer'"
    exit 1
fi

# Setup environment
echo "🔧 Setting up environment..."
unset JAVA_HOME JAVA_PATH JDK_HOME
module purge
module load miniconda3/24.11.1
module load nextflow/24.10.3

# Activate conda environment
source activate "${SCRATCH}/conda-envs/fade"

# Set output directory
OUTPUT_DIR="${SCRIPT_DIR}/results_$(date +'%Y%m%d_%H%M%S')"
mkdir -p "$OUTPUT_DIR"

echo "🚀 Running F.A.D.E Nextflow Pipeline..."
echo "   Query: $QUERY"
echo "   Output: $OUTPUT_DIR"
echo "   API Key: ${GEMINI_API_KEY:0:10}..."

# Change to nextflow directory
cd "${SCRIPT_DIR}/nextflow"

# Run Nextflow with explicit API key parameter
nextflow run main.nf \
    --query "$QUERY" \
    --output_dir "$OUTPUT_DIR" \
    --max_molecules 10 \
    --gemini_api_key "$GEMINI_API_KEY" \
    --fade_env_path "${SCRATCH}/conda-envs/fade" \
    -profile local \
    -work-dir "${SCRIPT_DIR}/work_api_test"

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Pipeline completed successfully!"
    echo "📂 Results in: $OUTPUT_DIR"
    
    # Show quick summary if available
    if [ -f "$OUTPUT_DIR/07_final_report/analysis_summary.md" ]; then
        echo ""
        echo "📋 Quick Summary:"
        head -20 "$OUTPUT_DIR/07_final_report/analysis_summary.md"
    fi
else
    echo ""
    echo "❌ Pipeline failed. Check logs for details."
    echo "📂 Work directory: ${SCRIPT_DIR}/work_api_test"
fi
