#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load .env file
if [ -f "${SCRIPT_DIR}/.env" ]; then
    # shellcheck disable=SC1090
    source "${SCRIPT_DIR}/.env"
    echo "✅ Loaded .env file"
    echo " API key length: ${#GEMINI_API_KEY}"
else
    echo "❌ .env file not found: ${SCRIPT_DIR}/.env"
    exit 1
fi

# Validate API key
if [ -z "${GEMINI_API_KEY:-}" ]; then
    echo "❌ GEMINI_API_KEY not found in .env file"
    exit 1
fi

# Get query from command line
QUERY="${1:-}"
if [ -z "$QUERY" ]; then
    echo "Usage: $0 'Your drug discovery query'"
    echo "Example: $0 'Find molecules targeting EGFR for lung cancer'"
    exit 1
fi

# Setup environment
echo "🔧 Setting up environment..."
unset JAVA_HOME JAVA_PATH JDK_HOME || true

if command -v module >/dev/null 2>&1; then
    module purge
    module load miniconda3/24.11.1
    module load nextflow/24.10.3
else
    echo "⚠️ module command not found, skipping module load"
fi

# Activate conda environment (fall back safely if needed)
CONDA_ENV_PATH="${SCRATCH:-${SCRIPT_DIR}}/conda-envs/fade"
if command -v conda >/dev/null 2>&1; then
    # Ensure conda functions are available
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV_PATH}" || source activate "${CONDA_ENV_PATH}" || echo "⚠️ Failed to activate conda env"
else
    if [ -f "${CONDA_ENV_PATH}/bin/activate" ]; then
        # shellcheck disable=SC1090
        source "${CONDA_ENV_PATH}/bin/activate"
    else
        echo "⚠️ Could not activate conda environment at ${CONDA_ENV_PATH}"
    fi
fi

# Prepare output directory
OUTPUT_DIR="${SCRIPT_DIR}/results_$(date +'%Y%m%d_%H%M%S')"
mkdir -p "$OUTPUT_DIR"

echo "🚀 Running F.A.D.E Nextflow Pipeline..."
echo " Query: $QUERY"
echo " Output: $OUTPUT_DIR"
echo " API Key: ${GEMINI_API_KEY:0:10}..."

# Change to nextflow directory
if [ -d "${SCRIPT_DIR}/nextflow" ]; then
    cd "${SCRIPT_DIR}/nextflow"
else
    echo "❌ nextflow directory not found: ${SCRIPT_DIR}/nextflow"
    exit 1
fi

# Run Nextflow with explicit API key parameter
nextflow run main.nf \
    --query "$QUERY" \
    --output_dir "$OUTPUT_DIR" \
    --max_molecules 10 \
    --gemini_api_key "$GEMINI_API_KEY" \
    --fade_env_path "${CONDA_ENV_PATH}" \
    -work-dir "${SCRIPT_DIR}/work_api_test"

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "🎉 Pipeline completed successfully!"
    echo "📂 Results in: $OUTPUT_DIR"

    # Show quick summary if available
    if [ -f "$OUTPUT_DIR/07_final_report/analysis_summary.md" ]; then
        echo ""
        echo "📋 Quick Summary:"
        head -n 20 "$OUTPUT_DIR/07_final_report/analysis_summary.md"
    fi
else
    echo ""
    echo "❌ Pipeline failed. Check logs for details."
    echo "📂 Work directory: ${SCRIPT_DIR}/work_api_test"
    exit $EXIT_CODE
fi
```# filepath: untitled:Untitled-1
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load .env file
if [ -f "${SCRIPT_DIR}/.env" ]; then
    # shellcheck disable=SC1090
    source "${SCRIPT_DIR}/.env"
    echo "✅ Loaded .env file"
    echo " API key length: ${#GEMINI_API_KEY}"
else
    echo "❌ .env file not found: ${SCRIPT_DIR}/.env"
    exit 1
fi

# Validate API key
if [ -z "${GEMINI_API_KEY:-}" ]; then
    echo "❌ GEMINI_API_KEY not found in .env file"
    exit 1
fi

# Get query from command line
QUERY="${1:-}"
if [ -z "$QUERY" ]; then
    echo "Usage: $0 'Your drug discovery query'"
    echo "Example: $0 'Find molecules targeting EGFR for lung cancer'"
    exit 1
fi

# Setup environment
echo "🔧 Setting up environment..."
unset JAVA_HOME JAVA_PATH JDK_HOME || true

if command -v module >/dev/null 2>&1; then
    module purge
    module load miniconda3/24.11.1
    module load nextflow/24.10.3
else
    echo "⚠️ module command not found, skipping module load"
fi

# Activate conda environment (fall back safely if needed)
CONDA_ENV_PATH="${SCRATCH:-${SCRIPT_DIR}}/conda-envs/fade"
if command -v conda >/dev/null 2>&1; then
    # Ensure conda functions are available
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV_PATH}" || source activate "${CONDA_ENV_PATH}" || echo "⚠️ Failed to activate conda env"
else
    if [ -f "${CONDA_ENV_PATH}/bin/activate" ]; then
        # shellcheck disable=SC1090
        source "${CONDA_ENV_PATH}/bin/activate"
    else
        echo "⚠️ Could not activate conda environment at ${CONDA_ENV_PATH}"
    fi
fi

# Prepare output directory
OUTPUT_DIR="${SCRIPT_DIR}/results_$(date +'%Y%m%d_%H%M%S')"
mkdir -p "$OUTPUT_DIR"

echo "🚀 Running F.A.D.E Nextflow Pipeline..."
echo " Query: $QUERY"
echo " Output: $OUTPUT_DIR"
echo " API Key: ${GEMINI_API_KEY:0:10}..."

# Change to nextflow directory
if [ -d "${SCRIPT_DIR}/nextflow" ]; then
    cd "${SCRIPT_DIR}/nextflow"
else
    echo "❌ nextflow directory not found: ${SCRIPT_DIR}/nextflow"
    exit 1
fi

# Run Nextflow with explicit API key parameter
nextflow run main.nf \
    --query "$QUERY" \
    --output_dir "$OUTPUT_DIR" \
    --max_molecules 10 \
    --gemini_api_key "$GEMINI_API_KEY" \
    --fade_env_path "${CONDA_ENV_PATH}" \
    -work-dir "${SCRIPT_DIR}/work_api_test"

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "🎉 Pipeline completed successfully!"
    echo "📂 Results in: $OUTPUT_DIR"

    # Show quick summary if available
    if [ -f "$OUTPUT_DIR/07_final_report/analysis_summary.md" ]; then
        echo ""
        echo "📋 Quick Summary:"
        head -n 20 "$OUTPUT_DIR/07_final_report/analysis_summary.md"
    fi
else
    echo ""
    echo "❌ Pipeline failed. Check logs for details."
    echo "📂 Work directory: ${SCRIPT_DIR}/work_api_test"
