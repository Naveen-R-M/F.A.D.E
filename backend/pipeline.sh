#!/bin/bash

# F.A.D.E Simple Pipeline - Direct execution
# Single script implementation of F.A.D.E workflow

set -e  # Exit on any error

# Configuration
QUERY="${1:-Find molecules targeting EGFR for lung cancer}"
OUTPUT_DIR="${2:-results_$(date +%Y%m%d_%H%M%S)}"

echo "================================================================"
echo "F.A.D.E PIPELINE"
echo "================================================================"
echo "Query: $QUERY"
echo "Output: $OUTPUT_DIR"
echo "Time: $(date)"
echo "================================================================"

# Check API key
if [[ -z "${GEMINI_API_KEY:-}" ]]; then
    echo "ERROR: GEMINI_API_KEY not set"
    exit 1
fi

# Get script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script directory: $SCRIPT_DIR"

# Setup output directory
mkdir -p "$OUTPUT_DIR/target_selection"
mkdir -p "$OUTPUT_DIR/logs"

# Set Python path
export PYTHONPATH="$SCRIPT_DIR:${PYTHONPATH:-}"

echo "Starting RCSB Target Selection..."

# Try the simple selector first
SIMPLE_SCRIPT="$SCRIPT_DIR/nextflow/bin/run_simple_rcsb_selector.py"
FULL_SCRIPT="$SCRIPT_DIR/nextflow/bin/run_rcsb_target_selector.py"

if [[ -f "$SIMPLE_SCRIPT" ]]; then
    SCRIPT_TO_USE="$SIMPLE_SCRIPT"
    echo "Using simple RCSB selector: $SIMPLE_SCRIPT"
elif [[ -f "$FULL_SCRIPT" ]]; then
    SCRIPT_TO_USE="$FULL_SCRIPT"
    echo "Using full RCSB selector: $FULL_SCRIPT"
else
    echo "ERROR: No target selector script found"
    echo "Looked for:"
    echo "  $SIMPLE_SCRIPT"
    echo "  $FULL_SCRIPT"
    exit 1
fi

# Change to output directory
cd "$OUTPUT_DIR/target_selection"

# Run target selection
echo "Executing target selection..."
python "$SCRIPT_TO_USE" \
    --query "$QUERY" \
    --output-dir . \
    --api-key "$GEMINI_API_KEY" \
    --model "models/gemini-2.5-flash" \
    --verbose

# Check outputs
echo "Checking outputs..."

if [[ -f "target_info.json" ]]; then
    echo "✓ target_info.json created"
    echo "Target info:"
    cat target_info.json | python3 -m json.tool 2>/dev/null || cat target_info.json
else
    echo "✗ target_info.json missing"
    exit 1
fi

if [[ -f "protein.fasta" ]]; then
    echo "✓ protein.fasta created"
    echo "FASTA preview:"
    head -2 protein.fasta
else
    echo "✗ protein.fasta missing"
    exit 1
fi

if [[ -f "structure.pdb" ]]; then
    echo "✓ structure.pdb created"
    echo "PDB info:"
    head -5 structure.pdb | grep -E "HEADER|TITLE|COMPND" || head -2 structure.pdb
else
    echo "✗ structure.pdb missing"
    exit 1
fi

if [[ -f "requirements.json" ]]; then
    echo "✓ requirements.json created"
else
    echo "✗ requirements.json missing"
    exit 1
fi

echo "================================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY!"
echo "Results location: $OUTPUT_DIR"
echo "================================================================"

# Show summary
echo ""
echo "SUMMARY:"
echo "--------"
echo "Query processed: $QUERY"
echo "Results saved to: $OUTPUT_DIR/target_selection/"
echo ""
echo "Key files:"
echo "- target_info.json (target metadata)"
echo "- structure.pdb (protein structure)"
echo "- protein.fasta (protein sequence)"
echo "- requirements.json (drug requirements)"
echo ""
echo "Next steps: Add molecule generation and docking steps"