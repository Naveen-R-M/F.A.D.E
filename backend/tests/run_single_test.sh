#!/bin/bash
# Script to run a single test with F.A.D.E

# Set the project directory
PROJECT_DIR="$SCRATCH/F.A.D.E"
cd $PROJECT_DIR

# Load miniconda module
module load miniconda3/24.11.1

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $SCRATCH/conda-envs/fade

# Create output directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="$PROJECT_DIR/test_results/single_test_$TIMESTAMP"
mkdir -p $OUTPUT_DIR

echo "==================================="
echo "F.A.D.E Single Test"
echo "Running query: $1"
echo "Results will be saved to: $OUTPUT_DIR"
echo "==================================="

# Run the query
python main.py --query "$1" --output-dir "$OUTPUT_DIR" --log-level INFO

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "==================================="
    echo "Test completed successfully!"
    echo "Results saved to: $OUTPUT_DIR"
    echo "Log files are in: $OUTPUT_DIR/logs"
    echo "==================================="
else
    echo "==================================="
    echo "Test failed!"
    echo "Check logs for errors: $OUTPUT_DIR/logs"
    echo "==================================="
fi
