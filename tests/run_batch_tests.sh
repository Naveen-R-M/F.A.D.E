#!/bin/bash
# Script to test F.A.D.E with different test queries

# Set the project directory
PROJECT_DIR="$SCRATCH/F.A.D.E"
cd $PROJECT_DIR

# Load miniconda module
module load miniconda3/24.11.1

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $SCRATCH/conda-envs/fade

# Create results directory with timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RESULTS_DIR="$PROJECT_DIR/test_results/$TIMESTAMP"
mkdir -p $RESULTS_DIR

echo "==================================="
echo "F.A.D.E Test Queries"
echo "Results will be saved to: $RESULTS_DIR"
echo "==================================="

# Log file for all test results
LOG_FILE="$RESULTS_DIR/test_queries.log"

# Function to run a test query
run_test_query() {
    local query_name="$1"
    local query="$2"
    local output_dir="$RESULTS_DIR/$query_name"
    
    echo "Running test: $query_name"
    echo "Query: $query"
    
    mkdir -p "$output_dir"
    
    # Log start time
    start_time=$(date +%s)
    
    # Run the query and capture output
    python main.py --query "$query" --output-dir "$output_dir" --log-level INFO > "$output_dir/stdout.log" 2> "$output_dir/stderr.log"
    exit_code=$?
    
    # Log end time and calculate duration
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    # Check if the command was successful
    if [ $exit_code -eq 0 ]; then
        status="SUCCESS"
    else
        status="FAILED"
    fi
    
    # Log the result
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $query_name: $status (${duration}s)" | tee -a "$LOG_FILE"
    echo "Output directory: $output_dir" | tee -a "$LOG_FILE"
    echo "Log files: $output_dir/logs" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Return the exit code
    return $exit_code
}

# Test queries
echo "Starting test queries..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Test 1: Original query
run_test_query "kras_g12d_basic" "Find molecules targeting KRAS G12D with good BBB permeability"

# Test 2: Detailed query
run_test_query "kras_detailed" "I'm looking for potential drug candidates that could target the KRAS G12D mutant protein, which is implicated in pancreatic cancer. The molecules should be able to cross the blood-brain barrier, have low toxicity, and follow Lipinski's rule of five for drug-likeness."

# Test 3: Different target with mutation
run_test_query "braf_v600e" "Find molecules that inhibit BRAF V600E with good BBB permeability and low toxicity."

# Test 4: Dual targeting
run_test_query "dual_target" "Develop dual inhibitors for both EGFR and HER2 that can cross the blood-brain barrier."

# Test 5: ACE2 target
run_test_query "ace2_target" "Design compounds that block human ACE2 receptor binding to viral spike proteins while maintaining BBB permeability."

# Summary
echo "" | tee -a "$LOG_FILE"
echo "==================================" | tee -a "$LOG_FILE"
echo "Test Queries Completed" | tee -a "$LOG_FILE"
echo "Results saved to: $RESULTS_DIR" | tee -a "$LOG_FILE"
echo "==================================" | tee -a "$LOG_FILE"

# Create a consolidated logs directory
mkdir -p "$RESULTS_DIR/consolidated_logs"
find "$RESULTS_DIR" -name "logs" -type d | while read log_dir; do
    query_name=$(basename $(dirname "$log_dir"))
    cp -r "$log_dir"/* "$RESULTS_DIR/consolidated_logs/"
    echo "Copied logs from $query_name to consolidated logs directory" | tee -a "$LOG_FILE"
done

echo "Consolidated logs available at: $RESULTS_DIR/consolidated_logs" | tee -a "$LOG_FILE"
