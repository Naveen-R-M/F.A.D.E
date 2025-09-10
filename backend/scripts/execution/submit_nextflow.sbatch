#!/bin/bash
#SBATCH --job-name=fade-nextflow
#SBATCH --partition=gpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --output=fade-nextflow-%j.log
#SBATCH --error=fade-nextflow-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rajagopalmohanraj.n@northeastern.edu

#############################################
# F.A.D.E Nextflow Pipeline SLURM Submission
#############################################

# Load required modules
module load nextflow/24.10.3
module load python/3.13.5

# Set environment variables
export NXF_WORK=$SCRATCH/.nextflow/work
export NXF_TEMP=$SCRATCH/.nextflow/temp
export SINGULARITY_CACHEDIR=$SCRATCH/singularity_cache

# Load API keys from .env file if it exists
if [ -f "$SCRATCH/F.A.D.E/.env" ]; then
    export $(grep -v '^#' $SCRATCH/F.A.D.E/.env | xargs)
fi

# Parse arguments (passed from sbatch command)
QUERY="${1:-Design a selective EGFR inhibitor for lung cancer}"
OUTPUT_DIR="${2:-results_$(date +%Y%m%d_%H%M%S)}"
MAX_MOLECULES="${3:-1000}"

# Create necessary directories
mkdir -p $NXF_WORK $NXF_TEMP $SINGULARITY_CACHEDIR

# Change to nextflow directory
cd $SCRATCH/F.A.D.E/nextflow

# Log start information
echo "=========================================="
echo "F.A.D.E Nextflow Pipeline SLURM Job"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "Query: $QUERY"
echo "Output Directory: $OUTPUT_DIR"
echo "Max Molecules: $MAX_MOLECULES"
echo "=========================================="

# Run the Nextflow pipeline
nextflow run main.nf \
    -profile northeastern \
    --query "$QUERY" \
    --output_dir "$OUTPUT_DIR" \
    --max_molecules "$MAX_MOLECULES" \
    --max_iterations 3 \
    --docking_method vina \
    -with-report "$OUTPUT_DIR/pipeline_report.html" \
    -with-timeline "$OUTPUT_DIR/timeline.html" \
    -with-trace "$OUTPUT_DIR/trace.txt" \
    -ansi-log false

# Check exit status
EXIT_CODE=$?

echo "=========================================="
echo "End Time: $(date)"
echo "Exit Code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully!"
    echo "Results saved in: $OUTPUT_DIR"
    
    # Generate summary email content
    if [ -f "$OUTPUT_DIR/07_report/summary.txt" ]; then
        echo ""
        echo "PIPELINE SUMMARY:"
        cat "$OUTPUT_DIR/07_report/summary.txt"
    fi
else
    echo "Pipeline failed with exit code $EXIT_CODE"
    echo "Check the log files for details"
fi

echo "=========================================="

exit $EXIT_CODE
