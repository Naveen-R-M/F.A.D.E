#!/bin/bash

#############################################
# F.A.D.E Nextflow Pipeline Runner
# Using Java 22 with clean environment
#############################################

# Set defaults
QUERY=""
PROFILE="northeastern"
RESUME=""
OUTPUT_DIR=""
MAX_MOLECULES=1000
MAX_ITERATIONS=3
DOCKING_METHOD="vina"
HELP=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -q|--query)
            QUERY="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --resume)
            RESUME="-resume"
            shift
            ;;
        --max-molecules)
            MAX_MOLECULES="$2"
            shift 2
            ;;
        --max-iterations)
            MAX_ITERATIONS="$2"
            shift 2
            ;;
        --docking-method)
            DOCKING_METHOD="$2"
            shift 2
            ;;
        -h|--help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Show help if requested
if [ "$HELP" = true ]; then
    cat << HELPTEXT
F.A.D.E Nextflow Pipeline Runner

Usage: ./run_nextflow.sh [OPTIONS]

Options:
    -q, --query <text>          Drug discovery query (required)
    -p, --profile <profile>     Execution profile (default: northeastern)
                               Options: northeastern, local, aws, gcp, dev
    -o, --output <dir>          Output directory (default: results_TIMESTAMP)
    --resume                    Resume from previous run
    --max-molecules <num>       Maximum molecules to generate (default: 1000)
    --max-iterations <num>      Maximum optimization iterations (default: 3)
    --docking-method <method>   Docking tool to use (default: vina)
                               Options: vina, glide, autodock
    -h, --help                 Show this help message

Examples:
    # Basic run
    ./run_nextflow.sh -q "Design EGFR inhibitor for lung cancer"
    
    # Resume previous run
    ./run_nextflow.sh -q "Design EGFR inhibitor" --resume
    
    # Custom output and molecules
    ./run_nextflow.sh -q "Find BTK inhibitors" -o btk_results --max-molecules 5000
    
    # Development mode with tracing
    ./run_nextflow.sh -q "Test query" -p dev

HELPTEXT
    exit 0
fi

# Check if query is provided
if [ -z "$QUERY" ]; then
    echo "Error: Query is required. Use -q or --query to specify."
    echo "Use -h or --help for usage information."
    exit 1
fi

echo "Setting up environment..."

# Clean environment completely
unset JAVA_HOME JDK_HOME JRE_HOME CLASSPATH LD_LIBRARY_PATH

# Use Java 22 cleanly
export JAVA_HOME=/shared/EL9/explorer/OpenJDK/22.0.2/jdk-22.0.2
export PATH=$JAVA_HOME/bin:$PATH

# Set clean library path for Java 22 only
export LD_LIBRARY_PATH=$JAVA_HOME/lib:$JAVA_HOME/lib/server

# Use Nextflow directly
NEXTFLOW=/shared/EL9/explorer/nextflow/24.10.3/bin/nextflow

# Activate conda environment for Python
echo "Activating conda environment..."
source /home/rajagopalmohanraj.n/.bashrc
conda activate fade

# Set Nextflow environment variables
export NXF_WORK=$SCRATCH/.nextflow/work
export NXF_TEMP=$SCRATCH/.nextflow/temp
export SINGULARITY_CACHEDIR=$SCRATCH/singularity_cache
export GEMINI_API_KEY=${GEMINI_API_KEY:-}
export NXF_JAVA_HOME=$JAVA_HOME

# Nextflow Java options
export NXF_OPTS="-Xms1g -Xmx4g"

# Load .env file if it exists
if [ -f "$SCRATCH/F.A.D.E/.env" ]; then
    export $(grep -v '^#' $SCRATCH/F.A.D.E/.env | xargs)
fi

# Check API key
if [ -z "$GEMINI_API_KEY" ]; then
    echo "Warning: GEMINI_API_KEY not set. Some features may not work."
    echo "Set it with: export GEMINI_API_KEY='your-key'"
fi

# Create necessary directories
mkdir -p $NXF_WORK $NXF_TEMP $SINGULARITY_CACHEDIR

# Change to nextflow directory
cd $SCRATCH/F.A.D.E/nextflow

# Build the Nextflow command
NF_CMD="$NEXTFLOW run main.nf"
NF_CMD="$NF_CMD -profile $PROFILE"
NF_CMD="$NF_CMD --query \"$QUERY\""
NF_CMD="$NF_CMD --max_molecules $MAX_MOLECULES"
NF_CMD="$NF_CMD --max_iterations $MAX_ITERATIONS"
NF_CMD="$NF_CMD --docking_method $DOCKING_METHOD"

if [ ! -z "$OUTPUT_DIR" ]; then
    NF_CMD="$NF_CMD --output_dir $OUTPUT_DIR"
fi

if [ ! -z "$RESUME" ]; then
    NF_CMD="$NF_CMD $RESUME"
fi

# Add reporting options for dev profile
if [ "$PROFILE" = "dev" ]; then
    NF_CMD="$NF_CMD -with-report report.html"
    NF_CMD="$NF_CMD -with-timeline timeline.html"
    NF_CMD="$NF_CMD -with-trace trace.txt"
    NF_CMD="$NF_CMD -with-dag pipeline_dag.pdf"
fi

# Print configuration
echo "=========================================="
echo "F.A.D.E Nextflow Pipeline Configuration"
echo "=========================================="
echo "Query: $QUERY"
echo "Profile: $PROFILE"
echo "Max Molecules: $MAX_MOLECULES"
echo "Max Iterations: $MAX_ITERATIONS"
echo "Docking Method: $DOCKING_METHOD"
echo "Output Dir: ${OUTPUT_DIR:-auto-generated}"
echo "Resume: ${RESUME:-no}"
echo "=========================================="
echo ""

# Debug: Verify environment
echo "Environment Check:"
echo "JAVA_HOME: $JAVA_HOME"
echo "Nextflow: $NEXTFLOW"
echo -n "Java: "
$JAVA_HOME/bin/java -version 2>&1 | head -1
echo -n "Nextflow: "
$NEXTFLOW -version 2>&1 | grep -i version | head -1 || echo "Nextflow version check failed"
echo ""

# Run the pipeline
echo "Starting Nextflow pipeline..."
echo "Command: $NF_CMD"
echo ""

eval $NF_CMD

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Pipeline completed successfully!"
    echo "=========================================="
    
    # Display output location
    if [ -z "$OUTPUT_DIR" ]; then
        OUTPUT_DIR=$(ls -td results_* 2>/dev/null | head -1)
    fi
    
    if [ -d "$OUTPUT_DIR" ]; then
        echo "Results saved in: $OUTPUT_DIR"
        echo ""
        echo "Key output files:"
        [ -f "$OUTPUT_DIR/07_report/final_report.html" ] && echo "  - HTML Report: $OUTPUT_DIR/07_report/final_report.html"
        [ -f "$OUTPUT_DIR/07_report/summary.txt" ] && echo "  - Summary: $OUTPUT_DIR/07_report/summary.txt"
        [ -f "$OUTPUT_DIR/05_docking/top_hits.csv" ] && echo "  - Top Hits: $OUTPUT_DIR/05_docking/top_hits.csv"
    fi
else
    echo ""
    echo "=========================================="
    echo "Pipeline failed. Check the logs for details."
    echo "=========================================="
    echo "To resume from the last checkpoint, run:"
    echo "./run_nextflow.sh -q \"$QUERY\" --resume"
    exit 1
fi
