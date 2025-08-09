#!/bin/bash

# F.A.D.E Integrated Nextflow Runner
# This script sets up the environment and runs the Nextflow workflow with actual Python agents

set -e  # Exit on any error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEXTFLOW_DIR="${SCRIPT_DIR}/nextflow"
CONDA_ENV_PATH="${SCRATCH}/conda-envs/fade"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')] $1${NC}"
}

error() {
    echo -e "${RED}[ERROR] $1${NC}" >&2
}

warning() {
    echo -e "${YELLOW}[WARNING] $1${NC}"
}

success() {
    echo -e "${GREEN}[SUCCESS] $1${NC}"
}

# Function to setup environment properly
setup_environment() {
    log "Setting up environment for Nextflow..."
    
    # Clear Java conflicts
    unset JAVA_HOME
    unset JAVA_PATH
    unset JDK_HOME
    
    # Purge modules and reload clean
    module purge
    module load miniconda3/24.11.1
    module load nextflow/24.10.3
    
    # Load .env file if it exists
    if [ -f "${SCRIPT_DIR}/.env" ]; then
        source "${SCRIPT_DIR}/.env"
        log "Loaded environment variables from .env"
    fi
    
    # Verify Nextflow is available
    if ! command -v nextflow &> /dev/null; then
        error "Nextflow not available after module loading"
        exit 1
    fi
    
    # Verify API key
    if [ -z "$GEMINI_API_KEY" ]; then
        warning "GEMINI_API_KEY not set. Some agents may fail."
    else
        success "GEMINI_API_KEY is configured"
    fi
    
    # Set Python path
    export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
    log "Python path set to include: $SCRIPT_DIR"
    
    success "Environment setup complete"
}

# Function to check if conda environment exists
check_conda_env() {
    if [ ! -d "$CONDA_ENV_PATH" ]; then
        warning "Conda environment not found at: $CONDA_ENV_PATH"
        warning "Run ./setup_nextflow_integration.sh to create it"
        warning "For now, using system Python..."
    else
        success "Conda environment found: $CONDA_ENV_PATH"
    fi
}

# Function to print usage
usage() {
    cat << 'USAGE_EOF'
Usage: ./run_nextflow_integrated.sh [OPTIONS] "QUERY"

F.A.D.E Integrated Nextflow Pipeline Runner

ARGUMENTS:
    QUERY                   Natural language drug discovery query (required)

OPTIONS:
    -h, --help             Show this help message
    -o, --output-dir DIR   Output directory (default: auto-generated)
    -m, --max-molecules N  Maximum molecules to generate (default: 100)
    -i, --iterations N     Maximum optimization iterations (default: 1)
    -d, --docking METHOD   Docking method: vina|glide (default: vina)
    -p, --profile PROFILE  Nextflow profile: local|northeastern|debug|stub (default: northeastern)
    -r, --resume           Resume previous run
    -w, --work-dir DIR     Nextflow work directory (default: ./work)
    --dry-run              Run with stub processes for testing
    --trace                Enable detailed execution tracing
    --debug                Enable debug mode

EXAMPLES:
    # Basic usage
    ./run_nextflow_integrated.sh "Find molecules targeting KRAS G12D with good BBB permeability"
    
    # With custom parameters
    ./run_nextflow_integrated.sh -m 200 -i 3 -d glide "Design EGFR inhibitors for lung cancer"
    
    # Debug mode
    ./run_nextflow_integrated.sh --debug --trace "Test query for BRAF V600E"
    
    # Dry run (stub testing)
    ./run_nextflow_integrated.sh --dry-run "Test workflow with stubs"
USAGE_EOF
}

# Parse command line arguments
QUERY=""
OUTPUT_DIR=""
MAX_MOLECULES=100
ITERATIONS=1
DOCKING_METHOD="vina"
PROFILE="northeastern"
RESUME=""
WORK_DIR=""
DRY_RUN=false
TRACE=false
DEBUG=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -m|--max-molecules)
            MAX_MOLECULES="$2"
            shift 2
            ;;
        -i|--iterations)
            ITERATIONS="$2"
            shift 2
            ;;
        -d|--docking)
            DOCKING_METHOD="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -r|--resume)
            RESUME="-resume"
            shift
            ;;
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            PROFILE="stub"
            shift
            ;;
        --trace)
            TRACE=true
            shift
            ;;
        --debug)
            DEBUG=true
            shift
            ;;
        -*)
            error "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [ -z "$QUERY" ]; then
                QUERY="$1"
            else
                error "Multiple queries provided. Please provide only one query."
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate required arguments
if [ -z "$QUERY" ]; then
    error "Query is required"
    usage
    exit 1
fi

# Main execution
main() {
    log "Starting F.A.D.E Integrated Nextflow Pipeline"
    log "Query: $QUERY"
    
    # Setup environment (this includes Java conflict resolution)
    setup_environment
    check_conda_env
    
    # Set default output directory if not provided
    if [ -z "$OUTPUT_DIR" ]; then
        OUTPUT_DIR="${SCRIPT_DIR}/results_$(date +'%Y%m%d_%H%M%S')"
    fi
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    log "Output directory: $OUTPUT_DIR"
    
    # Set work directory
    if [ -z "$WORK_DIR" ]; then
        WORK_DIR="${SCRIPT_DIR}/work"
    fi
    
    # Build Nextflow command
    NEXTFLOW_CMD="nextflow run ${NEXTFLOW_DIR}/main.nf"
    NEXTFLOW_CMD="$NEXTFLOW_CMD --query '$QUERY'"
    NEXTFLOW_CMD="$NEXTFLOW_CMD --output_dir '$OUTPUT_DIR'"
    NEXTFLOW_CMD="$NEXTFLOW_CMD --max_molecules $MAX_MOLECULES"
    NEXTFLOW_CMD="$NEXTFLOW_CMD --max_iterations $ITERATIONS"
    NEXTFLOW_CMD="$NEXTFLOW_CMD --docking_method $DOCKING_METHOD"
    NEXTFLOW_CMD="$NEXTFLOW_CMD --fade_env_path '$CONDA_ENV_PATH'"
    NEXTFLOW_CMD="$NEXTFLOW_CMD -profile $PROFILE"
    NEXTFLOW_CMD="$NEXTFLOW_CMD -work-dir '$WORK_DIR'"
    
    # Add optional parameters
    if [ "$RESUME" = "-resume" ]; then
        NEXTFLOW_CMD="$NEXTFLOW_CMD -resume"
    fi
    
    if [ "$TRACE" = true ]; then
        NEXTFLOW_CMD="$NEXTFLOW_CMD -with-trace"
        NEXTFLOW_CMD="$NEXTFLOW_CMD --trace_enabled true"
    fi
    
    if [ "$DEBUG" = true ]; then
        NEXTFLOW_CMD="$NEXTFLOW_CMD -with-dag ${OUTPUT_DIR}/workflow_dag.svg"
        NEXTFLOW_CMD="$NEXTFLOW_CMD -with-timeline ${OUTPUT_DIR}/timeline.html"
        NEXTFLOW_CMD="$NEXTFLOW_CMD -with-report ${OUTPUT_DIR}/execution_report.html"
        NEXTFLOW_CMD="$NEXTFLOW_CMD --debug true"
    fi
    
    # Change to nextflow directory
    cd "$NEXTFLOW_DIR"
    
    # Log the command
    log "Running Nextflow workflow..."
    log "Command: nextflow run main.nf ..."
    
    # Execute workflow
    if [ "$DRY_RUN" = true ]; then
        warning "Running in DRY RUN mode (stub processes only)"
    fi
    
    log "Executing: $NEXTFLOW_CMD"
    eval "$NEXTFLOW_CMD" || {
        error "Nextflow pipeline failed"
        error "Check logs in: $OUTPUT_DIR"
        error "Work directory: $WORK_DIR"
        
        # Show recent error logs if available
        if [ -f ".nextflow.log" ]; then
            error "Recent Nextflow log:"
            tail -20 .nextflow.log
        fi
        
        exit 1
    }
    
    success "Pipeline completed successfully!"
    success "Results available in: $OUTPUT_DIR"
    
    # Display quick summary if final report exists
    if [ -f "$OUTPUT_DIR/07_final_report/analysis_summary.md" ]; then
        log "=== QUICK SUMMARY ==="
        head -20 "$OUTPUT_DIR/07_final_report/analysis_summary.md"
        log "=== END SUMMARY ==="
        log "Full report: $OUTPUT_DIR/07_final_report/analysis_summary.md"
    fi
}

# Run main function
main "$@"
