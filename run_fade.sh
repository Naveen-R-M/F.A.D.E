#!/bin/bash
# Simple wrapper script for running F.A.D.E jobs on HPC clusters

# Check if query parameter is provided
if [ -z "$1" ]; then
    echo "Usage: $0 \"your query here\" [--skip-structure-prediction]"
    echo "Example: $0 \"Find molecules targeting KRAS G12D with good BBB permeability\""
    exit 1
fi

# Create timestamped directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="results_${TIMESTAMP}"

# Check if --skip-structure-prediction flag is provided
SKIP_FLAG=""
if [ "$2" == "--skip-structure-prediction" ]; then
    SKIP_FLAG="--skip-structure-prediction"
    echo "Structure prediction will be skipped"
fi

# Create output directory
mkdir -p $OUTPUT_DIR

# Create log file
LOG_FILE="${OUTPUT_DIR}/fade.log"
echo "F.A.D.E Job started at $(date)" > $LOG_FILE
echo "Query: $1" >> $LOG_FILE
echo "Output directory: $OUTPUT_DIR" >> $LOG_FILE
echo "Skip structure prediction: ${SKIP_FLAG:-false}" >> $LOG_FILE
echo "" >> $LOG_FILE

# Function to check if the job is running
function is_job_running() {
    if [ -e "/proc/$1" ]; then
        return 0  # Job is running
    else
        return 1  # Job is not running
    fi
}

# Run F.A.D.E in the background
nohup python main.py --query "$1" --output-dir "$OUTPUT_DIR" $SKIP_FLAG --foreground >> $LOG_FILE 2>&1 &
PID=$!

# Save job info
JOB_INFO="${OUTPUT_DIR}/job_info.txt"
echo "Job ID: $PID" > $JOB_INFO
echo "Started at: $(date)" >> $JOB_INFO
echo "Query: $1" >> $JOB_INFO
echo "Output directory: $OUTPUT_DIR" >> $JOB_INFO
echo "Log file: $LOG_FILE" >> $JOB_INFO
echo "Skip structure prediction: ${SKIP_FLAG:-false}" >> $JOB_INFO

# Create status checker script
STATUS_SCRIPT="${OUTPUT_DIR}/check_status.sh"
cat > $STATUS_SCRIPT << 'EOF'
#!/bin/bash
# Simple status checker for F.A.D.E job

# Get job info
JOB_DIR=$(dirname "$0")
JOB_INFO="${JOB_DIR}/job_info.txt"
LOG_FILE="${JOB_DIR}/fade.log"
RESULTS_FILE="${JOB_DIR}/results.json"

# Extract PID from job info
PID=$(grep "Job ID:" "$JOB_INFO" | cut -d ":" -f2 | tr -d ' ')

# Check if job is running
if [ -e "/proc/$PID" ]; then
    echo "Status: RUNNING (PID: $PID)"
    
    # Get runtime
    START_TIME=$(grep "Started at:" "$JOB_INFO" | cut -d ":" -f2- | sed 's/^[ \t]*//')
    START_EPOCH=$(date -d "$START_TIME" +%s)
    CURRENT_EPOCH=$(date +%s)
    RUNTIME=$((CURRENT_EPOCH - START_EPOCH))
    
    # Format runtime as HH:MM:SS
    HOURS=$((RUNTIME / 3600))
    MINUTES=$(( (RUNTIME % 3600) / 60 ))
    SECONDS=$((RUNTIME % 60))
    
    echo "Runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
    
    echo ""
    echo "Last few log lines:"
    tail -n 10 "$LOG_FILE"
else
    if [ -f "$RESULTS_FILE" ]; then
        echo "Status: COMPLETED"
        echo "Results file: $RESULTS_FILE"
        echo "Results size: $(du -h "$RESULTS_FILE" | cut -f1)"
    else
        echo "Status: NOT RUNNING (possibly failed)"
    fi
    
    echo ""
    echo "Last few log lines:"
    tail -n 10 "$LOG_FILE"
fi
EOF

# Make the script executable
chmod +x $STATUS_SCRIPT

# Print status message
echo "F.A.D.E job started in background mode:"
echo "  Job ID: $PID"
echo "  Output directory: $OUTPUT_DIR"
echo "  Log file: $LOG_FILE"
echo ""
echo "To check status:"
echo "  1. Run the status script: $STATUS_SCRIPT"
echo "  2. View log file: tail -f $LOG_FILE"
echo "  3. Use the job monitor: python monitor_jobs.py $OUTPUT_DIR"
echo ""
echo "Your job is now running in the background. You can close this terminal if needed."

# Add a watchdog to detect job completion
(
    # Keep checking if job is running
    while is_job_running $PID; do
        sleep 60  # Check every minute
    done
    
    # Check if results file exists (job completed successfully)
    if [ -f "${OUTPUT_DIR}/results.json" ]; then
        echo "F.A.D.E job completed successfully at $(date)" >> $LOG_FILE
        echo "Results saved to: ${OUTPUT_DIR}/results.json" >> $LOG_FILE
    else
        echo "F.A.D.E job stopped without producing results at $(date)" >> $LOG_FILE
    fi
) &

exit 0
