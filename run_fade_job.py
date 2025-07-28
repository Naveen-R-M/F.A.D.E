#!/usr/bin/env python3
"""
Job submission script for F.A.D.E
Allows running the pipeline in the background
"""

import os
import sys
import argparse
import subprocess
import datetime
import time

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="F.A.D.E Job Submission")
    
    parser.add_argument(
        "--query", "-q",
        type=str,
        required=True,
        help="Natural language query describing the drug discovery goal"
    )
    
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        help="Directory where results will be saved (default: results_TIMESTAMP)"
    )
    
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level"
    )
    
    parser.add_argument(
        "--skip-structure-prediction",
        action="store_true",
        help="Skip the structure prediction step"
    )
    
    return parser.parse_args()

def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Generate timestamp for unique job identification
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Set output directory if not provided
    output_dir = args.output_dir if args.output_dir else f"results_{timestamp}"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create log file path
    log_file = os.path.join(output_dir, "fade.log")
    
    # Build command
    cmd = [
        "python", "main.py",
        "--query", args.query,
        "--output-dir", output_dir,
        "--log-level", args.log_level
    ]
    
    if args.skip_structure_prediction:
        cmd.append("--skip-structure-prediction")
    
    # Print job information
    print(f"Starting F.A.D.E job:")
    print(f"  Query: {args.query}")
    print(f"  Output directory: {output_dir}")
    print(f"  Log file: {log_file}")
    print(f"  Skip structure prediction: {args.skip_structure_prediction}")
    
    # Start the process
    with open(log_file, "w") as f:
        # Write header to log file
        f.write(f"F.A.D.E Job started at {datetime.datetime.now()}\n")
        f.write(f"Query: {args.query}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Skip structure prediction: {args.skip_structure_prediction}\n\n")
        f.flush()
        
        # Start the process in background
        process = subprocess.Popen(
            cmd,
            stdout=f,
            stderr=subprocess.STDOUT,
            close_fds=True
        )
    
    # Create a job info file
    job_info_file = os.path.join(output_dir, "job_info.txt")
    with open(job_info_file, "w") as f:
        f.write(f"Job ID: {process.pid}\n")
        f.write(f"Started at: {datetime.datetime.now()}\n")
        f.write(f"Query: {args.query}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Log file: {log_file}\n")
        f.write(f"Skip structure prediction: {args.skip_structure_prediction}\n")
    
    print(f"Job started with PID {process.pid}")
    print(f"You can monitor the progress with: tail -f {log_file}")
    print(f"Results will be saved to: {output_dir}")

if __name__ == "__main__":
    main()
