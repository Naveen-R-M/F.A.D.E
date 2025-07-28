#!/usr/bin/env python3
"""
Job status monitoring script for F.A.D.E
Allows checking the status of running jobs
"""

import os
import sys
import argparse
import glob
import json
import time
import datetime
import psutil

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="F.A.D.E Job Status Monitor")
    
    parser.add_argument(
        "--job-dir", "-d",
        type=str,
        help="Directory containing the job information (default: search all results_* directories)"
    )
    
    parser.add_argument(
        "--job-id", "-j",
        type=int,
        help="Process ID of the job to monitor"
    )
    
    parser.add_argument(
        "--watch", "-w",
        action="store_true",
        help="Watch the job status continuously (refreshes every 5 seconds)"
    )
    
    return parser.parse_args()

def find_job_dirs():
    """Find all job directories."""
    return glob.glob("results_*")

def get_job_status(job_dir):
    """Get status of a job from its directory."""
    job_info_file = os.path.join(job_dir, "job_info.txt")
    log_file = os.path.join(job_dir, "fade.log")
    results_file = os.path.join(job_dir, "results.json")
    
    status = {
        "directory": job_dir,
        "status": "Unknown",
        "pid": None,
        "started_at": None,
        "query": None,
        "log_file": log_file if os.path.exists(log_file) else None,
        "results_file": results_file if os.path.exists(results_file) else None
    }
    
    # Check if job info file exists
    if not os.path.exists(job_info_file):
        status["status"] = "No job info found"
        return status
    
    # Parse job info file
    with open(job_info_file, "r") as f:
        for line in f:
            if line.startswith("Job ID:"):
                status["pid"] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Started at:"):
                status["started_at"] = line.split(":", 1)[1].strip()
            elif line.startswith("Query:"):
                status["query"] = line.split(":", 1)[1].strip()
    
    # Check if process is running
    if status["pid"]:
        try:
            process = psutil.Process(status["pid"])
            if process.is_running():
                status["status"] = "Running"
                status["runtime"] = str(datetime.timedelta(seconds=int(time.time() - process.create_time())))
            else:
                status["status"] = "Not running"
        except psutil.NoSuchProcess:
            status["status"] = "Not running"
    
    # Check if results file exists
    if os.path.exists(results_file):
        status["status"] = "Completed"
        
        # Get file size
        status["results_size"] = f"{os.path.getsize(results_file) / 1024:.1f} KB"
        
        # Check results content
        try:
            with open(results_file, "r") as f:
                results = json.load(f)
                
            # Check if structure predictor results exist
            if "structure_predictor" in results:
                status["has_structures"] = True
                status["structure_count"] = len(results["structure_predictor"].get("structures", {}))
            else:
                status["has_structures"] = False
        except Exception as e:
            status["results_error"] = str(e)
    
    # Get the last few lines of the log file
    if os.path.exists(log_file):
        try:
            with open(log_file, "r") as f:
                log_lines = f.readlines()
                status["last_log_lines"] = log_lines[-5:] if log_lines else []
        except Exception as e:
            status["log_error"] = str(e)
    
    return status

def display_job_status(status):
    """Display job status in a readable format."""
    print(f"\nJob in directory: {status['directory']}")
    print(f"Status: {status['status']}")
    print(f"Process ID: {status['pid']}")
    print(f"Started at: {status['started_at']}")
    print(f"Query: {status['query']}")
    
    if "runtime" in status:
        print(f"Runtime: {status['runtime']}")
    
    if status["log_file"]:
        print(f"Log file: {status['log_file']}")
        
        if "last_log_lines" in status and status["last_log_lines"]:
            print("\nLast few log entries:")
            for line in status["last_log_lines"]:
                print(f"  {line.strip()}")
    
    if status["results_file"]:
        print(f"\nResults file: {status['results_file']}")
        
        if "results_size" in status:
            print(f"Results size: {status['results_size']}")
        
        if "has_structures" in status:
            if status["has_structures"]:
                print(f"Structure count: {status['structure_count']}")
            else:
                print("No structures in results (likely skipped structure prediction)")
    
    print("\n" + "-" * 50)

def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Handle job directory or job ID
    job_dirs = []
    
    if args.job_dir:
        if os.path.exists(args.job_dir):
            job_dirs = [args.job_dir]
        else:
            print(f"Error: Job directory '{args.job_dir}' not found.")
            sys.exit(1)
    elif args.job_id:
        # Search for job directory with matching PID
        for job_dir in find_job_dirs():
            job_info_file = os.path.join(job_dir, "job_info.txt")
            if os.path.exists(job_info_file):
                with open(job_info_file, "r") as f:
                    for line in f:
                        if line.startswith("Job ID:") and int(line.split(":", 1)[1].strip()) == args.job_id:
                            job_dirs = [job_dir]
                            break
        
        if not job_dirs:
            print(f"Error: No job directory found for PID {args.job_id}.")
            sys.exit(1)
    else:
        # Find all job directories
        job_dirs = find_job_dirs()
        
        if not job_dirs:
            print("No job directories found.")
            sys.exit(0)
    
    # Watch mode
    if args.watch:
        try:
            while True:
                os.system('cls' if os.name == 'nt' else 'clear')
                print(f"F.A.D.E Job Status Monitor - {datetime.datetime.now()}")
                print(f"Monitoring {len(job_dirs)} job(s)")
                
                for job_dir in job_dirs:
                    status = get_job_status(job_dir)
                    display_job_status(status)
                
                print("\nPress Ctrl+C to exit")
                time.sleep(5)
        except KeyboardInterrupt:
            print("\nMonitoring stopped.")
            sys.exit(0)
    else:
        # One-time status check
        for job_dir in job_dirs:
            status = get_job_status(job_dir)
            display_job_status(status)

if __name__ == "__main__":
    main()
