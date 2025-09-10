#!/usr/bin/env python3
"""
Monitor F.A.D.E jobs and display status in real-time

This script is used to monitor the status of F.A.D.E jobs running in the background.
"""

import os
import sys
import time
import json
import glob
import argparse
import subprocess
from datetime import datetime

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="F.A.D.E Job Monitor")
    
    parser.add_argument(
        "job_dir",
        nargs="?",
        help="Job directory to monitor (if not provided, monitors all recent jobs)"
    )
    
    parser.add_argument(
        "--refresh",
        type=int,
        default=5,
        help="Refresh interval in seconds (default: 5)"
    )
    
    return parser.parse_args()

def find_recent_jobs(max_jobs=5):
    """Find the most recent job directories."""
    # Find all results_* directories
    job_dirs = glob.glob("results_*")
    
    # Sort by modification time (most recent first)
    job_dirs.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    
    # Return the most recent jobs
    return job_dirs[:max_jobs]

def get_job_info(job_dir):
    """Get information about a job."""
    job_info_file = os.path.join(job_dir, "job_info.txt")
    log_file = os.path.join(job_dir, "fade.log")
    results_file = os.path.join(job_dir, "results.json")
    
    job_info = {
        "directory": job_dir,
        "pid": None,
        "started_at": None,
        "query": None,
        "log_file": log_file if os.path.exists(log_file) else None,
        "results_file": results_file if os.path.exists(results_file) else None
    }
    
    # Read job info file
    if os.path.exists(job_info_file):
        with open(job_info_file, "r") as f:
            for line in f:
                if line.startswith("Job ID:"):
                    job_info["pid"] = int(line.split(":", 1)[1].strip())
                elif line.startswith("Started at:"):
                    job_info["started_at"] = line.split(":", 1)[1].strip()
                elif line.startswith("Query:"):
                    job_info["query"] = line.split(":", 1)[1].strip()
    else:
        job_info["status"] = "Unknown (no job info file)"
        return job_info
    
    # Check if process is running
    if job_info["pid"]:
        try:
            # Cross-platform way to check if process is running
            if os.name == "nt":  # Windows
                running = subprocess.call(["tasklist", "/FI", f"PID eq {job_info['pid']}"], 
                                         stdout=subprocess.DEVNULL) == 0
            else:  # Unix
                running = os.path.exists(f"/proc/{job_info['pid']}")
                
            if running:
                job_info["status"] = "Running"
                
                # Calculate runtime
                if job_info["started_at"]:
                    try:
                        started_time = datetime.strptime(job_info["started_at"], 
                                                       "%Y-%m-%d %H:%M:%S.%f")
                        runtime = datetime.now() - started_time
                        job_info["runtime"] = str(runtime).split(".")[0]  # Remove microseconds
                    except Exception:
                        pass
            else:
                if os.path.exists(results_file):
                    job_info["status"] = "Completed"
                    
                    # Get file size
                    job_info["results_size"] = f"{os.path.getsize(results_file) / 1024:.1f} KB"
                    
                    # Get completion time
                    job_info["completed_at"] = time.ctime(os.path.getmtime(results_file))
                else:
                    job_info["status"] = "Failed or stopped"
        except Exception as e:
            job_info["status"] = f"Error checking status: {e}"
    
    # Get last few log lines
    if os.path.exists(log_file):
        try:
            with open(log_file, "r") as f:
                lines = f.readlines()
                job_info["last_log_lines"] = lines[-10:] if lines else []
        except Exception as e:
            job_info["log_error"] = str(e)
    
    return job_info

def display_job_status(job_info):
    """Display job status in a readable format."""
    status_color = "\033[0m"  # Default
    
    if job_info["status"] == "Running":
        status_color = "\033[32m"  # Green
    elif job_info["status"] == "Completed":
        status_color = "\033[34m"  # Blue
    elif job_info["status"] == "Failed or stopped":
        status_color = "\033[31m"  # Red
    
    print(f"\nJob: {job_info['directory']}")
    print(f"Status: {status_color}{job_info['status']}\033[0m")
    
    if "query" in job_info and job_info["query"]:
        print(f"Query: {job_info['query']}")
    
    if "runtime" in job_info:
        print(f"Runtime: {job_info['runtime']}")
    
    if "completed_at" in job_info:
        print(f"Completed at: {job_info['completed_at']}")
        
    if "results_size" in job_info:
        print(f"Results size: {job_info['results_size']}")
    
    if "last_log_lines" in job_info and job_info["last_log_lines"]:
        print("\nRecent log entries:")
        for line in job_info["last_log_lines"]:
            # Colorize INFO, WARNING, ERROR, etc.
            if "INFO:" in line:
                line = line.replace("INFO:", "\033[32mINFO:\033[0m")
            elif "WARNING:" in line:
                line = line.replace("WARNING:", "\033[33mWARNING:\033[0m")
            elif "ERROR:" in line:
                line = line.replace("ERROR:", "\033[31mERROR:\033[0m")
            
            print(f"  {line.strip()}")
    
    print("\n" + "-" * 50)

def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Determine job directories to monitor
    job_dirs = []
    
    if args.job_dir:
        if os.path.exists(args.job_dir):
            job_dirs = [args.job_dir]
        else:
            print(f"Error: Job directory '{args.job_dir}' not found.")
            sys.exit(1)
    else:
        job_dirs = find_recent_jobs()
        
        if not job_dirs:
            print("No recent jobs found.")
            sys.exit(0)
    
    # Monitor jobs
    try:
        while True:
            # Clear screen
            os.system('cls' if os.name == 'nt' else 'clear')
            
            print(f"F.A.D.E Job Monitor - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"Monitoring {len(job_dirs)} job(s), refreshing every {args.refresh} seconds")
            
            # Get status for each job
            for job_dir in job_dirs:
                job_info = get_job_info(job_dir)
                display_job_status(job_info)
            
            print("\nPress Ctrl+C to exit")
            
            # Wait for refresh interval
            time.sleep(args.refresh)
    except KeyboardInterrupt:
        print("\nMonitoring stopped.")
        sys.exit(0)

if __name__ == "__main__":
    main()
