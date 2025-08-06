#!/usr/bin/env python3
"""
Monitor all F.A.D.E jobs across different result directories.
This is a global monitoring script, different from the job-specific check_status.py files.
"""

import os
import sys
import glob
import json
import time
from datetime import datetime
from typing import List, Dict, Tuple

def find_fade_jobs(base_dir: str = ".") -> List[str]:
    """
    Find all F.A.D.E job directories containing job_info.txt files.
    
    Args:
        base_dir: Base directory to search from
        
    Returns:
        List of directories containing F.A.D.E jobs
    """
    job_dirs = []
    
    # Search patterns for job directories
    search_patterns = [
        "results*/*/job_info.txt",
        "results*/job_info.txt",
        "test*/*/job_info.txt",
        "test*/job_info.txt",
        "*/results*/job_info.txt"
    ]
    
    for pattern in search_patterns:
        for job_info_path in glob.glob(os.path.join(base_dir, pattern)):
            job_dir = os.path.dirname(job_info_path)
            if job_dir not in job_dirs:
                job_dirs.append(job_dir)
    
    return sorted(job_dirs)

def read_job_info(job_dir: str) -> Dict[str, str]:
    """
    Read job information from job_info.txt file.
    
    Args:
        job_dir: Directory containing the job
        
    Returns:
        Dictionary with job information
    """
    info = {
        "directory": job_dir,
        "pid": None,
        "start_time": None,
        "query": None,
        "status": "UNKNOWN"
    }
    
    job_info_file = os.path.join(job_dir, "job_info.txt")
    
    try:
        with open(job_info_file, "r") as f:
            for line in f:
                if line.startswith("Job ID:"):
                    info["pid"] = line.split(":", 1)[1].strip()
                elif line.startswith("Start time:"):
                    info["start_time"] = line.split(":", 1)[1].strip()
                elif line.startswith("Query:"):
                    info["query"] = line.split(":", 1)[1].strip()[:50] + "..."
    except Exception as e:
        info["status"] = f"ERROR: {e}"
    
    return info

def is_process_running(pid: str) -> bool:
    """
    Check if a process is running.
    
    Args:
        pid: Process ID
        
    Returns:
        True if process is running
    """
    try:
        pid_int = int(pid)
        if os.name == "nt":  # Windows
            import subprocess
            return subprocess.call(["tasklist", "/FI", f"PID eq {pid_int}"], 
                                 stdout=subprocess.DEVNULL) == 0
        else:  # Unix
            return os.path.exists(f"/proc/{pid_int}")
    except (ValueError, TypeError):
        return False

def check_job_status(job_dir: str) -> str:
    """
    Check the status of a F.A.D.E job.
    
    Args:
        job_dir: Directory containing the job
        
    Returns:
        Status string
    """
    # Check for results file
    results_file = os.path.join(job_dir, "results.json")
    if os.path.exists(results_file):
        try:
            with open(results_file, "r") as f:
                data = json.load(f)
                if data:
                    return "COMPLETED"
        except:
            pass
    
    # Check for error indicators in log
    log_file = os.path.join(job_dir, "fade.log")
    if os.path.exists(log_file):
        try:
            with open(log_file, "r") as f:
                content = f.read()
                if "ERROR" in content or "Traceback" in content:
                    return "ERROR"
                if "Complete" in content or "Finished" in content:
                    return "COMPLETED"
        except:
            pass
    
    return "RUNNING"

def get_last_log_lines(job_dir: str, n: int = 3) -> List[str]:
    """
    Get the last n lines from the job log.
    
    Args:
        job_dir: Directory containing the job
        n: Number of lines to retrieve
        
    Returns:
        List of last log lines
    """
    log_file = os.path.join(job_dir, "fade.log")
    
    if not os.path.exists(log_file):
        return []
    
    try:
        with open(log_file, "r") as f:
            lines = f.readlines()
            return [line.strip() for line in lines[-n:] if line.strip()]
    except:
        return []

def format_time_ago(time_str: str) -> str:
    """
    Format a time string to show how long ago it was.
    
    Args:
        time_str: Time string
        
    Returns:
        Formatted string like "2h ago"
    """
    try:
        # Parse various time formats
        for fmt in ["%Y-%m-%d %H:%M:%S", "%a %b %d %H:%M:%S %Y"]:
            try:
                start_time = datetime.strptime(time_str.strip(), fmt)
                break
            except:
                continue
        else:
            return time_str
        
        delta = datetime.now() - start_time
        
        if delta.days > 0:
            return f"{delta.days}d ago"
        elif delta.seconds > 3600:
            return f"{delta.seconds // 3600}h ago"
        elif delta.seconds > 60:
            return f"{delta.seconds // 60}m ago"
        else:
            return "just now"
    except:
        return time_str

def main():
    """Main entry point for the monitoring script."""
    print("=" * 80)
    print("F.A.D.E Job Monitor")
    print("=" * 80)
    print(f"Scan time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Find all job directories
    job_dirs = find_fade_jobs()
    
    if not job_dirs:
        print("No F.A.D.E jobs found.")
        return
    
    print(f"Found {len(job_dirs)} job(s):\n")
    
    # Check each job
    for i, job_dir in enumerate(job_dirs, 1):
        info = read_job_info(job_dir)
        
        # Determine status
        if info["pid"]:
            if is_process_running(info["pid"]):
                info["status"] = check_job_status(job_dir)
            else:
                # Process not running, check if completed
                if check_job_status(job_dir) == "COMPLETED":
                    info["status"] = "COMPLETED"
                else:
                    info["status"] = "STOPPED"
        
        # Format output
        print(f"Job #{i}: {os.path.basename(job_dir)}")
        print(f"  Directory: {job_dir}")
        print(f"  PID: {info['pid'] or 'N/A'}")
        
        if info["start_time"]:
            time_ago = format_time_ago(info["start_time"])
            print(f"  Started: {time_ago}")
        
        # Status with color coding (if terminal supports it)
        status = info["status"]
        if sys.stdout.isatty():  # Check if output is to terminal
            if status == "RUNNING":
                status_str = f"\033[92m● {status}\033[0m"  # Green
            elif status == "COMPLETED":
                status_str = f"\033[94m✓ {status}\033[0m"  # Blue
            elif status in ["ERROR", "STOPPED"]:
                status_str = f"\033[91m✗ {status}\033[0m"  # Red
            else:
                status_str = f"  {status}"
        else:
            status_str = status
        
        print(f"  Status: {status_str}")
        
        if info["query"]:
            print(f"  Query: {info['query']}")
        
        # Show last log lines for running jobs
        if info["status"] == "RUNNING":
            last_lines = get_last_log_lines(job_dir)
            if last_lines:
                print("  Recent activity:")
                for line in last_lines:
                    print(f"    > {line[:70]}...")
        
        print()
    
    print("-" * 80)
    print("Commands:")
    print("  View specific job:  cd <directory> && python check_status.py")
    print("  Clear old outputs:  python clear_outputs.py --days 7")
    print("  Check disk usage:   python check_disk_usage.py")
    print("=" * 80)

if __name__ == "__main__":
    main()
