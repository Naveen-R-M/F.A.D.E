#!/usr/bin/env python3
"""
Cleanup script for Nextflow outputs, logs, and work directories.
"""

import os
import shutil
import argparse
import time
import glob
import subprocess
from datetime import datetime, timedelta

def get_directory_age_days(directory):
    """Calculates the age of a directory in days based on its last modification time."""
    try:
        mtime = os.path.getmtime(directory)
        return (time.time() - mtime) / (24 * 3600)
    except OSError:
        return float('inf')

def get_file_age_days(file_path):
    """Calculates the age of a file in days based on its last modification time."""
    try:
        mtime = os.path.getmtime(file_path)
        return (time.time() - mtime) / (24 * 3600)
    except OSError:
        return float('inf')

def find_items(base_path, pattern, days=0, item_type='dir'):
    """Finds files or directories based on age."""
    items = []
    full_pattern = os.path.join(base_path, pattern)
    for item_path in glob.glob(full_pattern):
        is_dir = os.path.isdir(item_path)
        is_file = os.path.isfile(item_path)
        
        if (item_type == 'dir' and not is_dir) or \
           (item_type == 'file' and not is_file):
            continue

        age = get_directory_age_days(item_path) if is_dir else get_file_age_days(item_path)
        
        if age > days:
            items.append(item_path)
    return items

def run_cleanup(items_to_delete, dry_run=False):
    """Performs the actual deletion of files and directories."""
    if not items_to_delete:
        return 0
        
    deleted_count = 0
    if dry_run:
        print("\nDRY RUN active. No files or directories will be deleted.")
        return 0
    
    print("\nProceeding with deletion...")
    for item in items_to_delete:
        try:
            if os.path.isdir(item):
                shutil.rmtree(item)
                print(f"  Deleted directory: {item}")
            elif os.path.isfile(item):
                os.remove(item)
                print(f"  Deleted file: {item}")
            deleted_count += 1
        except OSError as e:
            print(f"  Error deleting {item}: {e}")
    return deleted_count

def run_nextflow_clean(nextflow_dir, days, dry_run=False):
    """Executes 'nextflow clean' to safely remove intermediate files."""
    print("\n--- Cleaning Nextflow Work Directory (using 'nextflow clean') ---")
    before_date = (datetime.now() - timedelta(days=days)).strftime('%Y-%m-%d')
    command = ["nextflow", "clean", "-f", "-before", before_date]
    
    print(f"This will execute the following command in '{nextflow_dir}':")
    print(f"  $ {' '.join(command)}")
    
    if dry_run:
        print("\nDRY RUN active. The 'nextflow clean' command will not be executed.")
        return

    try:
        result = subprocess.run(command, cwd=nextflow_dir, capture_output=True, text=True, check=True)
        print("\n'nextflow clean' executed successfully.")
        if result.stdout:
            print("Output:", result.stdout, sep='\n')
    except FileNotFoundError:
        print("\nError: 'nextflow' command not found. Is Nextflow installed and in your PATH?")
    except subprocess.CalledProcessError as e:
        print(f"\nError executing 'nextflow clean':\n{e.stderr}")

def main():
    """Main function to parse arguments and perform cleanup."""
    parser = argparse.ArgumentParser(
        description="Clean up Nextflow outputs. No args deletes all; a number deletes files older than that many days.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'days',
        type=int,
        nargs='?',
        default=0,
        help="Optional: Delete items older than this many days. If omitted, deletes ALL items."
    )
    parser.add_argument("--dry-run", action="store_true", help="Show what would be deleted without actually deleting.")
    parser.add_argument("--confirm", action="store_true", help="Skip confirmation prompt (use with caution!).")
    args = parser.parse_args()

    nextflow_dir = "nextflow"
    if not os.path.isdir(nextflow_dir):
        print(f"Error: Nextflow directory '{nextflow_dir}' not found.")
        return

    # --- Total Cleanup (days = 0) ---
    if args.days == 0:
        print("Cleanup mode: DELETE ALL Nextflow outputs.")
        results_to_delete = find_items(nextflow_dir, "results_*", -1, 'dir')
        logs_to_delete = find_items(nextflow_dir, ".nextflow.log.*", -1, 'file')
        work_dir_path = os.path.join(nextflow_dir, "work")
        work_to_delete = [work_dir_path] if os.path.exists(work_dir_path) else []
        
        items_to_delete = results_to_delete + logs_to_delete + work_to_delete
        
        if not items_to_delete:
            print("No Nextflow outputs found to delete.")
            return
            
        print("The following items will be PERMANENTLY DELETED:")
        for item in items_to_delete:
            print(f"  - {item}")

        if not args.dry_run and not args.confirm:
            response = input("\nARE YOU SURE you want to delete all listed items? This is irreversible. (yes/no): ")
            if response.lower() != 'yes':
                print("Cleanup cancelled by user.")
                return
        
        run_cleanup(items_to_delete, args.dry_run)

    # --- Partial Cleanup (days > 0) ---
    else:
        print(f"Cleanup mode: Deleting items older than {args.days} days.")
        results_to_delete = find_items(nextflow_dir, "results_*", args.days, 'dir')
        logs_to_delete = find_items(nextflow_dir, ".nextflow.log.*", args.days, 'file')
        items_to_delete = results_to_delete + logs_to_delete

        if items_to_delete:
            print("\n--- Cleaning Logs and Results ---")
            print("The following log and result items will be DELETED:")
            for item in items_to_delete:
                print(f"  - {item}")
            
            if not args.dry_run and not args.confirm:
                response = input("\nAre you sure you want to delete the logs and results listed above? (yes/no): ")
                if response.lower() != 'yes':
                    print("Cleanup of logs and results cancelled by user.")
                    items_to_delete = [] # Empty the list so cleanup isn't run
            
            if items_to_delete:
                run_cleanup(items_to_delete, args.dry_run)
        else:
            print("\nNo old logs or results found to delete.")
        
        run_nextflow_clean(nextflow_dir, args.days, args.dry_run)

    print("\nCleanup script finished.")

if __name__ == "__main__":
    main()
