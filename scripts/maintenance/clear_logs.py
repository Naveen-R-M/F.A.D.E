#!/usr/bin/env python3
"""
A script to clear log files from specified directories within the project.
"""

import os
import argparse
import time
from datetime import datetime, timedelta

def clear_logs(days: int = 0):
    """
    Deletes log files from predefined directories.

    Args:
        days (int): If greater than 0, only deletes files older than
                    this number of days. If 0, deletes all files.
    """
    log_dirs = [
        "logs",
        "results/logs",
        "current_test/logs",
        "summary_test/logs",
    ]

    now = time.time()
    deleted_count = 0
    checked_count = 0

    if days > 0:
        print(f"Deleting log files older than {days} days...")
        time_threshold = now - timedelta(days=days).total_seconds()
    else:
        print("Deleting all log files...")
        time_threshold = now

    for log_dir in log_dirs:
        if not os.path.isdir(log_dir):
            print(f"Directory not found, skipping: {log_dir}")
            continue

        print(f"Checking directory: {log_dir}")
        for filename in os.listdir(log_dir):
            file_path = os.path.join(log_dir, filename)

            if os.path.isfile(file_path):
                checked_count += 1
                file_mod_time = os.path.getmtime(file_path)

                if days == 0 or file_mod_time < time_threshold:
                    try:
                        os.remove(file_path)
                        print(f"  Deleted: {file_path}")
                        deleted_count += 1
                    except OSError as e:
                        print(f"  Error deleting {file_path}: {e}")

    print(f"\nLog cleanup complete.")
    print(f"Checked {checked_count} files and deleted {deleted_count} files.")

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Clear log files from project directories.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--days",
        type=int,
        default=0,
        help="""Delete logs older than this many days.
If 0 or not specified, all logs will be deleted."""
    )
    args = parser.parse_args()

    clear_logs(args.days)

if __name__ == "__main__":
    main()
