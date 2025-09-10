#!/usr/bin/env python3
"""
A script to check disk usage of output and log directories in the F.A.D.E project.
"""

import os
import argparse
from typing import Dict, List, Tuple

def get_directory_size(path: str) -> int:
    """
    Calculate the total size of a directory in bytes.
    
    Args:
        path: Directory path
        
    Returns:
        Total size in bytes
    """
    total_size = 0
    if not os.path.exists(path):
        return 0
        
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            if os.path.exists(filepath):
                try:
                    total_size += os.path.getsize(filepath)
                except OSError:
                    pass
    return total_size

def format_size(size_bytes: int) -> str:
    """
    Format size in bytes to human-readable string.
    
    Args:
        size_bytes: Size in bytes
        
    Returns:
        Formatted string (e.g., "1.5 GB")
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"

def count_files(path: str) -> int:
    """
    Count the number of files in a directory (recursively).
    
    Args:
        path: Directory path
        
    Returns:
        Number of files
    """
    if not os.path.exists(path):
        return 0
        
    count = 0
    for dirpath, dirnames, filenames in os.walk(path):
        count += len(filenames)
    return count

def check_disk_usage(detailed: bool = False) -> Dict[str, Tuple[int, int]]:
    """
    Check disk usage of various directories in the project.
    
    Args:
        detailed: If True, show detailed breakdown
        
    Returns:
        Dictionary mapping directory paths to (size, file_count) tuples
    """
    directories = {
        "Logs": [
            "logs",
            "results/logs",
            "current_test/logs",
            "summary_test/logs",
        ],
        "Outputs": [
            "agents/data/outputs",
            "data/outputs",
        ],
        "Results": [
            "results",
            "test_results",
            "current_test",
            "summary_test",
        ],
        "Data": [
            "data/inputs",
            "data/memory",
            "data/references",
            "agents/data/inputs",
        ],
    }
    
    usage_info = {}
    category_totals = {}
    
    print("=" * 70)
    print("F.A.D.E Project Disk Usage Report")
    print("=" * 70)
    print()
    
    grand_total_size = 0
    grand_total_files = 0
    
    for category, dir_list in directories.items():
        print(f"\n{category}:")
        print("-" * 40)
        
        category_size = 0
        category_files = 0
        
        for directory in dir_list:
            if os.path.exists(directory):
                size = get_directory_size(directory)
                files = count_files(directory)
                usage_info[directory] = (size, files)
                category_size += size
                category_files += files
                
                if detailed or size > 0:
                    print(f"  {directory:<35} {format_size(size):>12} ({files:,} files)")
                    
                    # Show subdirectories for large directories
                    if detailed and size > 10 * 1024 * 1024:  # > 10 MB
                        subdirs = []
                        try:
                            for item in os.listdir(directory):
                                item_path = os.path.join(directory, item)
                                if os.path.isdir(item_path):
                                    subdir_size = get_directory_size(item_path)
                                    if subdir_size > 1024 * 1024:  # > 1 MB
                                        subdirs.append((item, subdir_size))
                            
                            if subdirs:
                                subdirs.sort(key=lambda x: x[1], reverse=True)
                                for subdir_name, subdir_size in subdirs[:5]:
                                    print(f"    └─ {subdir_name:<31} {format_size(subdir_size):>12}")
                        except OSError:
                            pass
            else:
                if detailed:
                    print(f"  {directory:<35} {'(not found)':>12}")
        
        category_totals[category] = (category_size, category_files)
        grand_total_size += category_size
        grand_total_files += category_files
        
        if category_size > 0:
            print(f"  {'TOTAL':.<35} {format_size(category_size):.>12} ({category_files:,} files)")
    
    # Show summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    for category, (size, files) in category_totals.items():
        if size > 0:
            percentage = (size / grand_total_size * 100) if grand_total_size > 0 else 0
            print(f"{category:<20} {format_size(size):>15} ({percentage:5.1f}%) - {files:,} files")
    
    print("-" * 70)
    print(f"{'GRAND TOTAL':<20} {format_size(grand_total_size):>15} - {grand_total_files:,} files")
    
    # Find largest files
    print("\n" + "=" * 70)
    print("TOP 5 LARGEST FILES")
    print("=" * 70)
    
    largest_files = []
    for directory in sum(directories.values(), []):
        if os.path.exists(directory):
            for dirpath, dirnames, filenames in os.walk(directory):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)
                    try:
                        size = os.path.getsize(filepath)
                        if size > 1024 * 1024:  # Only consider files > 1 MB
                            largest_files.append((filepath, size))
                    except OSError:
                        pass
    
    largest_files.sort(key=lambda x: x[1], reverse=True)
    for filepath, size in largest_files[:5]:
        # Shorten the path for display
        display_path = filepath
        if len(display_path) > 50:
            display_path = "..." + filepath[-47:]
        print(f"{display_path:<55} {format_size(size):>12}")
    
    return usage_info

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Check disk usage of F.A.D.E project directories.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--detailed", "-d",
        action="store_true",
        help="Show detailed breakdown including subdirectories."
    )
    
    args = parser.parse_args()
    
    check_disk_usage(args.detailed)
    
    print("\n" + "=" * 70)
    print("To clean up disk space, you can use:")
    print("  python clear_logs.py      # Remove log files")
    print("  python clear_outputs.py   # Remove output files")
    print("  python clear_outputs.py --keep-best  # Keep best models")
    print("=" * 70)

if __name__ == "__main__":
    main()
