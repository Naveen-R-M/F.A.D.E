#!/usr/bin/env python3
"""
Regenerate AF3 job scripts with updated specifications
- Use GPU: #SBATCH --gres=gpu:h200:1
- Time: #SBATCH --time=04:00:00
- Memory: 32G (increased from 16G)
- Logs: Output/error to F.A.D.E/logs/ directory
- DATA_PATH and MODEL_PATH environment variables for container
"""

import os
import sys
import json
from pathlib import Path

# Add the F.A.D.E root directory to Python path
fade_root = Path(__file__).parent
sys.path.insert(0, str(fade_root))

from utils.config_generator import ConfigGenerator


def regenerate_job_scripts():
    """Regenerate all existing job scripts with updated specifications."""
    
    # Initialize config generator
    config_gen = ConfigGenerator()
    
    # Paths
    fade_root = Path(__file__).parent
    configs_dir = fade_root / "data" / "inputs" / "configs"
    job_scripts_dir = fade_root / "data" / "inputs" / "job_scripts"
    logs_dir = fade_root / "logs"
    
    # Ensure logs directory exists
    logs_dir.mkdir(exist_ok=True)
    
    print(f"Regenerating job scripts with updated specifications...")
    print(f"Logs directory: {logs_dir}")
    print(f"Config directory: {configs_dir}")
    print(f"Job scripts directory: {job_scripts_dir}")
    
    # Find all existing config files
    config_files = list(configs_dir.glob("*_alphafold.json"))
    
    for config_file in config_files:
        try:
            # Load config
            with open(config_file, 'r') as f:
                config = json.load(f)
            
            protein_name = config.get("protein_name")
            if not protein_name:
                print(f"Skipping {config_file} - no protein name found")
                continue
            
            # Define paths
            sequence_file = config.get("sequence_file")
            output_dir = config.get("output_dir")
            job_script_path = job_scripts_dir / f"{protein_name}_alphafold.sh"
            
            print(f"Regenerating job script for {protein_name}...")
            
            # Regenerate job script with new specifications
            config_gen.create_alphafold_job(
                protein_name=protein_name,
                sequence_file=sequence_file,
                output_dir=output_dir,
                config_file=str(config_file),
                job_script_path=str(job_script_path),
                partition="gpu",
                memory="32G",  # Increased memory
                cpus=8,
                gpus=1,
                time="04:00:00",  # Updated time
                logs_dir=str(logs_dir),  # Use F.A.D.E/logs directory
                email=None
            )
            
            print(f"✓ Generated: {job_script_path}")
            
        except Exception as e:
            print(f"✗ Error processing {config_file}: {e}")
    
    print("\nJob script regeneration completed!")
    print(f"All logs will now be written to: {logs_dir}")
    print("Updated specifications:")
    print("  - GPU: H200 (--gres=gpu:h200:1)")
    print("  - Time: 4 hours (--time=04:00:00)")
    print("  - Memory: 32G (increased from 16G)")
    print("  - DATA_PATH and MODEL_PATH environment variables set")
    print("  - Container execution with proper bind mounts")


if __name__ == "__main__":
    regenerate_job_scripts()
