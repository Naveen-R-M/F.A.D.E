#!/usr/bin/env python3
"""
Script to regenerate AlphaFold3 job scripts with the updated template
Run this after updating the job template to fix the --fasta_paths issue
"""

import os
import sys
import json
from pathlib import Path

# Add the project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from utils.config_generator import ConfigGenerator

def regenerate_alphafold_scripts():
    """Regenerate all AlphaFold3 job scripts with the updated template."""
    
    # Initialize config generator
    config_gen = ConfigGenerator()
    
    # Find existing config files
    config_dir = Path("data/inputs/configs")
    job_script_dir = Path("data/inputs/job_scripts")
    
    if not config_dir.exists():
        print(f"Config directory not found: {config_dir}")
        return
    
    # Find all AlphaFold config files
    alphafold_configs = list(config_dir.glob("*_alphafold.json"))
    
    if not alphafold_configs:
        print("No AlphaFold config files found")
        return
    
    print(f"Found {len(alphafold_configs)} AlphaFold config files to regenerate:")
    
    for config_file in alphafold_configs:
        print(f"Processing: {config_file}")
        
        # Load the config
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        protein_name = config.get('protein_name', config_file.stem.replace('_alphafold', ''))
        
        # Generate new job script
        try:
            job_script_path = os.path.join(job_script_dir, f"{protein_name}_alphafold.sh")
            config_gen.create_alphafold_job(
                protein_name=protein_name,
                sequence_file=config['sequence_file'],
                output_dir=config['output_dir'],
                config_file=str(config_file),
                job_script_path=job_script_path
            )
            
            print(f"  ✅ Generated: {job_script_path}")
            
        except Exception as e:
            print(f"  ❌ Error generating job script for {protein_name}: {e}")
    
    print("\nRegeneration complete!")
    print("\nNext steps:")
    print("1. Cancel any currently running AlphaFold3 jobs:")
    print("   scancel <job_id>")
    print("2. Resubmit the jobs with the new scripts:")
    print("   sbatch data/inputs/job_scripts/<protein>_alphafold.sh")

if __name__ == "__main__":
    regenerate_alphafold_scripts()
