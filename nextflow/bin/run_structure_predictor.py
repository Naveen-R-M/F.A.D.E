#!/usr/bin/env python3
"""
Nextflow wrapper script for Structure Predictor Agent
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Add the parent directory to Python path for imports
fade_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(fade_root))

from agents.structure_predictor.structure_predictor import StructurePredictor
from utils.logging import setup_logging, get_logger

def main():
    parser = argparse.ArgumentParser(description="Run Structure Predictor Agent")
    parser.add_argument("--target-info", required=True, help="Path to target_info.json")
    parser.add_argument("--fasta-file", required=True, help="Path to protein.fasta")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    parser.add_argument("--use-alphafold", action="store_true", default=True, help="Use AlphaFold3 for prediction")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.structure_predictor")
    
    try:
        # Load input data
        with open(args.target_info, 'r') as f:
            target_info = json.load(f)
        
        # Read FASTA file
        with open(args.fasta_file, 'r') as f:
            fasta_content = f.read()
        
        # Parse FASTA to get sequence
        lines = fasta_content.strip().split('\n')
        header = lines[0] if lines else ""
        sequence = ''.join(lines[1:]) if len(lines) > 1 else ""
        
        # Prepare input for structure predictor
        target_name = target_info.get("target", "unknown")
        input_data = {
            "protein_targets": [target_info],
            "sequences": {
                target_name: {
                    "uniprot_id": target_info.get("uniprot_id", ""),
                    "sequence": sequence,
                    "header": header
                }
            },
            "job_configs": {}  # Will be populated by the agent
        }
        
        # Initialize structure predictor
        structure_predictor = StructurePredictor(
            gemini_api_key=args.api_key,
            gemini_model=args.model
        )
        
        # Process the structure prediction
        logger.info(f"Processing structure prediction for target: {target_name}")
        result = structure_predictor.process(input_data)
        
        # Write outputs
        structures = result.get("structures", {})
        binding_sites = result.get("binding_sites", {})
        prepared_structures = result.get("prepared_structures", {})
        
        # Write structure.pdb (use first available structure)
        if structures:
            first_structure_info = list(structures.values())[0]
            structure_path = first_structure_info.get("structure_path")
            if structure_path and os.path.exists(structure_path):
                import shutil
                shutil.copy(structure_path, os.path.join(args.output_dir, "structure.pdb"))
            else:
                # Create a placeholder if structure prediction failed
                with open(os.path.join(args.output_dir, "structure.pdb"), "w") as f:
                    f.write("# Structure prediction failed or not available\n")
        
        # Write binding_sites.json
        with open(os.path.join(args.output_dir, "binding_sites.json"), "w") as f:
            json.dump(binding_sites, f, indent=2)
        
        # Write prepared_receptor.pdb (for docking)
        if prepared_structures:
            first_prepared = list(prepared_structures.values())[0]
            prepared_path = first_prepared.get("prepared_path")
            if prepared_path and os.path.exists(prepared_path):
                import shutil
                shutil.copy(prepared_path, os.path.join(args.output_dir, "prepared_receptor.pdb"))
            else:
                # Copy the original structure as fallback
                import shutil
                shutil.copy(os.path.join(args.output_dir, "structure.pdb"), 
                           os.path.join(args.output_dir, "prepared_receptor.pdb"))
        
        # Write full result for debugging
        with open(os.path.join(args.output_dir, "structure_predictor_full_result.json"), "w") as f:
            json.dump(result, f, indent=2)
            
        logger.info("Structure prediction completed successfully")
        
    except Exception as e:
        logger.error(f"Structure prediction failed: {str(e)}")
        # Write placeholder outputs
        with open(os.path.join(args.output_dir, "structure.pdb"), "w") as f:
            f.write("# Structure prediction failed\n")
        with open(os.path.join(args.output_dir, "binding_sites.json"), "w") as f:
            json.dump({"error": str(e)}, f, indent=2)
        with open(os.path.join(args.output_dir, "prepared_receptor.pdb"), "w") as f:
            f.write("# Structure preparation failed\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
