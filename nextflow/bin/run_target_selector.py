#!/usr/bin/env python3
"""
Nextflow wrapper script for Target Selector Agent
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Add the parent directory to Python path for imports
fade_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(fade_root))

from agents.target_selector.target_selector import TargetSelector
from utils.logging import setup_logging, get_logger

def main():
    parser = argparse.ArgumentParser(description="Run Target Selector Agent")
    parser.add_argument("--query", required=True, help="Natural language query")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.target_selector")
    
    try:
        # Initialize target selector
        target_selector = TargetSelector(
            gemini_api_key=args.api_key,
            gemini_model=args.model
        )
        
        # Process the query
        logger.info(f"Processing query: {args.query}")
        result = target_selector.process(args.query)
        
        # Extract outputs
        protein_targets = result.get("protein_targets", [])
        sequences = result.get("sequences", {})
        configs = result.get("config_files", {})
        
        # Write target_info.json
        if protein_targets:
            target_info = {
                "target": protein_targets[0].get("name", "unknown"),
                "uniprot_id": protein_targets[0].get("uniprot_id", "unknown"),
                "pdb_id": protein_targets[0].get("pdb_id", "none"),
                "description": protein_targets[0].get("description", ""),
                "mutation": protein_targets[0].get("mutation", ""),
                "all_targets": protein_targets
            }
        else:
            target_info = {"error": "No targets identified"}
            
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump(target_info, f, indent=2)
        
        # Write protein.fasta (use first target)
        if sequences and protein_targets:
            first_target = protein_targets[0].get("name")
            if first_target and first_target in sequences:
                seq_info = sequences[first_target]
                with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
                    f.write(f">{seq_info.get('uniprot_id', 'unknown')}|{first_target}\n")
                    f.write(f"{seq_info.get('sequence', '')}\n")
        
        # Write requirements.json
        requirements = result.get("drug_requirements", {})
        with open(os.path.join(args.output_dir, "requirements.json"), "w") as f:
            json.dump(requirements, f, indent=2)
        
        # Write full result for debugging
        with open(os.path.join(args.output_dir, "target_selector_full_result.json"), "w") as f:
            json.dump(result, f, indent=2)
            
        logger.info("Target selection completed successfully")
        
    except Exception as e:
        logger.error(f"Target selection failed: {str(e)}")
        # Write error information
        error_info = {"error": str(e), "target": "unknown"}
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump(error_info, f, indent=2)
        sys.exit(1)

if __name__ == "__main__":
    main()
