#!/usr/bin/env python3
"""
RCSB Target Selector Runner for Nextflow
Runs the RCSB-based target selector agent and outputs results in Nextflow-compatible format
"""

import sys
import os
import json
import argparse
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from agents.target_selector.rcsb_target_selector import RCSBTargetSelector
from utils.logging import setup_logging, get_logger


def main():
    """Main function for running RCSB target selector."""
    parser = argparse.ArgumentParser(description="RCSB Target Selector for F.A.D.E")
    parser.add_argument("--query", required=True, help="Natural language query")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level=log_level)
    logger = get_logger("fade.rcsb_target_selector_runner")
    
    logger.info("Starting RCSB Target Selector")
    logger.info(f"Query: {args.query}")
    logger.info(f"Output directory: {args.output_dir}")
    
    try:
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Initialize agent
        agent = RCSBTargetSelector(
            name="rcsb_target_selector",
            config={"data_dir": str(project_root / "data")},
            gemini_api_key=args.api_key,
            gemini_model=args.model
        )
        
        # Process the query
        logger.info("Processing query with RCSB Target Selector")
        result = agent.process(args.query)
        
        # Extract results
        parsed_data = result.get("parsed_data", {})
        structures = result.get("structures", {})
        
        # Prepare outputs for Nextflow
        if structures:
            # Use the first (primary) target
            primary_target = next(iter(structures.keys()))
            primary_structure = structures[primary_target]
            
            # Create target_info.json
            target_info = {
                "target": primary_target,
                "pdb_id": primary_structure.get("pdb_id", "unknown"),
                "source": "rcsb_pdb",
                "resolution": primary_structure.get("resolution"),
                "method": primary_structure.get("method", ""),
                "organism": primary_structure.get("organism", ""),
                "ligands": primary_structure.get("ligands", []),
                "structure_file": primary_structure.get("pdb_file", ""),
                "confidence": 1.0,  # High confidence for experimental structures
                "mutations": primary_structure.get("mutations", [])
            }
            
            target_info_path = os.path.join(args.output_dir, "target_info.json")
            with open(target_info_path, "w") as f:
                json.dump(target_info, f, indent=2)
            
            logger.info(f"Created target_info.json: {target_info_path}")
            
            # Create or copy FASTA file
            fasta_source = primary_structure.get("fasta_file")
            fasta_target = os.path.join(args.output_dir, "protein.fasta")
            
            if fasta_source and os.path.exists(fasta_source):
                # Copy existing FASTA
                import shutil
                shutil.copy2(fasta_source, fasta_target)
                logger.info(f"Copied FASTA file: {fasta_target}")
            else:
                # Create minimal FASTA from PDB info
                with open(fasta_target, "w") as f:
                    f.write(f">{primary_target}|{primary_structure.get('pdb_id', 'UNKNOWN')}|{primary_target} from RCSB PDB\n")
                    f.write("UNKNOWN\n")  # Placeholder - would need PDB parsing for actual sequence
                logger.info(f"Created placeholder FASTA: {fasta_target}")
            
            # Create PDB ID file for Nextflow (return actual PDB ID)
            pdb_id_path = os.path.join(args.output_dir, "pdb_id.txt")
            with open(pdb_id_path, "w") as f:
                f.write(primary_structure.get("pdb_id", "unknown"))
            
            # Copy PDB file to standard location
            pdb_source = primary_structure.get("pdb_file")
            pdb_target = os.path.join(args.output_dir, "structure.pdb")
            
            if pdb_source and os.path.exists(pdb_source):
                import shutil
                shutil.copy2(pdb_source, pdb_target)
                logger.info(f"Copied PDB structure: {pdb_target}")
            
        else:
            # No structures found - create error outputs
            logger.error("No structures found in RCSB - pipeline cannot continue")
            
            error_target_info = {
                "target": "unknown",
                "pdb_id": "none", 
                "source": "rcsb_failed",
                "error": "No suitable structures found in RCSB PDB"
            }
            
            target_info_path = os.path.join(args.output_dir, "target_info.json")
            with open(target_info_path, "w") as f:
                json.dump(error_target_info, f, indent=2)
            
            # Create error outputs
            with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
                f.write(">unknown|unknown|Unknown protein - no RCSB structure found\n")
                f.write("UNKNOWN\n")
            
            with open(os.path.join(args.output_dir, "requirements.json"), "w") as f:
                json.dump({"error": "No structures found"}, f, indent=2)
            
            # Exit with error - pipeline should stop here
            logger.error("RCSB target selection failed - no experimental structures available")
            sys.exit(1)
        
        # Create requirements.json from parsed data
        requirements = parsed_data.get("drug_requirements", {})
        if not requirements:
            requirements = {"binding_affinity": "< -8 kcal/mol"}
            
        requirements_path = os.path.join(args.output_dir, "requirements.json")
        with open(requirements_path, "w") as f:
            json.dump(requirements, f, indent=2)
        
        logger.info(f"Created requirements.json: {requirements_path}")
        
        # Save complete results for debugging
        complete_results_path = os.path.join(args.output_dir, "complete_results.json")
        with open(complete_results_path, "w") as f:
            json.dump(result, f, indent=2)
        
        logger.info("RCSB Target Selector completed successfully")
        
    except Exception as e:
        logger.error(f"Error in RCSB target selector: {e}")
        
        # Create error outputs for Nextflow compatibility
        error_target_info = {
            "target": "error",
            "pdb_id": "none",
            "source": "rcsb_error",
            "error": str(e)
        }
        
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump(error_target_info, f, indent=2)
            
        with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
            f.write(">error|error|Error in target selection\n")
            f.write("ERROR\n")
            
        with open(os.path.join(args.output_dir, "requirements.json"), "w") as f:
            json.dump({"error": "Target selection failed"}, f, indent=2)
            
        with open(os.path.join(args.output_dir, "pdb_id.txt"), "w") as f:
            f.write("none")
        
        sys.exit(1)


if __name__ == "__main__":
    main()
