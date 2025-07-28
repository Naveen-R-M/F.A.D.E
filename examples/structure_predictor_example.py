#!/usr/bin/env python3
"""
Example usage of the Structure Predictor Agent
"""

import os
import json
import logging
from dotenv import load_dotenv

from agents.target_selector import TargetSelector
from agents.structure_predictor import StructurePredictor
from utils.logging import setup_logging, get_logger

# Configure logging
setup_logging("INFO")
logger = get_logger("fade.example")

# Load environment variables from .env file
load_dotenv()

def main():
    """Run the example."""
    logger.info("Starting F.A.D.E Structure Predictor Example")
    
    # Create output directories
    output_dir = "examples/output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Example query
    query = """
    I'm looking for potential drug candidates that could target the KRAS G12D mutant protein,
    which is implicated in pancreatic cancer. The molecules should be able to cross the
    blood-brain barrier, have low toxicity, and follow Lipinski's rule of five for drug-likeness.
    """
    
    logger.info("Query: %s", query)
    
    # Step 1: Run Target Selector
    logger.info("Running Target Selector...")
    target_selector = TargetSelector()
    target_selector_results = target_selector.process(query)
    
    # Save Target Selector results
    with open(os.path.join(output_dir, "target_selector_results.json"), "w") as f:
        json.dump(target_selector_results, f, indent=2)
    
    logger.info("Target Selector completed")
    
    # Step 2: Run Structure Predictor
    logger.info("Running Structure Predictor...")
    structure_predictor = StructurePredictor()
    
    structure_predictor_input = {
        "sequences": target_selector_results.get("sequences", {}),
        "job_configs": target_selector_results.get("config_files", {})
    }
    
    structure_predictor_results = structure_predictor.process(structure_predictor_input)
    
    # Save Structure Predictor results
    with open(os.path.join(output_dir, "structure_predictor_results.json"), "w") as f:
        json.dump(structure_predictor_results, f, indent=2)
    
    logger.info("Structure Predictor completed")
    
    # Print summary
    logger.info("Example completed successfully")
    logger.info("Results saved to: %s", output_dir)
    
    # Print structure information
    for target_name, structure_info in structure_predictor_results.get("structures", {}).items():
        logger.info("Structure for %s:", target_name)
        logger.info("  PDB file: %s", structure_info.get("pdb_file", "N/A"))
        logger.info("  Residue count: %d", structure_info.get("residue_count", 0))
        logger.info("  Atom count: %d", structure_info.get("atom_count", 0))
        logger.info("  Overall confidence: %.2f", structure_info.get("confidence_scores", {}).get("overall", 0.0))
    
    # Print binding site information
    for target_name, binding_sites in structure_predictor_results.get("binding_sites", {}).items():
        logger.info("Binding sites for %s:", target_name)
        for i, site in enumerate(binding_sites):
            logger.info("  Site %d: %s", i+1, site.get("name", "Unknown"))
            logger.info("    Type: %s", site.get("type", "Unknown"))
            logger.info("    Score: %.2f", site.get("score", 0.0))
            logger.info("    Residue count: %d", len(site.get("residue_ids", [])))


if __name__ == "__main__":
    main()
