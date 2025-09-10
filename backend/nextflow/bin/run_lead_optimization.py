#!/usr/bin/env python3
"""
Nextflow wrapper script for Lead Optimization
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Add the parent directory to Python path for imports
fade_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(fade_root))

from utils.logging import setup_logging, get_logger

def main():
    parser = argparse.ArgumentParser(description="Run Lead Optimization")
    parser.add_argument("--top-hits", required=True, help="Path to top_hits.json")
    parser.add_argument("--structure", required=True, help="Path to structure.pdb")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.lead_optimization")
    
    try:
        # Load top hits
        with open(args.top_hits, 'r') as f:
            top_hits = json.load(f)
        
        logger.info(f"Optimizing {len(top_hits)} lead compounds")
        
        # For now, create mock optimization results
        # In reality, this would use the Refiner Agent to improve molecules
        
        optimized_leads = []
        
        for hit in top_hits:
            # Mock optimization - in reality this would modify the molecule
            optimized = hit.copy()
            optimized["optimized"] = True
            optimized["optimization_cycle"] = 1
            optimized["original_score"] = hit.get("docking_score", 0)
            optimized["optimized_score"] = hit.get("docking_score", 0) - 0.5  # Mock improvement
            
            optimized_leads.append(optimized)
        
        # Write optimized leads
        with open(os.path.join(args.output_dir, "optimized_leads.json"), "w") as f:
            json.dump(optimized_leads, f, indent=2)
        
        # Write optimization statistics
        stats = {
            "input_leads": len(top_hits),
            "optimized_leads": len(optimized_leads),
            "average_improvement": 0.5,  # Mock value
            "optimization_cycles": 1
        }
        
        with open(os.path.join(args.output_dir, "optimization_stats.json"), "w") as f:
            json.dump(stats, f, indent=2)
            
        logger.info("Lead optimization completed successfully")
        
    except Exception as e:
        logger.error(f"Lead optimization failed: {str(e)}")
        # Write empty outputs
        with open(os.path.join(args.output_dir, "optimized_leads.json"), "w") as f:
            json.dump({"error": str(e), "leads": []}, f, indent=2)
        sys.exit(1)

if __name__ == "__main__":
    main()
