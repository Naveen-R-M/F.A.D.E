#!/usr/bin/env python3
"""
Nextflow wrapper script for Molecular Docking
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
    parser = argparse.ArgumentParser(description="Run Molecular Docking")
    parser.add_argument("--molecules", required=True, help="Path to molecules.sdf")
    parser.add_argument("--receptor", required=True, help="Path to prepared_receptor.pdb")
    parser.add_argument("--binding-sites", required=True, help="Path to binding_sites.json")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--method", default="vina", choices=["vina", "glide"], help="Docking method")
    parser.add_argument("--max-poses", type=int, default=10, help="Maximum poses per molecule")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.docking")
    
    try:
        # Load input data
        with open(args.binding_sites, 'r') as f:
            binding_sites = json.load(f)
        
        # Load molecules
        with open(args.molecules.replace('.sdf', '.json'), 'r') as f:
            molecules_data = json.load(f)
        
        logger.info(f"Running {args.method} docking for {len(molecules_data)} molecules")
        
        # For now, create mock docking results since the Docking Agent isn't fully implemented
        # In a real implementation, this would call the actual docking software
        
        docking_results = []
        top_hits = []
        
        for i, mol in enumerate(molecules_data):
            # Mock docking score (in reality, this would come from Vina/Glide)
            mock_score = -7.5 + (i * 0.1)  # Gradually worse scores
            
            result = {
                "molecule_id": mol.get("id", f"MOL_{i+1:03d}"),
                "smiles": mol.get("smiles", ""),
                "docking_score": mock_score,
                "binding_site": "site_1",
                "poses": [{
                    "pose_id": 1,
                    "score": mock_score,
                    "interactions": ["hydrogen_bond", "hydrophobic"]
                }],
                "properties": mol.get("properties", {})
            }
            
            docking_results.append(result)
            
            # Keep top hits (score < -8.0)
            if mock_score < -8.0:
                top_hits.append(result)
        
        # Sort by docking score
        docking_results.sort(key=lambda x: x["docking_score"])
        top_hits.sort(key=lambda x: x["docking_score"])
        
        # Write all docking results
        with open(os.path.join(args.output_dir, "docking_results.json"), "w") as f:
            json.dump(docking_results, f, indent=2)
        
        # Write top hits (best scoring molecules)
        with open(os.path.join(args.output_dir, "top_hits.json"), "w") as f:
            json.dump(top_hits[:10], f, indent=2)  # Top 10 hits
        
        # Write docking statistics
        stats = {
            "total_molecules": len(molecules_data),
            "successful_docking": len(docking_results),
            "top_hits_count": len(top_hits),
            "best_score": docking_results[0]["docking_score"] if docking_results else None,
            "method": args.method
        }
        
        with open(os.path.join(args.output_dir, "docking_stats.json"), "w") as f:
            json.dump(stats, f, indent=2)
            
        logger.info(f"Docking completed. {len(top_hits)} top hits identified.")
        
    except Exception as e:
        logger.error(f"Docking failed: {str(e)}")
        # Write empty outputs
        with open(os.path.join(args.output_dir, "docking_results.json"), "w") as f:
            json.dump({"error": str(e), "results": []}, f, indent=2)
        with open(os.path.join(args.output_dir, "top_hits.json"), "w") as f:
            json.dump({"error": str(e), "hits": []}, f, indent=2)
        sys.exit(1)

if __name__ == "__main__":
    main()
