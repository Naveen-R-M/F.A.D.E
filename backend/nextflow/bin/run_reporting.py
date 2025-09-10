#!/usr/bin/env python3
"""
Nextflow wrapper script for Final Reporting
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime

# Add the parent directory to Python path for imports
fade_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(fade_root))

from utils.logging import setup_logging, get_logger

def main():
    parser = argparse.ArgumentParser(description="Generate Final Report")
    parser.add_argument("--target-info", required=True, help="Path to target_info.json")
    parser.add_argument("--structure", required=True, help="Path to structure.pdb")
    parser.add_argument("--docking-results", required=True, help="Path to docking_results.json")
    parser.add_argument("--optimized-leads", required=True, help="Path to optimized_leads.json")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.reporting")
    
    try:
        # Load all input data
        with open(args.target_info, 'r') as f:
            target_info = json.load(f)
        
        with open(args.docking_results, 'r') as f:
            docking_results = json.load(f)
        
        with open(args.optimized_leads, 'r') as f:
            optimized_leads = json.load(f)
        
        logger.info("Generating final report")
        
        # Create comprehensive report
        report = {
            "timestamp": datetime.now().isoformat(),
            "target": target_info,
            "summary": {
                "total_molecules_evaluated": len(docking_results),
                "top_candidates": len(optimized_leads),
                "best_score": min([r.get("docking_score", 0) for r in docking_results]) if docking_results else None
            },
            "top_candidates": optimized_leads[:5],  # Top 5 candidates
            "methodology": {
                "structure_prediction": "AlphaFold3",
                "molecule_generation": "LLM-guided + RDKit",
                "docking": "AutoDock Vina",
                "optimization": "Iterative refinement"
            }
        }
        
        # Write JSON report
        with open(os.path.join(args.output_dir, "final_report.json"), "w") as f:
            json.dump(report, f, indent=2)
        
        # Generate human-readable report
        target_name = target_info.get("target", "Unknown")
        best_candidates = optimized_leads[:3] if optimized_leads else []
        
        markdown_report = f"""# F.A.D.E Drug Discovery Report
        
## Target: {target_name}

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

### Summary
- **Total molecules evaluated:** {len(docking_results)}
- **Top candidates identified:** {len(optimized_leads)}
- **Best docking score:** {report['summary']['best_score']:.2f} kcal/mol

### Top 3 Candidates

"""
        
        for i, candidate in enumerate(best_candidates, 1):
            score = candidate.get("optimized_score", candidate.get("docking_score", "N/A"))
            smiles = candidate.get("smiles", "N/A")
            mol_id = candidate.get("id", f"Candidate_{i}")
            
            markdown_report += f"""#### {mol_id}
- **SMILES:** `{smiles}`
- **Docking Score:** {score} kcal/mol
- **Properties:** {json.dumps(candidate.get('properties', {}), indent=2)}

"""
        
        markdown_report += f"""
### Methodology
- **Structure Prediction:** {report['methodology']['structure_prediction']}
- **Molecule Generation:** {report['methodology']['molecule_generation']}
- **Docking:** {report['methodology']['docking']}
- **Optimization:** {report['methodology']['optimization']}

### Files Generated
- `final_report.json` - Detailed results in JSON format
- `top_candidates.sdf` - 3D structures of top candidates
- `analysis_summary.txt` - This human-readable summary
"""
        
        with open(os.path.join(args.output_dir, "analysis_summary.md"), "w") as f:
            f.write(markdown_report)
        
        # Create SDF file with top candidates
        sdf_content = ""
        for candidate in best_candidates:
            mol_id = candidate.get("id", "MOL_001")
            smiles = candidate.get("smiles", "")
            score = candidate.get("optimized_score", candidate.get("docking_score", 0))
            
            sdf_content += f"{mol_id}\n"
            sdf_content += f"  F.A.D.E Generated - Score: {score}\n"
            sdf_content += f"  SMILES: {smiles}\n"
            sdf_content += f"$$$$\n"
        
        with open(os.path.join(args.output_dir, "top_candidates.sdf"), "w") as f:
            f.write(sdf_content)
            
        logger.info("Final report generated successfully")
        
    except Exception as e:
        logger.error(f"Report generation failed: {str(e)}")
        # Write error report
        error_report = {
            "error": str(e),
            "timestamp": datetime.now().isoformat(),
            "status": "failed"
        }
        with open(os.path.join(args.output_dir, "final_report.json"), "w") as f:
            json.dump(error_report, f, indent=2)
        sys.exit(1)

if __name__ == "__main__":
    main()
