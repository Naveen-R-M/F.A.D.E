#!/usr/bin/env python3
"""
Simple RCSB Target Selector Runner - No LLM Dependencies
Quick test runner for RCSB functionality
"""

import sys
import os
import json
import argparse
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from utils.simple_rcsb_selector import SimpleRCSBTargetSelector


def main():
    """Main function for simple RCSB target selector."""
    parser = argparse.ArgumentParser(description="Simple RCSB Target Selector")
    parser.add_argument("--query", required=True, help="Natural language query")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key (not used)")
    parser.add_argument("--model", help="Gemini model (not used)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    print("=== SIMPLE RCSB TARGET SELECTOR ===")
    print(f"Query: {args.query}")
    print(f"Output directory: {args.output_dir}")
    
    try:
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Initialize simple selector
        selector = SimpleRCSBTargetSelector()
        
        # Process the query
        result = selector.process(args.query)
        
        # Extract results
        parsed_data = result.get("parsed_data", {})
        structures = result.get("structures", {})
        
        if structures:
            # Use the first target
            primary_target = next(iter(structures.keys()))
            primary_structure = structures[primary_target]
            
            # Create target_info.json
            target_info = {
                "target": primary_target,
                "pdb_id": primary_structure.get("pdb_id", "unknown"),
                "source": "simple_rcsb",
                "resolution": primary_structure.get("resolution"),
                "method": primary_structure.get("method", ""),
                "organism": primary_structure.get("organism", ""),
                "confidence": 1.0
            }
            
            target_info_path = os.path.join(args.output_dir, "target_info.json")
            with open(target_info_path, "w") as f:
                json.dump(target_info, f, indent=2)
            
            print(f"✓ Created target_info.json")
            
            # Create simple FASTA (placeholder for now)
            fasta_path = os.path.join(args.output_dir, "protein.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{primary_target}|{primary_structure.get('pdb_id', 'UNKNOWN')}|{primary_target} from RCSB PDB\\n")
                f.write("PLACEHOLDER_SEQUENCE\\n")  # Will extract from PDB later
            
            print(f"✓ Created protein.fasta")
            
            # Copy PDB structure
            pdb_source = primary_structure.get("pdb_file")
            if pdb_source and os.path.exists(pdb_source):
                import shutil
                pdb_target = os.path.join(args.output_dir, "structure.pdb")
                shutil.copy2(pdb_source, pdb_target)
                print(f"✓ Copied PDB structure: {pdb_target}")
                
        else:
            raise Exception("No structures found")
        
        # Create requirements.json
        requirements = parsed_data.get("drug_requirements", {})
        requirements_path = os.path.join(args.output_dir, "requirements.json")
        with open(requirements_path, "w") as f:
            json.dump(requirements, f, indent=2)
        
        print(f"✓ Created requirements.json")
        print("✓ Simple RCSB Target Selector completed successfully")
        
    except Exception as e:
        print(f"✗ Error: {e}")
        
        # Create error outputs
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump({"error": str(e), "target": "unknown"}, f)
            
        with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
            f.write(">error|error|Error\\n")
            f.write("ERROR\\n")
            
        with open(os.path.join(args.output_dir, "requirements.json"), "w") as f:
            json.dump({"error": str(e)}, f)
        
        sys.exit(1)


if __name__ == "__main__":
    main()
