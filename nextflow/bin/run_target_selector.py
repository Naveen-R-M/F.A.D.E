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

def load_env_file(env_path=None):
    """Load environment variables from .env file"""
    if env_path is None:
        env_path = fade_root / ".env"
    
    if env_path.exists():
        with open(env_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    os.environ[key.strip()] = value.strip()

def main():
    parser = argparse.ArgumentParser(description="Run Target Selector Agent")
    parser.add_argument("--query", required=True, help="Natural language query")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    
    args = parser.parse_args()
    
    # Load .env file first
    load_env_file()
    
    # Use API key from arguments, or environment, or .env file
    api_key = args.api_key or os.getenv('GEMINI_API_KEY')
    
    try:
        from utils.logging import setup_logging, get_logger
        
        # Setup logging
        setup_logging(log_level="INFO")
        logger = get_logger("nextflow.target_selector")
        
        # Validate API key
        if not api_key:
            logger.error("Gemini API key not provided. Set GEMINI_API_KEY environment variable or provide it directly.")
            sys.exit(1)
        
        logger.info(f"API key loaded successfully (length: {len(api_key)})")
        
        # Initialize target selector
        from agents.target_selector.target_selector import TargetSelector
        
        target_selector = TargetSelector(
            gemini_api_key=api_key,
            gemini_model=args.model
        )
        
        # Process the query
        logger.info(f"Processing query: {args.query}")
        result = target_selector.process(args.query)
        
        # Extract outputs - CORRECTED DATA STRUCTURE
        parsed_data = result.get("parsed_data", {})
        protein_targets = parsed_data.get("protein_targets", [])
        sequences = result.get("sequences", {})
        configs = result.get("config_files", {})
        
        print(f"DEBUG: Found {len(protein_targets)} protein targets")
        print(f"DEBUG: Found {len(sequences)} sequences")
        
        # Write target_info.json
        if protein_targets:
            first_target = protein_targets[0]
            target_name = first_target.get("name", "unknown")
            
            # Get UniProt accession from sequences data
            uniprot_id = "unknown"
            if target_name in sequences:
                seq_data = sequences[target_name]
                uniprot_id = seq_data.get("accession", "unknown")  # Use 'accession', not 'uniprot_id'
                print(f"DEBUG: Found accession {uniprot_id} for {target_name}")
            
            target_info = {
                "target": target_name,
                "uniprot_id": uniprot_id,
                "pdb_id": first_target.get("pdb_id", "none"),
                "description": first_target.get("description", ""),
                "mutation": first_target.get("mutation", ""),
                "mutations": first_target.get("mutations", []),
                "all_targets": protein_targets,
                "success": True
            }
        else:
            target_info = {"error": "No targets identified", "target": "unknown", "success": False}
            
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump(target_info, f, indent=2)
        
        # Write protein.fasta - FIXED FASTA HEADER GENERATION
        if sequences and protein_targets:
            first_target_name = protein_targets[0].get("name")
            print(f"DEBUG: Processing FASTA for target: {first_target_name}")
            
            if first_target_name and first_target_name in sequences:
                seq_info = sequences[first_target_name]
                print(f"DEBUG: Sequence info keys: {list(seq_info.keys())}")
                
                # Extract all the information correctly
                accession = seq_info.get('accession', 'unknown')
                sequence = seq_info.get('sequence', '')
                protein_name = seq_info.get('protein_name', first_target_name)
                gene_name = seq_info.get('gene_name', first_target_name)
                organism = seq_info.get('organism', '')
                
                print(f"DEBUG: Accession: {accession}, Protein: {protein_name}, Gene: {gene_name}")
                
                with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
                    # Format: >UniProt_ID|Gene_Name|Protein_Name [Organism]
                    if organism:
                        header = f">{accession}|{gene_name}|{protein_name} [{organism}]"
                    else:
                        header = f">{accession}|{gene_name}|{protein_name}"
                    
                    f.write(f"{header}\n")
                    f.write(f"{sequence}\n")
                
                print(f"DEBUG: FASTA written with header: {header}")
                print(f"DEBUG: Sequence length: {len(sequence)} amino acids")
            else:
                print(f"DEBUG: No sequence found for {first_target_name}")
                with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
                    f.write(f">unknown|{first_target_name}|No sequence available\n")
                    f.write("UNKNOWN\n")
        else:
            print("DEBUG: No sequences or targets found - writing empty FASTA")
            with open(os.path.join(args.output_dir, "protein.fasta"), "w") as f:
                f.write(">unknown|unknown|Unknown protein\n")
                f.write("UNKNOWN\n")
        
        # Write requirements.json - CORRECTED STRUCTURE
        molecule_properties = parsed_data.get("molecule_properties", {})
        requirements = {}
        
        # Extract requirements from molecule_properties if available
        if molecule_properties:
            for key, value in molecule_properties.items():
                if value is not None:
                    requirements[key] = value
        
        # Add default binding affinity
        if not requirements:
            requirements = {"binding_affinity": "< -8 kcal/mol"}
        elif "binding_affinity" not in requirements:
            requirements["binding_affinity"] = "< -8 kcal/mol"
        
        requirements["success"] = len(protein_targets) > 0
        
        with open(os.path.join(args.output_dir, "requirements.json"), "w") as f:
            json.dump(requirements, f, indent=2)
        
        # Write full result for debugging
        with open(os.path.join(args.output_dir, "target_selector_full_result.json"), "w") as f:
            json.dump(result, f, indent=2)
            
        logger.info("Target selection completed successfully")
        
    except ImportError as e:
        print(f"ERROR: Import failed: {str(e)}", file=sys.stderr)
        print("This might be due to missing dependencies in the conda environment", file=sys.stderr)
        # Write error outputs
        error_info = {"error": f"Import error: {str(e)}", "target": "unknown"}
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump(error_info, f, indent=2)
        sys.exit(1)
        
    except Exception as e:
        print(f"ERROR: Target selection failed: {str(e)}", file=sys.stderr)
        # Write error information
        error_info = {"error": str(e), "target": "unknown"}
        with open(os.path.join(args.output_dir, "target_info.json"), "w") as f:
            json.dump(error_info, f, indent=2)
        sys.exit(1)

if __name__ == "__main__":
    main()
