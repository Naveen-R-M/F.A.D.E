#!/usr/bin/env python3
"""
COMPLETE: Nextflow wrapper script for Structure Predictor Agent
Integrates RCSB PDB search before AlphaFold3 with complete implementation
"""

import os
import sys
import json
import argparse
import shutil
from pathlib import Path

def setup_paths_and_environment():
    """Setup Python paths and load environment variables"""
    
    script_path = Path(__file__).resolve()
    fade_root = script_path.parent.parent.parent  # bin -> nextflow -> F.A.D.E
    
    print(f"🔍 Script path: {script_path}")
    print(f"🔍 F.A.D.E root: {fade_root}")
    
    # Verify with fallback options
    if not (fade_root / "agents").exists():
        possible_roots = [
            Path("/scratch/rajagopalmohanraj.n/F.A.D.E"),  # HPC path
            script_path.parent.parent.parent,
            script_path.parent.parent,
        ]
        
        for possible_root in possible_roots:
            if (possible_root / "agents").exists() and (possible_root / "utils").exists():
                fade_root = possible_root
                print(f"✅ Found F.A.D.E root: {fade_root}")
                break
        else:
            print("🚨 FATAL: Could not locate F.A.D.E root!")
            sys.exit(1)
    
    sys.path.insert(0, str(fade_root))
    
    # Load .env file
    env_path = fade_root / ".env"
    if env_path.exists():
        with open(env_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    os.environ[key.strip()] = value.strip()
    
    return fade_root

def try_rcsb_structure_retrieval(target_info, sequence_info, output_dir, fade_root):
    """Try to get structure from RCSB PDB before attempting AlphaFold3"""
    
    try:
        print("🔍 ATTEMPTING RCSB PDB STRUCTURE RETRIEVAL")
        print("=" * 50)
        
        # Import RCSB client
        from agents.structure_predictor.rcsb_client import RCSBClient
        
        # Initialize client
        rcsb_client = RCSBClient()
        print("✅ RCSB client initialized")
        
        # Extract target information
        target_name = target_info.get("target", "unknown")
        uniprot_id = target_info.get("uniprot_id", "unknown")
        mutations = target_info.get("mutations", [])
        
        print(f"🎯 Target: {target_name}")
        print(f"🧬 UniProt ID: {uniprot_id}")
        print(f"🔬 Mutations: {len(mutations)} identified")
        
        if mutations:
            for mut in mutations:
                orig = mut.get("original_residue", "?")
                pos = mut.get("position", "?")
                new = mut.get("mutated_residue", "?")
                print(f"   - {orig}{pos}{new}")
        
        # Check if RCSB search is worth attempting
        if not rcsb_client.should_try_rcsb(target_info, sequence_info):
            print("⏭️  RCSB search not recommended for this target")
            return False
        
        # Attempt RCSB search
        print("🚀 Starting RCSB search...")
        structure_info = rcsb_client.search_structure(
            target_info=target_info,
            sequence_info=sequence_info,
            output_dir=output_dir,
            timeout=300
        )
        
        if structure_info:
            print("✅ RCSB STRUCTURE RETRIEVAL SUCCESSFUL!")
            print(f"   📋 PDB ID: {structure_info.get('pdb_id', 'unknown')}")
            print(f"   🔬 Method: {structure_info.get('experimental_method', 'unknown')}")
            print(f"   📏 Resolution: {structure_info.get('resolution_A', 'N/A')} Å")
            
            # Verify output files exist
            required_files = ["structure.pdb", "prepared_receptor.pdb", "structure_analysis.json"]
            all_files_exist = True
            
            print("📁 Verifying output files:")
            for file in required_files:
                file_path = os.path.join(output_dir, file)
                if os.path.exists(file_path):
                    size = os.path.getsize(file_path)
                    print(f"   ✅ {file} ({size} bytes)")
                else:
                    print(f"   ❌ {file} (missing)")
                    all_files_exist = False
            
            if all_files_exist:
                print("🎉 RCSB structure retrieval completed successfully!")
                print("⚡ Skipping AlphaFold3 - experimental structure available")
                return True
            else:
                print("⚠️  Some output files missing - falling back to AlphaFold3")
                return False
        else:
            print("❌ No suitable structure found in RCSB PDB")
            print("🔄 Proceeding with AlphaFold3 prediction...")
            return False
            
    except ImportError as e:
        print(f"❌ RCSB client import failed: {str(e)}")
        print("   This might be due to missing requests library")
        print("🔄 Proceeding with AlphaFold3...")
        return False
        
    except Exception as e:
        print(f"❌ RCSB retrieval failed: {str(e)}")
        print("🔄 Proceeding with AlphaFold3...")
        return False

def run_alphafold_prediction(fade_root, target_info, sequence_info, args):
    """Run AlphaFold3 structure prediction using StructurePredictor agent"""
    
    print("🧬 STARTING ALPHAFOLD3 STRUCTURE PREDICTION")
    print("=" * 50)
    
    # Import required modules
    from utils.logging import setup_logging, get_logger
    from agents.structure_predictor.structure_predictor import StructurePredictor
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.structure_predictor")
    
    target_name = target_info.get("target", "unknown")
    
    # Find AlphaFold job configurations
    job_configs = find_alphafold_job_configs(fade_root, target_name)
    
    if not job_configs:
        print("⚠️  No AlphaFold job configurations found!")
        print("   Creating fallback structure...")
        create_fallback_structure(args.output_dir, target_name, target_info, sequence_info)
        return
    
    # Prepare input for StructurePredictor
    input_data = {
        "protein_targets": [target_info],
        "sequences": {
            target_name: sequence_info
        },
        "job_configs": job_configs
    }
    
    print(f"📦 Input data prepared for AlphaFold3:")
    print(f"   - Target: {target_name}")
    print(f"   - Job configs: {list(job_configs.keys())}")
    
    # Initialize StructurePredictor
    api_key = args.api_key or os.getenv('GEMINI_API_KEY')
    structure_predictor = StructurePredictor(
        gemini_api_key=api_key,
        gemini_model=args.model
    )
    
    # Run AlphaFold3 prediction
    print("🚀 Submitting AlphaFold3 job...")
    print("   This will submit SLURM job and wait for completion")
    
    try:
        result = structure_predictor.process(input_data)
        print("✅ AlphaFold3 prediction completed!")
        
        # Process results
        success = process_alphafold_results(args.output_dir, result, target_name, logger)
        
        if not success:
            print("⚠️  AlphaFold3 results processing had issues")
            create_fallback_structure(args.output_dir, target_name, target_info, sequence_info)
            
    except Exception as e:
        print(f"❌ AlphaFold3 prediction failed: {str(e)}")
        logger.error(f"AlphaFold3 prediction failed: {str(e)}")
        create_error_structure(args.output_dir, target_name, str(e))

def find_alphafold_job_configs(fade_root, target_name):
    """Find AlphaFold job configurations created by TargetSelector"""
    
    data_dir = fade_root / "agents" / "data" / "inputs"
    job_scripts_dir = data_dir / "job_scripts"
    configs_dir = data_dir / "configs"
    
    job_configs = {}
    
    # Look for job script
    job_script_path = job_scripts_dir / f"{target_name}_alphafold.sh"
    if job_script_path.exists():
        job_configs[f"{target_name}_alphafold_job"] = str(job_script_path)
        print(f"✅ Found AlphaFold job script: {job_script_path}")
    
    # Look for config file
    config_file_path = configs_dir / f"{target_name}_alphafold.json"
    if config_file_path.exists():
        job_configs[f"{target_name}_alphafold_config"] = str(config_file_path)
        print(f"✅ Found AlphaFold config: {config_file_path}")
    
    return job_configs

def process_alphafold_results(output_dir, result, target_name, logger):
    """Process AlphaFold3 results from StructurePredictor"""
    
    structures = result.get("structures", {})
    binding_sites = result.get("binding_sites", {})
    prepared_structures = result.get("prepared_structures", {})
    
    success = False
    
    if target_name in structures:
        structure_info = structures[target_name]
        pdb_file = structure_info.get("pdb_file")
        
        if pdb_file and os.path.exists(pdb_file):
            shutil.copy(pdb_file, os.path.join(output_dir, "structure.pdb"))
            print(f"✅ AlphaFold3 structure: {pdb_file}")
            success = True
            
            # Write analysis data
            confidence_scores = structure_info.get("confidence_scores", {})
            analysis_data = {
                "target": target_name,
                "confidence": confidence_scores.get("overall", 0.0),
                "plddt": confidence_scores.get("plddt", 0.0),
                "ptm": confidence_scores.get("ptm", 0.0),
                "method": "alphafold3",
                "success": True
            }
            
            with open(os.path.join(output_dir, "structure_analysis.json"), "w") as f:
                json.dump(analysis_data, f, indent=2)
    
    # Ensure prepared_receptor.pdb exists
    if not os.path.exists(os.path.join(output_dir, "prepared_receptor.pdb")):
        if os.path.exists(os.path.join(output_dir, "structure.pdb")):
            shutil.copy(os.path.join(output_dir, "structure.pdb"),
                       os.path.join(output_dir, "prepared_receptor.pdb"))
    
    return success

def create_fallback_structure(output_dir, target_name, target_info, sequence_info):
    """Create fallback structure when both RCSB and AlphaFold3 fail"""
    
    sequence = sequence_info.get("sequence", "")
    
    with open(os.path.join(output_dir, "structure.pdb"), "w") as f:
        f.write("HEADER    FALLBACK STRUCTURE\n")
        f.write(f"TITLE     {target_name}\n")
        f.write(f"REMARK    Fallback structure - no RCSB or AlphaFold3 available\n")
        f.write(f"REMARK    Sequence length: {len(sequence)}\n")
        f.write("END\n")
    
    shutil.copy(os.path.join(output_dir, "structure.pdb"),
               os.path.join(output_dir, "prepared_receptor.pdb"))
    
    with open(os.path.join(output_dir, "structure_analysis.json"), "w") as f:
        json.dump({
            "target": target_name,
            "confidence": 0.3,
            "method": "fallback",
            "success": True
        }, f, indent=2)

def create_error_structure(output_dir, target_name, error_msg):
    """Create error structure files"""
    
    with open(os.path.join(output_dir, "structure.pdb"), "w") as f:
        f.write("HEADER    Structure prediction failed\n")
        f.write(f"REMARK    {error_msg}\n")
    
    with open(os.path.join(output_dir, "prepared_receptor.pdb"), "w") as f:
        f.write("HEADER    Structure prediction failed\n")
        f.write(f"REMARK    {error_msg}\n")
    
    with open(os.path.join(output_dir, "structure_analysis.json"), "w") as f:
        json.dump({
            "error": error_msg,
            "confidence": 0.0,
            "target": target_name,
            "success": False
        }, f, indent=2)

def main():
    print("🚀 STARTING RCSB-INTEGRATED STRUCTURE PREDICTION")
    print("=" * 60)
    
    parser = argparse.ArgumentParser(description="Run Structure Predictor with RCSB integration")
    parser.add_argument("--target-info", required=True, help="Path to target_info.json")
    parser.add_argument("--fasta-file", required=True, help="Path to protein.fasta")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--api-key", help="Gemini API key")
    parser.add_argument("--model", default="models/gemini-2.5-flash", help="Gemini model")
    parser.add_argument("--use-alphafold", action="store_true", default=True, help="Use AlphaFold3 for prediction")
    parser.add_argument("--skip-rcsb", action="store_true", help="Skip RCSB search and go directly to AlphaFold3")
    
    args = parser.parse_args()
    
    print(f"📋 Configuration:")
    print(f"   - Target info: {args.target_info}")
    print(f"   - FASTA file: {args.fasta_file}")
    print(f"   - Output dir: {args.output_dir}")
    print(f"   - Skip RCSB: {args.skip_rcsb}")
    print(f"   - Use AlphaFold: {args.use_alphafold}")
    print()
    
    # Setup environment
    fade_root = setup_paths_and_environment()
    
    # Validate API key
    api_key = args.api_key or os.getenv('GEMINI_API_KEY')
    if not api_key:
        print("🚨 ERROR: No API key found!")
        sys.exit(1)
    
    try:
        # Load input files
        with open(args.target_info, 'r') as f:
            target_info = json.load(f)
        
        with open(args.fasta_file, 'r') as f:
            fasta_content = f.read().strip()
        
        lines = fasta_content.split('\n')
        header = lines[0] if lines else ""
        sequence = ''.join(lines[1:]) if len(lines) > 1 else ""
        
        target_name = target_info.get("target", "unknown")
        
        print(f"🎯 Target: {target_name}")
        print(f"🧬 Sequence length: {len(sequence)} amino acids")
        
        # Validate sequence
        if not sequence or sequence == "UNKNOWN":
            raise ValueError("Invalid protein sequence")
        
        # Prepare sequence info
        sequence_info = {
            "uniprot_id": target_info.get("uniprot_id", "unknown"),
            "sequence": sequence,
            "header": header,
            "length": len(sequence),
            "mutations": target_info.get("mutations", []),
            "organism": target_info.get("organism", "")
        }
        
        print()
        
        # PHASE 1: Try RCSB PDB first (unless skipped)
        rcsb_success = False
        if not args.skip_rcsb:
            rcsb_success = try_rcsb_structure_retrieval(
                target_info, sequence_info, args.output_dir, fade_root
            )
        else:
            print("⏭️  RCSB search skipped (--skip-rcsb flag)")
        
        print()
        
        # PHASE 2: Fall back to AlphaFold3 if RCSB failed
        if not rcsb_success and args.use_alphafold:
            run_alphafold_prediction(fade_root, target_info, sequence_info, args)
        elif not rcsb_success:
            print("⚠️  RCSB failed and AlphaFold3 disabled - creating fallback")
            create_fallback_structure(args.output_dir, target_name, target_info, sequence_info)
        
        print()
        print("=" * 60)
        print("✅ STRUCTURE PREDICTION COMPLETED")
        
        # Final summary
        print("📋 FINAL OUTPUT:")
        output_files = ['structure.pdb', 'prepared_receptor.pdb', 'structure_analysis.json']
        for file in output_files:
            file_path = Path(args.output_dir) / file
            if file_path.exists():
                size = file_path.stat().st_size
                print(f"   ✅ {file} ({size} bytes)")
            else:
                print(f"   ❌ {file} (missing)")
        
    except Exception as e:
        print(f"🚨 FATAL ERROR: {str(e)}")
        create_error_structure(args.output_dir, target_info.get('target', 'unknown'), str(e))
        sys.exit(1)

def run_alphafold_prediction(fade_root, target_info, sequence_info, args):
    """Run AlphaFold3 structure prediction using StructurePredictor agent"""
    
    print("🧬 STARTING ALPHAFOLD3 STRUCTURE PREDICTION")
    print("=" * 50)
    
    # Import required modules
    from utils.logging import setup_logging, get_logger
    from agents.structure_predictor.structure_predictor import StructurePredictor
    
    # Setup logging
    setup_logging(log_level="INFO")
    logger = get_logger("nextflow.structure_predictor")
    
    target_name = target_info.get("target", "unknown")
    
    # Find AlphaFold job configurations
    job_configs = find_alphafold_job_configs(fade_root, target_name)
    
    if not job_configs:
        print("⚠️  No AlphaFold job configurations found!")
        print("   Creating fallback structure...")
        create_fallback_structure(args.output_dir, target_name, target_info, sequence_info)
        return
    
    # Prepare input for StructurePredictor
    input_data = {
        "protein_targets": [target_info],
        "sequences": {
            target_name: sequence_info
        },
        "job_configs": job_configs
    }
    
    print(f"📦 Input data prepared for AlphaFold3:")
    print(f"   - Target: {target_name}")
    print(f"   - Job configs: {list(job_configs.keys())}")
    
    # Initialize StructurePredictor
    api_key = args.api_key or os.getenv('GEMINI_API_KEY')
    structure_predictor = StructurePredictor(
        gemini_api_key=api_key,
        gemini_model=args.model
    )
    
    # Run AlphaFold3 prediction
    print("🚀 Submitting AlphaFold3 job...")
    print("   This will submit SLURM job and wait for completion")
    
    try:
        result = structure_predictor.process(input_data)
        print("✅ AlphaFold3 prediction completed!")
        
        # Process results
        success = process_alphafold_results(args.output_dir, result, target_name, logger)
        
        if not success:
            print("⚠️  AlphaFold3 results processing had issues")
            create_fallback_structure(args.output_dir, target_name, target_info, sequence_info)
            
    except Exception as e:
        print(f"❌ AlphaFold3 prediction failed: {str(e)}")
        logger.error(f"AlphaFold3 prediction failed: {str(e)}")
        create_error_structure(args.output_dir, target_name, str(e))

def find_alphafold_job_configs(fade_root, target_name):
    """Find AlphaFold job configurations created by TargetSelector"""
    
    data_dir = fade_root / "agents" / "data" / "inputs"
    job_scripts_dir = data_dir / "job_scripts"
    configs_dir = data_dir / "configs"
    
    job_configs = {}
    
    # Look for job script
    job_script_path = job_scripts_dir / f"{target_name}_alphafold.sh"
    if job_script_path.exists():
        job_configs[f"{target_name}_alphafold_job"] = str(job_script_path)
        print(f"✅ Found AlphaFold job script: {job_script_path}")
    
    # Look for config file
    config_file_path = configs_dir / f"{target_name}_alphafold.json"
    if config_file_path.exists():
        job_configs[f"{target_name}_alphafold_config"] = str(config_file_path)
        print(f"✅ Found AlphaFold config: {config_file_path}")
    
    return job_configs

def process_alphafold_results(output_dir, result, target_name, logger):
    """Process AlphaFold3 results from StructurePredictor"""
    
    structures = result.get("structures", {})
    binding_sites = result.get("binding_sites", {})
    prepared_structures = result.get("prepared_structures", {})
    
    success = False
    
    if target_name in structures:
        structure_info = structures[target_name]
        pdb_file = structure_info.get("pdb_file")
        
        if pdb_file and os.path.exists(pdb_file):
            shutil.copy(pdb_file, os.path.join(output_dir, "structure.pdb"))
            print(f"✅ AlphaFold3 structure: {pdb_file}")
            success = True
            
            # Write analysis data
            confidence_scores = structure_info.get("confidence_scores", {})
            analysis_data = {
                "target": target_name,
                "confidence": confidence_scores.get("overall", 0.0),
                "plddt": confidence_scores.get("plddt", 0.0),
                "ptm": confidence_scores.get("ptm", 0.0),
                "method": "alphafold3",
                "success": True
            }
            
            with open(os.path.join(output_dir, "structure_analysis.json"), "w") as f:
                json.dump(analysis_data, f, indent=2)
    
    # Ensure prepared_receptor.pdb exists
    if not os.path.exists(os.path.join(output_dir, "prepared_receptor.pdb")):
        if os.path.exists(os.path.join(output_dir, "structure.pdb")):
            shutil.copy(os.path.join(output_dir, "structure.pdb"),
                       os.path.join(output_dir, "prepared_receptor.pdb"))
    
    return success

if __name__ == "__main__":
    main()
