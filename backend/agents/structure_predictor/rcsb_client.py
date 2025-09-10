"""
RCSB PDB Client for F.A.D.E Structure Predictor
Provides integration between RCSB PDB search and F.A.D.E structure prediction workflow
"""

import os
import json
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any
from utils.logging import get_logger

class RCSBClient:
    """
    Client for retrieving protein structures from RCSB PDB
    Integrates with the F.A.D.E structure prediction workflow
    """
    
    def __init__(self, script_path: Optional[str] = None):
        """
        Initialize RCSB client
        
        Args:
            script_path: Path to get_structures_from_rcsb.py script
        """
        self.logger = get_logger("fade.rcsb_client")
        
        # Find the RCSB script
        if script_path and Path(script_path).exists():
            self.script_path = Path(script_path)
        else:
            # Look in the same directory as this file
            current_dir = Path(__file__).parent
            self.script_path = current_dir / "get_structures_from_rcsb.py"
            
            if not self.script_path.exists():
                raise FileNotFoundError(f"RCSB script not found at {self.script_path}")
        
        self.logger.info(f"RCSB Client initialized with script: {self.script_path}")
    
    def should_try_rcsb(self, target_info: Dict[str, Any], sequence_info: Dict[str, Any]) -> bool:
        """
        Determine if RCSB lookup is worth attempting based on target information
        
        Args:
            target_info: Target information from TargetSelector
            sequence_info: Sequence information
            
        Returns:
            Boolean indicating if RCSB search should be attempted
        """
        
        # High priority cases for RCSB lookup
        uniprot_id = target_info.get("uniprot_id") or sequence_info.get("uniprot_id")
        target_name = target_info.get("target", "").upper()
        organism = sequence_info.get("organism", "")
        
        # Skip RCSB for obviously invalid cases
        if not uniprot_id and not target_name:
            self.logger.info("No UniProt ID or target name - skipping RCSB")
            return False
        
        if uniprot_id == "unknown" or target_name == "UNKNOWN":
            self.logger.info("Unknown target/UniProt - skipping RCSB")
            return False
        
        # Check for well-known proteins that likely have structures
        well_known_targets = [
            "KRAS", "EGFR", "TP53", "BRAF", "PIK3CA", "PTEN", "MYC", "RAS",
            "BRCA1", "BRCA2", "ATM", "CHEK2", "MDM2", "BCL2", "VEGFA"
        ]
        
        if any(known in target_name for known in well_known_targets):
            self.logger.info(f"Well-known target {target_name} - attempting RCSB")
            return True
        
        # Always try if we have a UniProt ID
        if uniprot_id and uniprot_id != "unknown":
            self.logger.info(f"Valid UniProt ID {uniprot_id} - attempting RCSB")
            return True
        
        # Try for human/mouse proteins
        if organism and any(org in organism.lower() for org in ["homo sapiens", "human", "mouse", "mus musculus"]):
            self.logger.info(f"Model organism {organism} - attempting RCSB")
            return True
        
        self.logger.info("Target doesn't meet RCSB criteria - skipping")
        return False
    
    def search_structure(
        self,
        target_info: Dict[str, Any],
        sequence_info: Dict[str, Any],
        output_dir: str,
        timeout: int = 300
    ) -> Optional[Dict[str, Any]]:
        """
        Search for protein structure in RCSB PDB
        
        Args:
            target_info: Target information from TargetSelector
            sequence_info: Sequence information
            output_dir: Directory to save structure files
            timeout: Timeout for RCSB search in seconds
            
        Returns:
            Structure information dict if successful, None otherwise
        """
        
        if not self.should_try_rcsb(target_info, sequence_info):
            return None
        
        target_name = target_info.get("target", "PROTEIN")
        uniprot_id = target_info.get("uniprot_id") or sequence_info.get("uniprot_id")
        organism = sequence_info.get("organism")
        
        self.logger.info(f"Searching RCSB for {target_name} (UniProt: {uniprot_id})")
        
        # Build command for RCSB script
        cmd = ["python", str(self.script_path)]
        
        # Add search parameters in order of preference
        if uniprot_id and uniprot_id != "unknown":
            cmd.extend(["--uniprot", uniprot_id])
            self.logger.info(f"Using UniProt ID: {uniprot_id}")
        elif target_name and target_name != "UNKNOWN":
            cmd.extend(["--gene", target_name])
            self.logger.info(f"Using gene name: {target_name}")
        else:
            self.logger.warning("No valid search parameters for RCSB")
            return None
        
        # Add organism if available
        if organism and organism not in ["unknown", ""]:
            cmd.extend(["--organism", organism])
            self.logger.info(f"Using organism: {organism}")
        
        # Add output parameters
        cmd.extend([
            "--protein-tag", target_name,
            "--output-dir", output_dir,
            "--limit", "50"  # Reasonable limit for speed
        ])
        
        # Handle mutations - for now, get wild-type structure
        mutations = target_info.get("mutations", [])
        if mutations:
            self.logger.info(f"Target has {len(mutations)} mutations - searching for wild-type structure")
            # Could add mutation-specific search logic here in the future
        
        try:
            self.logger.info(f"Executing RCSB search: {' '.join(cmd[2:])}")  # Don't log full python path
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=Path(self.script_path).parent  # Run from script directory
            )
            
            if result.returncode == 0:
                self.logger.info("RCSB search completed successfully")
                
                # Parse the output and convert to F.A.D.E format
                structure_info = self._convert_rcsb_output_to_fade_format(
                    output_dir, target_name, result.stdout, target_info, mutations
                )
                
                if structure_info:
                    self.logger.info(f"Successfully retrieved structure: {structure_info.get('pdb_id', 'unknown')}")
                    return structure_info
                else:
                    self.logger.warning("RCSB search completed but no usable structure found")
                    return None
            else:
                self.logger.warning(f"RCSB search failed: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            self.logger.warning(f"RCSB search timed out after {timeout} seconds")
            return None
        except Exception as e:
            self.logger.error(f"RCSB search error: {str(e)}")
            return None
    
    def _convert_rcsb_output_to_fade_format(
        self,
        output_dir: str,
        target_name: str,
        rcsb_output: str,
        target_info: Dict[str, Any],
        mutations: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """
        Convert RCSB output to F.A.D.E StructurePredictor format
        
        Args:
            output_dir: Output directory where RCSB files were saved
            target_name: Target protein name
            rcsb_output: Output from RCSB script
            target_info: Original target information
            mutations: List of mutations in the target
            
        Returns:
            Structure information in F.A.D.E format
        """
        
        try:
            # Find the files created by RCSB script
            best_model_cif = Path(output_dir) / f"{target_name}_best_model.cif"
            summary_json_path = Path(output_dir) / f"{target_name}_af3_input" / target_name.lower() / f"{target_name.lower()}_summary_confidences.json"
            
            if not best_model_cif.exists():
                self.logger.error(f"RCSB structure file not found: {best_model_cif}")
                return None
            
            # Convert CIF to PDB format for compatibility
            pdb_file = self._convert_cif_to_pdb(best_model_cif, output_dir)
            
            # Load RCSB metadata
            rcsb_metadata = {}
            if summary_json_path.exists():
                with open(summary_json_path, 'r') as f:
                    rcsb_metadata = json.load(f)
            
            # Create F.A.D.E-compatible structure info
            structure_info = {
                "target_name": target_name,
                "pdb_file": pdb_file,
                "source": "rcsb_pdb",
                "pdb_id": rcsb_metadata.get("pdb_id", "unknown"),
                "experimental_method": rcsb_metadata.get("experimental_method", "unknown"),
                "resolution_A": rcsb_metadata.get("resolution_A"),
                "confidence_scores": {
                    "overall": 1.0,  # High confidence for experimental structures
                    "plddt": 95.0,   # Mock high pLDDT for experimental
                    "ptm": 0.9,      # Mock PTM score
                    "iptm": None
                },
                "validation_results": {
                    "overall_score": 1.0,
                    "method": "experimental",
                    "source": "rcsb_pdb"
                },
                "chain_info": {},  # Would need PDB parsing to fill this
                "residue_count": 0,  # Would need PDB parsing
                "atom_count": 0,     # Would need PDB parsing
                "secondary_structure": {},
                "mutations": mutations,  # Pass through mutations for annotation
                "structure_type": "experimental_wildtype" if mutations else "experimental"
            }
            
            # Create structure analysis JSON for Nextflow compatibility
            analysis_data = {
                "target": target_name,
                "confidence": 1.0,
                "plddt": 95.0,
                "ptm": 0.9,
                "method": "experimental_rcsb",
                "pdb_id": rcsb_metadata.get("pdb_id", "unknown"),
                "resolution_A": rcsb_metadata.get("resolution_A"),
                "experimental_method": rcsb_metadata.get("experimental_method", "unknown"),
                "success": True,
                "mutations": mutations,
                "structure_type": "experimental_wildtype" if mutations else "experimental"
            }
            
            # Write structure analysis for Nextflow
            with open(os.path.join(output_dir, "structure_analysis.json"), "w") as f:
                json.dump(analysis_data, f, indent=2)
            
            # Create prepared receptor (copy of structure for now)
            prepared_receptor_path = os.path.join(output_dir, "prepared_receptor.pdb")
            if pdb_file and os.path.exists(pdb_file):
                shutil.copy(pdb_file, prepared_receptor_path)
            
            self.logger.info(f"Successfully converted RCSB structure for {target_name}")
            return structure_info
            
        except Exception as e:
            self.logger.error(f"Error converting RCSB output: {str(e)}")
            return None
    
    def _convert_cif_to_pdb(self, cif_path: Path, output_dir: str) -> str:
        """
        Convert CIF file to PDB format
        
        Args:
            cif_path: Path to CIF file
            output_dir: Output directory
            
        Returns:
            Path to converted PDB file
        """
        
        pdb_path = os.path.join(output_dir, "structure.pdb")
        
        try:
            # Try using BioPython if available
            try:
                from Bio import PDB
                
                # Parse CIF
                parser = PDB.MMCIFParser(QUIET=True)
                structure = parser.get_structure("structure", str(cif_path))
                
                # Write PDB
                io = PDB.PDBIO()
                io.set_structure(structure)
                io.save(pdb_path)
                
                self.logger.info(f"Converted CIF to PDB using BioPython: {pdb_path}")
                return pdb_path
                
            except ImportError:
                self.logger.warning("BioPython not available for CIF conversion")
                pass
            
            # Fallback: simple copy and rename (most PDB readers can handle CIF)
            shutil.copy(str(cif_path), pdb_path)
            self.logger.info(f"Copied CIF as PDB (fallback): {pdb_path}")
            return pdb_path
            
        except Exception as e:
            self.logger.error(f"CIF conversion failed: {str(e)}")
            # Last resort: just return the CIF path
            return str(cif_path)
