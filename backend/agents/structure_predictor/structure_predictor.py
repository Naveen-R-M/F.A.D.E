"""
Structure Predictor Agent for F.A.D.E - Enhanced with RCSB Integration

This agent handles the generation and validation of 3D protein structures
from sequence data, and identifies potential binding sites for drug targeting.

ENHANCED: Now tries RCSB PDB first before AlphaFold3 for faster results.
"""

import os
import json
import time
from typing import Any, Dict, List, Optional, Tuple, Union

from agents.base.base_agent import BaseAgent
from agents.base.agentic_mixin import AgenticMixin
from agents.structure_predictor.pdb_processor import PDBProcessor
from agents.structure_predictor.structure_validator import StructureValidator
from agents.structure_predictor.binding_site_detector import BindingSiteDetector
from utils.gemini_client import GeminiClient
from utils.slurm_client import SlurmClient
from utils.alphafold_client import AlphaFoldClient
from utils.config_generator import ConfigGenerator
from utils.logging import get_logger

# NEW: Import RCSB client
try:
    from agents.structure_predictor.rcsb_client import RCSBClient
except ImportError:
    RCSBClient = None


class StructurePredictor(BaseAgent, AgenticMixin):
    """
    Agent for predicting, validating, and analyzing protein structures.
    
    ENHANCED WORKFLOW:
    1. First tries RCSB PDB for existing experimental structures (fast)
    2. Falls back to AlphaFold3 prediction if no suitable structure found
    3. Validates and processes structures
    4. Identifies potential binding sites
    5. Prepares structures for docking
    """
    
    def __init__(
        self,
        name: str = "structure_predictor",
        config: Optional[Dict[str, Any]] = None,
        gemini_api_key: Optional[str] = None,
        gemini_model: Optional[str] = None,
        enable_rcsb: bool = True,
    ) -> None:
        """
        Initialize the Structure Predictor agent with RCSB integration.
        
        Args:
            name: Unique identifier for the agent.
            config: Optional configuration parameters.
            gemini_api_key: API key for the Gemini model.
            gemini_model: Model name for Gemini.
            enable_rcsb: Whether to enable RCSB PDB lookup (default: True)
        """
        # Initialize base agent
        BaseAgent.__init__(self, name, config)
        
        # Initialize clients
        self.gemini_client = GeminiClient(api_key=gemini_api_key, model=gemini_model)
        self.slurm_client = SlurmClient()
        self.alphafold_client = AlphaFoldClient()
        self.config_generator = ConfigGenerator()
        
        # NEW: Initialize RCSB client if available and enabled
        self.enable_rcsb = enable_rcsb and RCSBClient is not None
        if self.enable_rcsb:
            try:
                self.rcsb_client = RCSBClient()
                self.logger.info("RCSB client initialized - will try PDB first")
            except Exception as e:
                self.logger.warning(f"RCSB client initialization failed: {e}")
                self.enable_rcsb = False
                self.rcsb_client = None
        else:
            self.rcsb_client = None
            if RCSBClient is None:
                self.logger.info("RCSB client not available - using AlphaFold3 only")
            else:
                self.logger.info("RCSB client disabled - using AlphaFold3 only")
        
        # Initialize agentic components
        AgenticMixin.initialize_agentic_components(self, llm_client=self.gemini_client)
        
        # Create specialized components
        self.pdb_processor = PDBProcessor()
        self.structure_validator = StructureValidator(llm_client=self.gemini_client)
        self.binding_site_detector = BindingSiteDetector(llm_client=self.gemini_client)
        
        # Set up memory file
        data_dir = self._get_data_dir()
        memory_dir = os.path.join(data_dir, "memory")
        os.makedirs(memory_dir, exist_ok=True)
        self.memory_file = os.path.join(memory_dir, f"{name}_memory.json")
    
    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process protein sequence data to generate and analyze structures.
        
        ENHANCED: Now tries RCSB PDB first, then AlphaFold3 if needed.
        
        Args:
            input_data: Dictionary containing:
                - protein_targets: List of protein target information
                - sequences: Dictionary mapping target names to sequence info
                - job_configs: Dictionary of job configuration file paths
                
        Returns:
            Dictionary containing:
                - structures: Dictionary mapping target names to structure info
                - binding_sites: Dictionary mapping target names to binding site info
                - prepared_structures: Dictionary mapping target names to prepared structure paths
        """
        self.logger.info("Processing structure prediction for %d targets", 
                         len(input_data.get("sequences", {})))
        
        sequences = input_data.get("sequences", {})
        job_configs = input_data.get("job_configs", {})
        protein_targets = input_data.get("protein_targets", [])
        
        structures = {}
        binding_sites = {}
        prepared_structures = {}
        
        for target_name, sequence_info in sequences.items():
            self.logger.info(f"Processing target: {target_name}")
            
            # Get corresponding target info
            target_info = None
            for target in protein_targets:
                if target.get("target") == target_name:
                    target_info = target
                    break
            
            if not target_info:
                self.logger.warning(f"No target info found for {target_name}")
                target_info = {"target": target_name}
            
            # ENHANCED: Try RCSB first, then AlphaFold3
            structure_info = self.get_structure_rcsb_or_alphafold(
                target_name, target_info, sequence_info, job_configs
            )
            
            if structure_info:
                structures[target_name] = structure_info
                
                # Identify binding sites
                self.logger.info(f"Identifying binding sites for {target_name}")
                target_binding_sites = self.identify_binding_sites(
                    target_name, structure_info, sequence_info
                )
                binding_sites[target_name] = target_binding_sites
                
                # Prepare structure for docking
                self.logger.info(f"Preparing structure for docking: {target_name}")
                prepared_info = self.prepare_structure_for_docking(
                    target_name, structure_info, target_binding_sites
                )
                if prepared_info:
                    prepared_structures[target_name] = prepared_info
            else:
                self.logger.warning(f"Failed to get structure for {target_name}")
        
        return {
            "structures": structures,
            "binding_sites": binding_sites,
            "prepared_structures": prepared_structures
        }
    
    def get_structure_rcsb_or_alphafold(
        self,
        target_name: str,
        target_info: Dict[str, Any],
        sequence_info: Dict[str, Any],
        job_configs: Dict[str, str]
    ) -> Optional[Dict[str, Any]]:
        """
        Get protein structure using RCSB first, then AlphaFold3 if needed.
        
        Args:
            target_name: Name of the target protein
            target_info: Target information from TargetSelector
            sequence_info: Sequence information
            job_configs: AlphaFold job configurations
            
        Returns:
            Structure information dict or None if both methods fail
        """
        
        # PHASE 1: Try RCSB PDB first
        if self.enable_rcsb and self.rcsb_client:
            self.logger.info(f"Attempting RCSB structure retrieval for {target_name}")
            
            try:
                data_dir = self._get_data_dir()
                rcsb_output_dir = os.path.join(data_dir, "outputs", "structures", "rcsb", target_name)
                os.makedirs(rcsb_output_dir, exist_ok=True)
                
                rcsb_structure = self.rcsb_client.search_structure(
                    target_info=target_info,
                    sequence_info=sequence_info,
                    output_dir=rcsb_output_dir,
                    timeout=300
                )
                
                if rcsb_structure:
                    self.logger.info(f"Successfully retrieved RCSB structure for {target_name}")
                    
                    # Enhance RCSB structure info with additional processing
                    enhanced_structure = self.enhance_rcsb_structure(rcsb_structure, sequence_info)
                    
                    # Log success and return
                    self.learn_from_interaction({
                        "interaction_type": "structure_acquisition",
                        "key": f"rcsb_success:{target_name}",
                        "target_name": target_name,
                        "method": "rcsb_pdb",
                        "pdb_id": rcsb_structure.get("pdb_id", "unknown"),
                        "success": True
                    })
                    
                    return enhanced_structure
                else:
                    self.logger.info(f"No suitable RCSB structure found for {target_name}")
                    
            except Exception as e:
                self.logger.warning(f"RCSB structure retrieval failed for {target_name}: {e}")
        
        # PHASE 2: Fall back to AlphaFold3
        self.logger.info(f"Attempting AlphaFold3 prediction for {target_name}")
        
        job_config_key = f"{target_name}_alphafold_job"
        job_script_path = job_configs.get(job_config_key)
        
        if not job_script_path:
            self.logger.warning(f"No AlphaFold job script found for {target_name}")
            return None
        
        alphafold_structure = self.predict_structure(target_name, sequence_info, job_script_path)
        
        if alphafold_structure:
            self.logger.info(f"Successfully predicted AlphaFold3 structure for {target_name}")
            
            # Log success
            self.learn_from_interaction({
                "interaction_type": "structure_acquisition",
                "key": f"alphafold_success:{target_name}",
                "target_name": target_name,
                "method": "alphafold3",
                "confidence": alphafold_structure.get("confidence_scores", {}).get("overall", 0.0),
                "success": True
            })
            
            return alphafold_structure
        else:
            self.logger.error(f"Both RCSB and AlphaFold3 failed for {target_name}")
            
            # Log failure
            self.learn_from_interaction({
                "interaction_type": "structure_acquisition",
                "key": f"structure_failure:{target_name}",
                "target_name": target_name,
                "method": "both_failed",
                "success": False
            })
            
            return None
    
    def enhance_rcsb_structure(
        self, 
        rcsb_structure: Dict[str, Any], 
        sequence_info: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Enhance RCSB structure information with additional processing
        
        Args:
            rcsb_structure: Structure info from RCSB client
            sequence_info: Original sequence information
            
        Returns:
            Enhanced structure information
        """
        
        try:
            pdb_file = rcsb_structure.get("pdb_file")
            if not pdb_file or not os.path.exists(pdb_file):
                self.logger.warning("RCSB structure file not found for enhancement")
                return rcsb_structure
            
            # Process PDB file to get detailed structure data
            structure_data = self.pdb_processor.parse_pdb(pdb_file)
            
            # Validate structure against sequence
            validation_results = self.structure_validator.validate(
                pdb_file, 
                sequence_info.get("sequence", "")
            )
            
            # Enhance the structure info
            enhanced_structure = {
                **rcsb_structure,  # Keep all original RCSB info
                "chain_info": structure_data.get("chains", {}),
                "residue_count": structure_data.get("residue_count", 0),
                "atom_count": structure_data.get("atom_count", 0),
                "secondary_structure": structure_data.get("secondary_structure", {}),
                "validation_results": validation_results,
                "structure_source": "rcsb_pdb_enhanced"
            }
            
            self.logger.info(f"Enhanced RCSB structure: {enhanced_structure['residue_count']} residues, {enhanced_structure['atom_count']} atoms")
            return enhanced_structure
            
        except Exception as e:
            self.logger.warning(f"Failed to enhance RCSB structure: {e}")
            return rcsb_structure
    
    def predict_structure(
        self,
        target_name: str,
        sequence_info: Dict[str, Any],
        job_script_path: str
    ) -> Optional[Dict[str, Any]]:
        """
        Run AlphaFold to predict the structure of a protein.
        
        Args:
            target_name: Name of the target protein.
            sequence_info: Dictionary containing sequence information.
            job_script_path: Path to the SLURM job script.
            
        Returns:
            Dictionary containing structure information or None if prediction failed.
        """
        try:
            job_id = self.execute_with_retry(
                self.slurm_client.submit_job,
                job_script_path,
                operation_name=f"Submit AlphaFold job for {target_name}"
            )
            
            self.logger.info("Submitted AlphaFold job for %s (Job ID: %s)", 
                             target_name, job_id)
            
            # Monitor job
            job_status = self.slurm_client.monitor_job(job_id)
            
            if job_status != "COMPLETED":
                self.logger.error("AlphaFold job for %s failed with status: %s", 
                                 target_name, job_status)
                return None
                
            # Get output path
            data_dir = self._get_data_dir()
            alphafold_base_output_dir = os.path.join(data_dir, "outputs", "structures", target_name)
            
            # Get PDB file
            pdb_file = self.alphafold_client.get_best_model_path(alphafold_base_output_dir, target_name)
            
            if not pdb_file or not os.path.exists(pdb_file):
                self.logger.error("No PDB file found for %s", target_name)
                return None
                
            # Get confidence scores
            confidence_scores_af = self.alphafold_client.get_confidence_scores(alphafold_base_output_dir, target_name)
            multimer_info_af = self.alphafold_client.get_multimer_info(alphafold_base_output_dir, target_name)

            # Process PDB file
            structure_data = self.pdb_processor.parse_pdb(pdb_file)
            
            # Validate structure
            validation_results = self.structure_validator.validate(
                pdb_file, 
                sequence_info.get("sequence", "")
            )
            
            # Create structure info
            structure_info = {
                "target_name": target_name,
                "pdb_file": pdb_file,
                "source": "alphafold3",
                "confidence_scores": {
                    "overall": confidence_scores_af.get("ranking_score", validation_results.get("overall_score", 0.0)) if confidence_scores_af else validation_results.get("overall_score", 0.0),
                    "plddt": confidence_scores_af.get("avg_plddt", validation_results.get("plddt", 0.0)) if confidence_scores_af else validation_results.get("plddt", 0.0),
                    "ptm": confidence_scores_af.get("ptm", validation_results.get("ptm", 0.0)) if confidence_scores_af else validation_results.get("ptm", 0.0),
                    "iptm": multimer_info_af.get("iptm", None) if multimer_info_af else None
                },
                "validation_results": validation_results,
                "chain_info": structure_data.get("chains", {}),
                "residue_count": structure_data.get("residue_count", 0),
                "atom_count": structure_data.get("atom_count", 0),
                "secondary_structure": structure_data.get("secondary_structure", {}),
                "structure_source": "alphafold3_predicted"
            }
            
            self.logger.info("AlphaFold3 prediction for %s: %d residues, %d atoms, confidence: %.2f", 
                             target_name, 
                             structure_info["residue_count"],
                             structure_info["atom_count"],
                             structure_info["confidence_scores"]["overall"])
            
            return structure_info
            
        except Exception as e:
            self.logger.error("Error predicting structure for %s: %s", target_name, e)
            
            self.learn_from_interaction({
                "interaction_type": "structure_prediction",
                "key": f"prediction:{target_name}",
                "target_name": target_name,
                "error_message": str(e),
                "success": False
            })
            
            return None
    
    def identify_binding_sites(
        self,
        target_name: str,
        structure_info: Dict[str, Any],
        sequence_info: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Identify potential binding sites in the protein structure.
        
        Args:
            target_name: Name of the target protein.
            structure_info: Dictionary containing structure information.
            sequence_info: Dictionary containing sequence information.
            
        Returns:
            List of dictionaries describing binding sites.
        """
        try:
            pdb_file = structure_info.get("pdb_file")
            
            if not pdb_file:
                self.logger.warning("No PDB file found for %s", target_name)
                return []
                
            # Extract known binding sites from sequence info
            known_binding_sites = []
            for site in sequence_info.get("binding_sites", []):
                known_binding_sites.append({
                    "name": site.get("name", "Unknown"),
                    "residues": site.get("residues", []),
                    "type": site.get("type", "Unknown"),
                    "ligand": site.get("ligand", "Unknown")
                })
            
            # Identify binding sites using detector
            detected_sites = self.binding_site_detector.detect(
                pdb_file,
                known_binding_sites=known_binding_sites
            )
            
            # Process mutations as potential binding sites
            mutations = sequence_info.get("mutations", [])
            if mutations:
                for mutation in mutations:
                    position = mutation.get("position")
                    mutated = mutation.get("mutated_residue")
                    original = mutation.get("original_residue")
                    
                    if position and mutated and original:
                        # Check if mutation is within existing binding site
                        in_binding_site = any(
                            position in site.get("residue_ids", []) 
                            for site in detected_sites
                        )
                        
                        if not in_binding_site:
                            mutation_site = {
                                "name": f"Mutation_{original}{position}{mutated}",
                                "type": "mutation",
                                "center": {"x": 0, "y": 0, "z": 0},
                                "radius": 10.0,
                                "residue_ids": [position],
                                "score": 0.9,
                                "description": f"Mutation site {original}{position}{mutated}"
                            }
                            detected_sites.append(mutation_site)
            
            self.logger.info("Detected %d binding sites for %s", 
                             len(detected_sites), target_name)
            
            return detected_sites
            
        except Exception as e:
            self.logger.error("Error identifying binding sites for %s: %s", 
                             target_name, e)
            
            self.learn_from_interaction({
                "interaction_type": "binding_site_detection",
                "key": f"binding_site:{target_name}",
                "target_name": target_name,
                "error_message": str(e),
                "success": False
            })
            
            return []
    
    def prepare_structure_for_docking(
        self,
        target_name: str,
        structure_info: Dict[str, Any],
        binding_sites: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """
        Prepare the protein structure for docking.
        
        Args:
            target_name: Name of the target protein.
            structure_info: Dictionary containing structure information.
            binding_sites: List of dictionaries describing binding sites.
            
        Returns:
            Dictionary containing prepared structure information or None if preparation failed.
        """
        try:
            pdb_file = structure_info.get("pdb_file")
            
            if not pdb_file:
                self.logger.warning("No PDB file found for %s", target_name)
                return None
                
            # Get output directory
            data_dir = self._get_data_dir()
            output_dir = os.path.join(data_dir, "outputs", "docking", target_name)
            os.makedirs(output_dir, exist_ok=True)
            
            # Prepare structure for each binding site
            prepared_structures = []
            
            for i, site in enumerate(binding_sites):
                site_name = site.get("name", f"site_{i+1}")
                site_radius = site.get("radius", 10.0)
                site_center = site.get("center", {"x": 0, "y": 0, "z": 0})
                
                # Prepare structure for this site
                prepared_file = self.pdb_processor.prepare_for_docking(
                    pdb_file,
                    os.path.join(output_dir, f"{target_name}_{site_name}_prepared.pdb"),
                    center=[site_center["x"], site_center["y"], site_center["z"]],
                    radius=site_radius
                )
                
                if prepared_file:
                    prepared_structures.append({
                        "site_name": site_name,
                        "file_path": prepared_file,
                        "center": site_center,
                        "radius": site_radius
                    })
            
            # Create master prepared structure info
            prepared_info = {
                "target_name": target_name,
                "original_pdb": pdb_file,
                "prepared_structures": prepared_structures,
                "structure_source": structure_info.get("source", "unknown"),
                "binding_sites_count": len(binding_sites)
            }
            
            # Write master prepared structure file
            master_prepared_file = os.path.join(output_dir, f"{target_name}_prepared.json")
            with open(master_prepared_file, "w") as f:
                json.dump(prepared_info, f, indent=2)
                
            prepared_info["prepared_path"] = master_prepared_file
                
            self.logger.info("Prepared %d structures for %s", 
                             len(prepared_structures), target_name)
            
            return prepared_info
            
        except Exception as e:
            self.logger.error("Error preparing structure for %s: %s", 
                             target_name, e)
            
            self.learn_from_interaction({
                "interaction_type": "structure_preparation",
                "key": f"preparation:{target_name}",
                "target_name": target_name,
                "error_message": str(e),
                "success": False
            })
            
            return None
    
    def _get_data_dir(self) -> str:
        """
        Get the data directory path.
        
        Returns:
            Path to the data directory.
        """
        data_dir = self.config.get("data_dir")
        
        if not data_dir:
            # Try to determine the root directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            root_dir = os.path.dirname(os.path.dirname(current_dir))
            data_dir = os.path.join(root_dir, "data")
            
        return data_dir