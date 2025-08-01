"""
Structure Predictor Agent for F.A.D.E

This agent handles the generation and validation of 3D protein structures
from sequence data, and identifies potential binding sites for drug targeting.
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


class StructurePredictor(BaseAgent, AgenticMixin):
    """
    Agent for predicting, validating, and analyzing protein structures.
    
    This agent:
    1. Takes protein sequence information from Target Selector
    2. Runs structure prediction using AlphaFold3
    3. Validates the predicted structures
    4. Identifies potential binding sites
    5. Prepares the structures for docking
    """
    
    def __init__(
        self,
        name: str = "structure_predictor",
        config: Optional[Dict[str, Any]] = None,
        gemini_api_key: Optional[str] = None,
        gemini_model: Optional[str] = None,
    ) -> None:
        """
        Initialize the Structure Predictor agent.
        
        Args:
            name: Unique identifier for the agent.
            config: Optional configuration parameters.
            gemini_api_key: API key for the Gemini model. If not provided,
                            it will be loaded from environment variables.
            gemini_model: Model name for Gemini. If not provided,
                          it will be loaded from environment variables.
        """
        # Initialize base agent
        BaseAgent.__init__(self, name, config)
        
        # Initialize clients
        self.gemini_client = GeminiClient(api_key=gemini_api_key, model=gemini_model)
        self.slurm_client = SlurmClient()
        self.alphafold_client = AlphaFoldClient()
        self.config_generator = ConfigGenerator()
        
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
        
        # Results containers
        structures = {}
        binding_sites = {}
        prepared_structures = {}
        
        # Process each target
        for target_name, sequence_info in sequences.items():
            # Get AlphaFold job config
            job_config_key = f"{target_name}_alphafold_job"
            job_script_path = job_configs.get(job_config_key)
            
            if not job_script_path:
                self.logger.warning("No job script found for %s, skipping", target_name)
                continue
                
            # Run structure prediction job
            self.logger.info("Running structure prediction for %s", target_name)
            
            structure_info = self.predict_structure(target_name, sequence_info, job_script_path)
            
            if not structure_info:
                self.logger.warning("Failed to predict structure for %s", target_name)
                continue
                
            structures[target_name] = structure_info
            
            # Identify binding sites
            self.logger.info("Identifying binding sites for %s", target_name)
            
            target_binding_sites = self.identify_binding_sites(
                target_name, 
                structure_info, 
                sequence_info
            )
            
            binding_sites[target_name] = target_binding_sites
            
            # Prepare structure for docking
            self.logger.info("Preparing structure for docking: %s", target_name)
            
            prepared_structure_path = self.prepare_structure_for_docking(
                target_name,
                structure_info,
                target_binding_sites
            )
            
            if prepared_structure_path:
                prepared_structures[target_name] = prepared_structure_path
        
        return {
            "structures": structures,
            "binding_sites": binding_sites,
            "prepared_structures": prepared_structures
        }
    
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
        # Submit job to SLURM
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
                
            # Get output path from config
            data_dir = self._get_data_dir()
            output_dir = os.path.join(data_dir, "outputs", "structures", target_name)
            
            # Get PDB file
            pdb_file = self.alphafold_client.get_best_model_path(output_dir)
            
            if not pdb_file or not os.path.exists(pdb_file):
                self.logger.error("No PDB file found for %s", target_name)
                return None
                
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
                "confidence_scores": {
                    "overall": validation_results.get("overall_score", 0.0),
                    "plddt": validation_results.get("plddt", 0.0),
                    "ptm": validation_results.get("ptm", 0.0)
                },
                "validation_results": validation_results,
                "chain_info": structure_data.get("chains", {}),
                "residue_count": structure_data.get("residue_count", 0),
                "atom_count": structure_data.get("atom_count", 0),
                "secondary_structure": structure_data.get("secondary_structure", {})
            }
            
            # Log structure summary
            self.logger.info("Structure prediction for %s: %d residues, %d atoms, confidence: %.2f", 
                             target_name, 
                             structure_info["residue_count"],
                             structure_info["atom_count"],
                             structure_info["confidence_scores"]["overall"])
            
            return structure_info
            
        except Exception as e:
            self.logger.error("Error predicting structure for %s: %s", target_name, e)
            
            # Learn from failure
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
            
            # Identify binding sites using our detector
            detected_sites = self.binding_site_detector.detect(
                pdb_file,
                known_binding_sites=known_binding_sites
            )
            
            # Process mutations if present
            mutations = sequence_info.get("mutations", [])
            if mutations:
                # Add mutation sites as potential binding sites
                for mutation in mutations:
                    position = mutation.get("position")
                    mutated = mutation.get("mutated_residue")
                    original = mutation.get("original_residue")
                    
                    if position and mutated and original:
                        # Check if this mutation is within an existing binding site
                        in_binding_site = False
                        for site in detected_sites:
                            if position in site.get("residue_ids", []):
                                in_binding_site = True
                                break
                        
                        # If not in an existing site, add it as a new site
                        if not in_binding_site:
                            mutation_site = {
                                "name": f"Mutation_{original}{position}{mutated}",
                                "type": "mutation",
                                "center": {"x": 0, "y": 0, "z": 0},  # Will be calculated by the detector
                                "radius": 10.0,  # Default radius around mutation
                                "residue_ids": [position],
                                "score": 0.9,  # High score for known mutations
                                "description": f"Mutation site {original}{position}{mutated}"
                            }
                            detected_sites.append(mutation_site)
            
            # Log detected sites
            self.logger.info("Detected %d binding sites for %s", 
                             len(detected_sites), target_name)
            
            return detected_sites
            
        except Exception as e:
            self.logger.error("Error identifying binding sites for %s: %s", 
                             target_name, e)
            
            # Learn from failure
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
    ) -> Optional[str]:
        """
        Prepare the protein structure for docking.
        
        Args:
            target_name: Name of the target protein.
            structure_info: Dictionary containing structure information.
            binding_sites: List of dictionaries describing binding sites.
            
        Returns:
            Path to the prepared structure file or None if preparation failed.
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
            
            # Create master prepared structure file with all sites
            master_prepared_file = os.path.join(output_dir, f"{target_name}_prepared.json")
            
            with open(master_prepared_file, "w") as f:
                json.dump({
                    "target_name": target_name,
                    "original_pdb": pdb_file,
                    "prepared_structures": prepared_structures
                }, f, indent=2)
                
            # Log results
            self.logger.info("Prepared %d structures for %s", 
                             len(prepared_structures), target_name)
            
            return master_prepared_file
            
        except Exception as e:
            self.logger.error("Error preparing structure for %s: %s", 
                             target_name, e)
            
            # Learn from failure
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
