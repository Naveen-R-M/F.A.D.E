"""
Molecule Generator Agent for F.A.D.E

This agent generates potential drug candidates based on protein target information
and identified binding sites.
"""

import os
import json
from typing import Any, Dict, List, Optional

from agents.base.base_agent import BaseAgent
from agents.base.agentic_mixin import AgenticMixin
from utils.gemini_client import GeminiClient
from utils.logging import get_logger

from agents.molecule_generator.molecule_builder import MoleculeBuilder
from agents.molecule_generator.property_calculator import PropertyCalculator
from agents.molecule_generator.molecule_filter import MoleculeFilter


class MoleculeGenerator(BaseAgent, AgenticMixin):
    """
    Agent for generating, calculating properties for, and filtering potential drug molecules.
    """

    def __init__(
        self,
        name: str = "molecule_generator",
        config: Optional[Dict[str, Any]] = None,
        gemini_api_key: Optional[str] = None,
        gemini_model: Optional[str] = None,
    ) -> None:
        """
        Initialize the Molecule Generator agent.

        Args:
            name: Unique identifier for the agent.
            config: Optional configuration parameters.
            gemini_api_key: API key for the Gemini model. If not provided,
                            it will be loaded from environment variables.
            gemini_model: Model name for Gemini. If not provided,
                          it will be loaded from environment variables.
        """
        BaseAgent.__init__(self, name, config)

        self.gemini_client = GeminiClient(api_key=gemini_api_key, model=gemini_model)
        AgenticMixin.initialize_agentic_components(self, llm_client=self.gemini_client)

        self.molecule_builder = MoleculeBuilder(llm_client=self.gemini_client)
        self.property_calculator = PropertyCalculator()
        self.molecule_filter = MoleculeFilter()

        # Set up memory file
        data_dir = self._get_data_dir()
        memory_dir = os.path.join(data_dir, "memory")
        os.makedirs(memory_dir, exist_ok=True)
        self.memory_file = os.path.join(memory_dir, f"{name}_memory.json")
        self.logger.info(f"Initialized {self.name} agent with memory file: {self.memory_file}")

    def process(self, input_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate molecules based on protein structure and binding site information.

        Args:
            input_data: Dictionary containing:
                - prepared_structures: Dictionary mapping target names to prepared structure paths
                - binding_sites: Dictionary mapping target names to binding site info

        Returns:
            Dictionary containing generated molecules and their properties.
        """
        self.logger.info("Starting molecule generation process.")

        prepared_structures = input_data.get("prepared_structures", {})
        binding_sites_info = input_data.get("binding_sites", {})
        parsed_query_data = input_data.get("parsed_query_data", {})
        
        # Extract overall molecular property requirements from the parsed query
        molecular_requirements = parsed_query_data.get("molecule_properties", {})

        generated_molecules = {}

        for target_name, structure_path_info in prepared_structures.items():
            self.logger.info(f"Generating molecules for target: {target_name}")
            
            # Get binding sites specific to this target
            target_binding_sites = binding_sites_info.get(target_name, [])
            
            if not target_binding_sites:
                self.logger.warning(f"No binding sites found for {target_name}, skipping molecule generation.")
                continue

            target_candidates = []

            for i, site in enumerate(target_binding_sites):
                self.logger.info(f"Generating molecules for {target_name} - Binding Site {site.get('name', f'#{i+1}')}")
                
                # 1. LLM-guided scaffold/fragment generation
                smiles_suggestions = self.execute_with_retry(
                    self.molecule_builder.generate_smiles_from_binding_site,
                    site,
                    operation_name=f"LLM SMILES generation for {target_name} site {i+1}"
                )

                if not smiles_suggestions:
                    self.logger.warning(f"No SMILES suggestions from LLM for {target_name} site {i+1}.")
                    continue

                for smiles in smiles_suggestions:
                    # 2. RDKit-based molecule construction and 3D generation
                    mol = self.execute_with_retry(
                        self.molecule_builder.build_molecule_from_smiles,
                        smiles,
                        operation_name=f"Build RDKit Mol from SMILES {smiles}"
                    )
                    
                    if mol is None:
                        continue # Skip invalid SMILES

                    mol_3d = self.execute_with_retry(
                        self.molecule_builder.generate_3d_conformation,
                        mol,
                        operation_name=f"Generate 3D conformation for {smiles}"
                    )

                    if mol_3d is None:
                        continue # Skip if 3D generation fails

                    # 3. Property calculation
                    properties = self.property_calculator.calculate_properties(mol_3d)

                    # Store original SMILES and RDKit Mol for filtering
                    candidate_data = {
                        "smiles": smiles,
                        "mol": mol_3d, # Keep RDKit Mol object for filtering
                        "properties": properties,
                        "binding_site_name": site.get('name', f'Site {i+1}')
                    }
                    target_candidates.append(candidate_data)
            
            # 4. Filtering using MoleculeFilter
            # Pass the overall molecular requirements from the parsed query
            filtered_candidates = self.molecule_filter.filter_molecules(
                target_candidates, 
                molecular_requirements
            )

            # Convert filtered RDKit Mol objects to SDF strings for output
            final_candidates_for_output = []
            for candidate in filtered_candidates:
                sdf_string = self.molecule_builder.mol_to_sdf_string(
                    candidate["mol"], 
                    mol_name=f"{target_name}_{candidate['binding_site_name']}_{candidate['smiles']}"
                )
                if sdf_string:
                    final_candidates_for_output.append({
                        "smiles": candidate["smiles"],
                        "properties": candidate["properties"],
                        "sdf_string": sdf_string,
                        "binding_site_name": candidate["binding_site_name"]
                    })
            
            generated_molecules[target_name] = {
                "status": "molecules_generated",
                "candidates": final_candidates_for_output,
                "num_generated": len(target_candidates),
                "num_filtered": len(final_candidates_for_output)
            }
            self.logger.info(f"Generated {len(target_candidates)} and filtered to {len(final_candidates_for_output)} molecules for {target_name}.")

        return {"generated_molecules": generated_molecules}

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