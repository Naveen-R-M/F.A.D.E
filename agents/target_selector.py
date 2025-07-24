"""
Target Selector Agent for F.A.D.E framework.

This agent processes natural language queries to identify protein targets
and drug requirements, and generates configuration files for downstream processes.
"""

import os
import json
import logging
from typing import Any, Dict, List, Optional, Tuple

from agents.base.base_agent import BaseAgent
from utils.gemini_client import GeminiClient
from utils.uniprot_client import UniProtClient
from utils.config_generator import ConfigGenerator


class TargetSelector(BaseAgent):
    """
    Agent for processing natural language queries to identify protein targets
    and requirements for drug discovery.
    """
    
    def __init__(
        self, 
        name: str = "target_selector",
        config: Optional[Dict[str, Any]] = None,
        gemini_api_key: Optional[str] = None,
    ) -> None:
        """
        Initialize the Target Selector agent.
        
        Args:
            name: Unique identifier for the agent.
            config: Optional configuration parameters.
            gemini_api_key: API key for the Gemini model. If not provided,
                            it will be loaded from environment variables.
        """
        super().__init__(name, config)
        
        self.gemini_client = GeminiClient(api_key=gemini_api_key)
        self.uniprot_client = UniProtClient()
        self.config_generator = ConfigGenerator()
        
        # Set up logging
        self.logger = logging.getLogger(f"fade.agent.{name}")
        
    def process(self, input_data: str) -> Dict[str, Any]:
        """
        Process a natural language query to identify protein targets and requirements.
        
        Args:
            input_data: Natural language query describing the drug discovery goal.
            
        Returns:
            Dictionary containing parsed information, retrieved sequences,
            and generated configuration files.
        """
        self.logger.info("Processing query: %s", input_data)
        
        # Extract structured data from query
        parsed_data = self.parse_query(input_data)
        
        # Retrieve protein sequences
        sequences = {}
        for target in parsed_data.get("protein_targets", []):
            target_name = target.get("name")
            if not target_name:
                continue
                
            # Get sequence
            sequence_info = self.fetch_protein_sequence(target)
            if sequence_info:
                sequences[target_name] = sequence_info
        
        # Generate configurations
        config_files = self.generate_configs(parsed_data, sequences)
        
        return {
            "parsed_data": parsed_data,
            "sequences": sequences,
            "config_files": config_files
        }
    
    def parse_query(self, query: str) -> Dict[str, Any]:
        """
        Parse a natural language query to extract protein targets and requirements.
        
        Args:
            query: Natural language query describing the drug discovery goal.
            
        Returns:
            Structured data extracted from the query.
        """
        self.logger.info("Parsing query: %s", query)
        
        # Use Gemini to extract structured data
        parsed_data = self.gemini_client.extract_protein_info(query)
        
        self.logger.info("Extracted targets: %s", 
                         [t.get("name") for t in parsed_data.get("protein_targets", [])])
        
        return parsed_data
    
    def fetch_protein_sequence(self, target_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Fetch protein sequence information from UniProt.
        
        Args:
            target_info: Dictionary containing target protein information.
            
        Returns:
            Dictionary with sequence information or None if not found.
        """
        target_name = target_info.get("name")
        if not target_name:
            self.logger.warning("Target name not provided")
            return None
            
        self.logger.info("Fetching sequence for %s", target_name)
        
        # Try to get protein by gene name first
        organism = target_info.get("organism")
        protein = self.uniprot_client.get_protein_by_gene_name(target_name, organism)
        
        # If not found, try searching directly
        if not protein:
            search_results = self.uniprot_client.search_protein(target_name, limit=1)
            protein = search_results[0] if search_results else None
        
        if not protein:
            self.logger.warning("Protein not found: %s", target_name)
            return None
            
        # Get accession and sequence
        accession = protein.get("primaryAccession")
        sequence = protein.get("sequence", {}).get("value")
        
        if not sequence:
            self.logger.warning("Sequence not found for %s", target_name)
            return None
            
        # Apply mutations if specified
        mutations = target_info.get("mutations", [])
        mutated_sequence = sequence
        
        for mutation in mutations:
            original = mutation.get("original_residue")
            position = mutation.get("position")
            mutated = mutation.get("mutated_residue")
            
            if original and position and mutated:
                mutation_desc = f"{original}{position}{mutated}"
                try:
                    mutated_sequence = self.uniprot_client.apply_mutation(
                        mutated_sequence, mutation_desc
                    )
                    self.logger.info("Applied mutation %s", mutation_desc)
                except ValueError as e:
                    self.logger.error("Failed to apply mutation %s: %s", mutation_desc, e)
        
        # Create sequence info dictionary
        sequence_info = {
            "accession": accession,
            "original_sequence": sequence,
            "sequence": mutated_sequence,
            "length": len(mutated_sequence),
            "mutations": mutations,
            "organism": protein.get("organism", {}).get("scientificName"),
            "protein_name": protein.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", target_name),
            "gene_name": target_name
        }
        
        return sequence_info
    
    def generate_configs(
        self, 
        parsed_data: Dict[str, Any], 
        sequences: Dict[str, Dict[str, Any]]
    ) -> Dict[str, str]:
        """
        Generate configuration files for downstream processes.
        
        Args:
            parsed_data: Structured data extracted from the query.
            sequences: Dictionary mapping target names to sequence information.
            
        Returns:
            Dictionary mapping configuration types to file paths.
        """
        self.logger.info("Generating configuration files")
        
        config_files = {}
        data_dir = self._get_data_dir()
        
        # Create input directories
        sequences_dir = os.path.join(data_dir, "inputs", "sequences")
        os.makedirs(sequences_dir, exist_ok=True)
        
        configs_dir = os.path.join(data_dir, "inputs", "configs")
        os.makedirs(configs_dir, exist_ok=True)
        
        job_scripts_dir = os.path.join(data_dir, "inputs", "job_scripts")
        os.makedirs(job_scripts_dir, exist_ok=True)
        
        # Generate FASTA files and configurations for each target
        for target_name, sequence_info in sequences.items():
            # Save sequence to FASTA file
            fasta_path = os.path.join(sequences_dir, f"{target_name}.fasta")
            
            description = f"{sequence_info.get('protein_name')} [{sequence_info.get('organism')}]"
            if sequence_info.get('mutations'):
                mutation_str = ", ".join([
                    f"{m.get('original_residue')}{m.get('position')}{m.get('mutated_residue')}"
                    for m in sequence_info.get('mutations', [])
                ])
                description += f" Mutations: {mutation_str}"
                
            self.uniprot_client.save_fasta(
                fasta_path,
                sequence_info.get("accession", target_name),
                sequence_info.get("sequence"),
                description
            )
            
            config_files[f"{target_name}_fasta"] = fasta_path
            
            # Generate AlphaFold configuration and job script
            output_dir = os.path.join(data_dir, "outputs", "structures", target_name)
            
            alphafold_config_path = os.path.join(configs_dir, f"{target_name}_alphafold.json")
            alphafold_job_path = os.path.join(job_scripts_dir, f"{target_name}_alphafold.sh")
            
            self.config_generator.create_alphafold_job(
                protein_name=target_name,
                sequence_file=fasta_path,
                output_dir=output_dir,
                config_file=alphafold_config_path,
                job_script_path=alphafold_job_path
            )
            
            config_files[f"{target_name}_alphafold_config"] = alphafold_config_path
            config_files[f"{target_name}_alphafold_job"] = alphafold_job_path
        
        # Save parsed data
        parsed_data_path = os.path.join(configs_dir, "parsed_query.json")
        with open(parsed_data_path, "w") as f:
            json.dump(parsed_data, f, indent=2)
            
        config_files["parsed_query"] = parsed_data_path
        
        return config_files
    
    def _get_data_dir(self) -> str:
        """
        Get the data directory path.
        
        Returns:
            Path to the data directory.
        """
        data_dir = self.config.get("data_dir")
        
        if not data_dir:
            # Try to determine the framework root directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            root_dir = os.path.dirname(current_dir)
            data_dir = os.path.join(root_dir, "data")
            
        return data_dir
