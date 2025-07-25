"""
Target Selector Agent for F.A.D.E

This agent processes natural language queries to identify protein targets
and drug requirements, and generates configuration files for downstream processes.
"""

import os
import json
import time
from typing import Any, Dict, List, Optional, Tuple

from agents.base.base_agent import BaseAgent
from agents.base.agentic_mixin import AgenticMixin
from agents.target_selector.error_analyzer import ErrorAnalyzer
from agents.target_selector.query_reformulator import QueryReformulator
from agents.target_selector.sequence_validator import SequenceValidator
from agents.target_selector.strategy_selector import SearchStrategySelector
from utils.gemini_client import GeminiClient
from utils.uniprot_client import UniProtClient
from utils.config_generator import ConfigGenerator


class TargetSelector(BaseAgent, AgenticMixin):
    """
    Agent for processing natural language queries to identify protein targets
    and requirements for drug discovery.
    """
    
    def __init__(
        self, 
        name: str = "target_selector",
        config: Optional[Dict[str, Any]] = None,
        gemini_api_key: Optional[str] = None,
        gemini_model: Optional[str] = None,
    ) -> None:
        """
        Initialize the Target Selector agent.
        
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
        self.uniprot_client = UniProtClient()
        self.config_generator = ConfigGenerator()
        
        # Initialize agentic components
        AgenticMixin.initialize_agentic_components(self, llm_client=self.gemini_client)
        
        # Create specialized components
        self.error_analyzer = ErrorAnalyzer(llm_client=self.gemini_client)
        self.query_reformulator = QueryReformulator(llm_client=self.gemini_client)
        self.sequence_validator = SequenceValidator(llm_client=self.gemini_client)
        self.strategy_selector = SearchStrategySelector(llm_client=self.gemini_client)
        
        # Set up memory file
        data_dir = self._get_data_dir()
        memory_dir = os.path.join(data_dir, "memory")
        os.makedirs(memory_dir, exist_ok=True)
        self.memory_file = os.path.join(memory_dir, f"{name}_memory.json")
        
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
        parsed_data = self.execute_with_retry(
            self.gemini_client.extract_protein_info,
            query,
            operation_name="Query parsing"
        )
        
        self.logger.info("Extracted targets: %s", 
                         [t.get("name") for t in parsed_data.get("protein_targets", [])])
        
        return parsed_data
    
    def fetch_protein_sequence(self, target_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Fetch protein sequence information from UniProt with agentic error recovery.
        
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
        
        # Track search attempts
        search_attempts = []
        max_attempts = 5
        
        # Initial search strategy
        current_strategy = self.strategy_selector.select_strategy(
            target_info, 
            search_attempts
        )
        
        # Attempt retrieval with automatic recovery
        for attempt in range(max_attempts):
            try:
                # Apply current strategy
                if current_strategy["method"] == "gene_name":
                    protein = self.uniprot_client.get_protein_by_gene_name(
                        current_strategy.get("gene_name", target_name), 
                        current_strategy.get("organism")
                    )
                elif current_strategy["method"] == "search":
                    search_query = current_strategy.get("query", target_name)
                    search_results = self.uniprot_client.search_protein(
                        search_query, 
                        limit=current_strategy.get("limit", 1)
                    )
                    protein = search_results[0] if search_results else None
                elif current_strategy["method"] == "accession":
                    protein = self.uniprot_client.get_protein_by_accession(
                        current_strategy.get("accession")
                    )
                else:
                    protein = None
                    
                # Record attempt
                search_attempts.append({
                    "strategy": current_strategy,
                    "success": protein is not None,
                    "error": None,
                    "timestamp": time.time()
                })
                
                # If found, validate sequence
                if protein:
                    accession = protein.get("primaryAccession")
                    sequence = protein.get("sequence", {}).get("value")
                    
                    if not sequence:
                        self.logger.warning("Sequence not found for %s", target_name)
                        # Record validation failure
                        search_attempts[-1]["validation_error"] = "No sequence found"
                        
                        # Get a new strategy
                        current_strategy = self.strategy_selector.select_strategy(
                            target_info, 
                            search_attempts
                        )
                        continue
                    
                    # Prepare metadata for validation
                    metadata = {
                        "accession": accession,
                        "protein_name": protein.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", target_name),
                        "organism": protein.get("organism", {}).get("scientificName"),
                        "mutations": target_info.get("mutations", [])
                    }
                    
                    # Validate sequence
                    is_valid, reason, validation_results = self.sequence_validator.validate(
                        target_name, 
                        sequence,
                        metadata
                    )
                    
                    # Record validation results
                    search_attempts[-1]["validation"] = {
                        "is_valid": is_valid,
                        "reason": reason,
                        "score": validation_results.get("overall_score", 0.0)
                    }
                    
                    if is_valid:
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
                                        mutated_sequence, 
                                        mutation_desc
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
                            "gene_name": target_name,
                            "validation_results": validation_results
                        }
                        
                        # Update search attempts with final outcome
                        search_attempts[-1]["final_outcome"] = "success"
                        
                        # Learn from successful interaction
                        self.learn_from_interaction({
                            "interaction_type": "protein_retrieval",
                            "key": f"retrieval:{target_name}",
                            "target_info": target_info,
                            "successful_strategy": current_strategy,
                            "attempts": len(search_attempts),
                            "success": True,
                            "sequence_info": {
                                "accession": accession,
                                "length": len(sequence),
                                "organism": metadata.get("organism")
                            }
                        })
                        
                        return sequence_info
                    else:
                        self.logger.warning("Invalid sequence: %s", reason)
                        
                        # Record validation failure
                        search_attempts[-1]["validation_error"] = reason
                
                # If we reach here, we need a new strategy
                current_strategy = self.strategy_selector.select_strategy(
                    target_info, 
                    search_attempts
                )
                
            except Exception as e:
                # Handle error
                error_message = str(e)
                self.logger.error("Error during search: %s", error_message)
                
                # Record error
                search_attempts.append({
                    "strategy": current_strategy,
                    "success": False,
                    "error": error_message,
                    "timestamp": time.time()
                })
                
                # Analyze error
                error_analysis = self.error_analyzer.analyze(
                    error_message, 
                    {"operation": "protein_retrieval", "target": target_name, "strategy": current_strategy}
                )
                
                # Generate new strategies
                alternative_strategies = self.query_reformulator.reformulate(
                    current_strategy, 
                    error_analysis,
                    target_info
                )
                
                if alternative_strategies:
                    current_strategy = alternative_strategies[0]
                else:
                    # Fall back to simple search if no alternatives
                    current_strategy = {
                        "method": "search",
                        "query": target_name,
                        "limit": 5,
                        "description": "Fallback search"
                    }
                
                # Learn from failed interaction
                self.learn_from_interaction({
                    "interaction_type": "error_recovery",
                    "key": f"error:{target_name}",
                    "error_message": error_message,
                    "error_analysis": error_analysis,
                    "next_strategy": current_strategy,
                    "success": False
                })
        
        # If all attempts fail, return None and record failure
        self.logger.warning("Failed to find valid sequence for %s after %d attempts", 
                           target_name, max_attempts)
        
        # Learn from the failed retrieval
        self.learn_from_interaction({
            "interaction_type": "protein_retrieval",
            "key": f"retrieval:{target_name}",
            "target_info": target_info,
            "attempts": len(search_attempts),
            "success": False,
            "reason": "Max attempts reached"
        })
        
        return None
    
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
                
            try:
                # Use error-handled version for critical operations
                self.execute_with_retry(
                    self.uniprot_client.save_fasta,
                    fasta_path,
                    sequence_info.get("accession", target_name),
                    sequence_info.get("sequence"),
                    description,
                    operation_name=f"Save FASTA for {target_name}"
                )
                
                config_files[f"{target_name}_fasta"] = fasta_path
                
                # Generate AlphaFold configuration and job script
                output_dir = os.path.join(data_dir, "outputs", "structures", target_name)
                
                alphafold_config_path = os.path.join(configs_dir, f"{target_name}_alphafold.json")
                alphafold_job_path = os.path.join(job_scripts_dir, f"{target_name}_alphafold.sh")
                
                self.execute_with_retry(
                    self.config_generator.create_alphafold_job,
                    protein_name=target_name,
                    sequence_file=fasta_path,
                    output_dir=output_dir,
                    config_file=alphafold_config_path,
                    job_script_path=alphafold_job_path,
                    operation_name=f"Generate AlphaFold config for {target_name}"
                )
                
                config_files[f"{target_name}_alphafold_config"] = alphafold_config_path
                config_files[f"{target_name}_alphafold_job"] = alphafold_job_path
                
            except Exception as e:
                self.logger.error(f"Failed to generate config for {target_name}: {e}")
                # Use error analysis to understand the issue
                error_analysis = self.error_analyzer.analyze(
                    str(e), 
                    {"operation": "config_generation", "target": target_name}
                )
                
                # Record the error
                self.learn_from_interaction({
                    "interaction_type": "config_generation",
                    "key": f"config:{target_name}",
                    "error_message": str(e),
                    "error_analysis": error_analysis,
                    "success": False
                })
        
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
            # Try to determine the root directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            root_dir = os.path.dirname(current_dir)
            data_dir = os.path.join(root_dir, "data")
            
        return data_dir
