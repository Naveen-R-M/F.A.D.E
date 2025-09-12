"""
Enhanced Target Selector Agent for F.A.D.E using RCSB PDB
This version uses RCSB PDB instead of UniProt for target selection
"""

import os
import json
import time
from typing import Any, Dict, List, Optional, Tuple

from agents.base.base_agent import BaseAgent
from agents.base.agentic_mixin import AgenticMixin
from agents.target_selector.error_analyzer import ErrorAnalyzer
from agents.target_selector.query_reformulator import QueryReformulator
from utils.gemini_client import GeminiClient
from utils.rcsb_client import RCSBClient
from utils.config_generator import ConfigGenerator


class RCSBTargetSelector(BaseAgent, AgenticMixin):
    """
    Enhanced Target Selector agent that uses RCSB PDB instead of UniProt.
    Processes natural language queries to identify protein targets and retrieve
    experimental structures directly from RCSB PDB.
    """
    
    def __init__(
        self, 
        name: str = "rcsb_target_selector",
        config: Optional[Dict[str, Any]] = None,
        gemini_api_key: Optional[str] = None,
        gemini_model: Optional[str] = None,
    ) -> None:
        """
        Initialize the RCSB Target Selector agent.
        
        Args:
            name: Unique identifier for the agent.
            config: Optional configuration parameters.
            gemini_api_key: API key for the Gemini model.
            gemini_model: Model name for Gemini.
        """
        # Initialize base agent
        BaseAgent.__init__(self, name, config)
        
        # Initialize clients
        self.gemini_client = GeminiClient(api_key=gemini_api_key, model=gemini_model)
        self.rcsb_client = RCSBClient()
        self.config_generator = ConfigGenerator()
        
        # Initialize agentic components
        AgenticMixin.initialize_agentic_components(self, llm_client=self.gemini_client)
        
        # Create specialized components
        self.error_analyzer = ErrorAnalyzer(llm_client=self.gemini_client)
        self.query_reformulator = QueryReformulator(llm_client=self.gemini_client)
        
        # Set up memory file
        data_dir = self._get_data_dir()
        memory_dir = os.path.join(data_dir, "memory")
        os.makedirs(memory_dir, exist_ok=True)
        self.memory_file = os.path.join(memory_dir, f"{name}_memory.json")
        
    def process(self, input_data: str) -> Dict[str, Any]:
        """
        Process a natural language query to identify protein targets and retrieve structures.
        
        Args:
            input_data: Natural language query describing the drug discovery goal.
            
        Returns:
            Dictionary containing parsed information, retrieved structures,
            and generated configuration files.
        """
        self.logger.info("Processing query with RCSB: %s", input_data)
        
        # Extract structured data from query
        parsed_data = self.parse_query(input_data)
        
        # Retrieve protein structures from RCSB
        structures = {}
        for target in parsed_data.get("protein_targets", []):
            target_name = target.get("name")
            if not target_name:
                continue
                
            # Get structure from RCSB
            structure_info = self.fetch_protein_structure(input_data, target, parsed_data)
            if structure_info:
                structures[target_name] = structure_info
        
        # Generate configurations
        config_files = self.generate_configs(parsed_data, structures)
        
        return {
            "parsed_data": parsed_data,
            "structures": structures,
            "config_files": config_files,
            "source": "rcsb_pdb"
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
        
        # Enhanced prompt for RCSB-focused parsing
        enhanced_prompt = f"""
        Parse this drug discovery query and extract structured information.
        Focus on identifying specific proteins that likely have experimental structures in PDB.
        
        Query: "{query}"
        
        Extract and return JSON with:
        {{
            "protein_targets": [
                {{
                    "name": "official protein name or gene symbol",
                    "aliases": ["alternative names"],
                    "organism": "organism name (prefer human)",
                    "mutations": [
                        {{
                            "original_residue": "single letter",
                            "position": number,
                            "mutated_residue": "single letter",
                            "notation": "standard notation like G12D"
                        }}
                    ]
                }}
            ],
            "disease_context": "disease or condition",
            "drug_requirements": {{
                "binding_affinity": "target affinity",
                "selectivity": "selectivity requirements",
                "pharmacokinetics": "ADMET requirements",
                "blood_brain_barrier": "BBB permeability needs"
            }},
            "binding_sites": ["specific sites if mentioned"],
            "molecule_preferences": {{
                "molecular_weight": "MW constraints",
                "lipophilicity": "LogP preferences",
                "rule_of_five": "Lipinski compliance"
            }}
        }}
        
        Prioritize well-studied proteins that likely have experimental structures.
        """
        
        # Use Gemini to extract structured data with enhanced prompt
        parsed_data = self.execute_with_retry(
            lambda: self.gemini_client.extract_protein_info(query),
            operation_name="Enhanced query parsing"
        )
        
        self.logger.info("Extracted targets: %s", 
                         [t.get("name") for t in parsed_data.get("protein_targets", [])])
        
        return parsed_data

    def _extract_with_enhanced_prompt(self, prompt: str) -> Dict[str, Any]:
        """
        Extract information using enhanced prompt.
        
        Args:
            prompt: Enhanced prompt for extraction
            
        Returns:
            Parsed data dictionary
        """
        response = self.gemini_client.generate_text(prompt)
        
        try:
            # Try to parse as JSON
            parsed_data = json.loads(response)
            return parsed_data
        except json.JSONDecodeError:
            # Fallback parsing if JSON parsing fails
            self.logger.warning("JSON parsing failed, using fallback extraction")
            return self.gemini_client.extract_protein_info(prompt.split('Query: "')[1].split('"')[0])
    
    def fetch_protein_structure(self, original_query: str, target_info: Dict[str, Any], 
                              full_parsed_data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Fetch protein structure from RCSB PDB with agentic error recovery.
        
        Args:
            original_query: Original user query
            target_info: Dictionary containing target protein information
            full_parsed_data: Complete parsed data for context
            
        Returns:
            Dictionary with structure information or None if not found
        """
        target_name = target_info.get("name")
        if not target_name:
            self.logger.warning("Target name not provided")
            return None
            
        self.logger.info("Fetching structure for %s from RCSB", target_name)
        
        # Set up output directory
        data_dir = self._get_data_dir()
        output_dir = os.path.join(data_dir, "outputs", "structures", target_name)
        os.makedirs(output_dir, exist_ok=True)
        
        # Track search attempts for learning
        search_attempts = []
        max_attempts = 3
        
        for attempt in range(max_attempts):
            try:
                # Use RCSB client to process the query
                structure_info = self.rcsb_client.process_target_query(
                    self.gemini_client,
                    original_query,
                    full_parsed_data,
                    output_dir
                )
                
                if structure_info:
                    # Enhance structure info with F.A.D.E-specific data
                    enhanced_info = self._enhance_structure_info(
                        structure_info, target_info, output_dir
                    )
                    
                    # Record successful attempt
                    search_attempts.append({
                        "attempt": attempt + 1,
                        "success": True,
                        "pdb_id": structure_info.get("pdb_id"),
                        "method": structure_info.get("source"),
                        "timestamp": time.time()
                    })
                    
                    # Learn from successful interaction
                    self.learn_from_interaction({
                        "interaction_type": "structure_retrieval",
                        "key": f"rcsb:{target_name}",
                        "target_info": target_info,
                        "attempts": len(search_attempts),
                        "success": True,
                        "structure_info": {
                            "pdb_id": structure_info.get("pdb_id"),
                            "source": structure_info.get("source"),
                            "resolution": structure_info.get("resolution"),
                            "organism": structure_info.get("organism")
                        }
                    })
                    
                    return enhanced_info
                else:
                    # Record failed attempt
                    search_attempts.append({
                        "attempt": attempt + 1,
                        "success": False,
                        "error": "No structure found",
                        "timestamp": time.time()
                    })
                    
            except Exception as e:
                error_message = str(e)
                self.logger.error("Error during RCSB search: %s", error_message)
                
                # Record error
                search_attempts.append({
                    "attempt": attempt + 1,
                    "success": False,
                    "error": error_message,
                    "timestamp": time.time()
                })
                
                # Analyze error for potential recovery
                error_analysis = self.error_analyzer.analyze(
                    error_message,
                    {"operation": "rcsb_structure_retrieval", "target": target_name}
                )
                
                # Learn from error
                self.learn_from_interaction({
                    "interaction_type": "error_recovery",
                    "key": f"rcsb_error:{target_name}",
                    "error_message": error_message,
                    "error_analysis": error_analysis,
                    "success": False
                })
        
        # If all attempts fail
        self.logger.warning("Failed to find structure for %s after %d attempts", 
                           target_name, max_attempts)
        
        # Learn from complete failure
        self.learn_from_interaction({
            "interaction_type": "structure_retrieval",
            "key": f"rcsb:{target_name}",
            "target_info": target_info,
            "attempts": len(search_attempts),
            "success": False,
            "reason": "Max attempts reached"
        })
        
        return None
    
    def _enhance_structure_info(self, structure_info: Dict[str, Any], 
                              target_info: Dict[str, Any], output_dir: str) -> Dict[str, Any]:
        """
        Enhance structure information with F.A.D.E-specific data.
        
        Args:
            structure_info: Basic structure info from RCSB client
            target_info: Target information
            output_dir: Output directory
            
        Returns:
            Enhanced structure information
        """
        enhanced = dict(structure_info)
        
        # Add F.A.D.E-specific fields
        enhanced.update({
            "target_name": target_info.get("name"),
            "gene_name": target_info.get("name"), 
            "mutations": target_info.get("mutations", []),
            "organism": structure_info.get("organism", "unknown"),
            "structure_type": "experimental_pdb",
            "confidence_scores": {
                "overall": 1.0,  # High confidence for experimental structures
                "method": "experimental",
                "source": "rcsb_pdb"
            },
            "binding_sites": [],  # Will be filled by structure analysis
            "druggable_sites": []  # Will be filled by binding site detection
        })
        
        # Create FASTA file from PDB if needed
        pdb_file = structure_info.get("pdb_file")
        if pdb_file and os.path.exists(pdb_file):
            fasta_file = self._extract_sequence_from_pdb(pdb_file, output_dir, target_info.get("name"))
            if fasta_file:
                enhanced["fasta_file"] = fasta_file
        
        return enhanced
    
    def _extract_sequence_from_pdb(self, pdb_file: str, output_dir: str, target_name: str) -> Optional[str]:
        """
        Extract protein sequence from PDB file and save as FASTA.
        
        Args:
            pdb_file: Path to PDB file
            output_dir: Output directory
            target_name: Target name for FASTA header
            
        Returns:
            Path to FASTA file or None
        """
        try:
            # Try using BioPython if available
            try:
                from Bio import PDB
                from Bio import SeqIO
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord
                
                parser = PDB.PDBParser(QUIET=True)
                structure = parser.get_structure("structure", pdb_file)
                
                sequences = []
                for model in structure:
                    for chain in model:
                        seq = ""
                        for residue in chain:
                            if PDB.is_aa(residue):
                                try:
                                    seq += PDB.protein_letters_3to1[residue.get_resname()]
                                except KeyError:
                                    seq += "X"  # Unknown residue
                        
                        if seq:
                            sequences.append((chain.id, seq))
                
                if sequences:
                    # Use the longest sequence (usually chain A)
                    longest_chain, longest_seq = max(sequences, key=lambda x: len(x[1]))
                    
                    fasta_path = os.path.join(output_dir, f"{target_name}.fasta")
                    
                    record = SeqRecord(
                        Seq(longest_seq),
                        id=f"{target_name}_{longest_chain}",
                        description=f"{target_name} from PDB structure (chain {longest_chain})"
                    )
                    
                    with open(fasta_path, "w") as f:
                        SeqIO.write(record, f, "fasta")
                    
                    self.logger.info(f"Extracted sequence from PDB: {fasta_path}")
                    return fasta_path
                    
            except ImportError:
                self.logger.warning("BioPython not available for PDB sequence extraction")
                pass
            
            # Fallback: simple PDB parsing
            return self._simple_pdb_sequence_extraction(pdb_file, output_dir, target_name)
            
        except Exception as e:
            self.logger.error(f"Failed to extract sequence from PDB: {e}")
            return None
    
    def _simple_pdb_sequence_extraction(self, pdb_file: str, output_dir: str, target_name: str) -> Optional[str]:
        """
        Simple PDB sequence extraction without BioPython.
        
        Args:
            pdb_file: Path to PDB file
            output_dir: Output directory
            target_name: Target name
            
        Returns:
            Path to FASTA file or None
        """
        # Simple 3-letter to 1-letter amino acid code mapping
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        
        try:
            sequences = {}
            
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') and line[12:16].strip() == 'CA':  # Alpha carbon only
                        chain_id = line[21]
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26].strip())
                        
                        if chain_id not in sequences:
                            sequences[chain_id] = {}
                        
                        if res_name in aa_map:
                            sequences[chain_id][res_num] = aa_map[res_name]
            
            if sequences:
                # Find the longest chain
                longest_chain = max(sequences.keys(), key=lambda x: len(sequences[x]))
                chain_seq = sequences[longest_chain]
                
                # Build sequence in order
                sorted_positions = sorted(chain_seq.keys())
                sequence = ''.join([chain_seq[pos] for pos in sorted_positions])
                
                if sequence:
                    fasta_path = os.path.join(output_dir, f"{target_name}.fasta")
                    
                    with open(fasta_path, 'w') as f:
                        f.write(f">{target_name}_{longest_chain}|PDB|{target_name} from PDB structure\n")
                        # Write sequence in 80-character lines
                        for i in range(0, len(sequence), 80):
                            f.write(sequence[i:i+80] + '\n')
                    
                    self.logger.info(f"Extracted sequence (simple method): {fasta_path}")
                    return fasta_path
            
            return None
            
        except Exception as e:
            self.logger.error(f"Simple PDB parsing failed: {e}")
            return None
    
    def generate_configs(
        self, 
        parsed_data: Dict[str, Any], 
        structures: Dict[str, Dict[str, Any]]
    ) -> Dict[str, str]:
        """
        Generate configuration files for downstream processes using PDB structures.
        
        Args:
            parsed_data: Structured data extracted from the query.
            structures: Dictionary mapping target names to structure information.
            
        Returns:
            Dictionary mapping configuration types to file paths.
        """
        self.logger.info("Generating configuration files for RCSB structures")
        
        config_files = {}
        data_dir = self._get_data_dir()
        
        # Create input directories
        configs_dir = os.path.join(data_dir, "inputs", "configs")
        os.makedirs(configs_dir, exist_ok=True)
        
        # Generate configurations for each target with structure
        for target_name, structure_info in structures.items():
            try:
                pdb_file = structure_info.get("pdb_file")
                fasta_file = structure_info.get("fasta_file")
                
                if not pdb_file or not os.path.exists(pdb_file):
                    self.logger.warning(f"PDB file not found for {target_name}")
                    continue
                
                # Create target-specific configuration
                target_config = {
                    "target_name": target_name,
                    "pdb_id": structure_info.get("pdb_id"),
                    "pdb_file": pdb_file,
                    "fasta_file": fasta_file,
                    "source": "rcsb_pdb",
                    "structure_info": {
                        "resolution": structure_info.get("resolution"),
                        "method": structure_info.get("method"),
                        "organism": structure_info.get("organism"),
                        "ligands": structure_info.get("ligands", [])
                    },
                    "requirements": parsed_data.get("drug_requirements", {}),
                    "binding_sites": structure_info.get("binding_sites", []),
                    "mutations": structure_info.get("mutations", [])
                }
                
                # Save target configuration
                config_path = os.path.join(configs_dir, f"{target_name}_rcsb_config.json")
                with open(config_path, "w") as f:
                    json.dump(target_config, f, indent=2)
                
                config_files[f"{target_name}_config"] = config_path
                
                # Link PDB file to standard location for downstream processes
                config_files[f"{target_name}_pdb"] = pdb_file
                if fasta_file:
                    config_files[f"{target_name}_fasta"] = fasta_file
                
            except Exception as e:
                self.logger.error(f"Failed to generate config for {target_name}: {e}")
        
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
            root_dir = os.path.dirname(os.path.dirname(current_dir))
            data_dir = os.path.join(root_dir, "data")
            
        return data_dir
