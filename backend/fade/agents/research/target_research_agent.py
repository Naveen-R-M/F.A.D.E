"""
Target Research Agent for identifying and researching protein targets.

This agent processes natural language queries to identify protein targets,
retrieve their sequences, and find known compounds.
"""

import re
from typing import Dict, Any, List, Optional
from datetime import datetime

from langchain.chat_models import init_chat_model
from langchain_core.messages import HumanMessage, SystemMessage

from fade.agents.base_agent import BaseAgent
from fade.state import DrugDiscoveryState, ProteinTarget, KnownCompound
from fade.tools import get_uniprot_client, get_rcsb_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("agent.target_research")


class TargetResearchAgent(BaseAgent):
    """
    Agent responsible for researching protein targets and known compounds.
    
    This agent:
    1. Parses natural language queries to extract target information
    2. Searches UniProt for protein details
    3. Retrieves known compounds from ChEMBL
    4. Checks for existing PDB structures
    """
    
    def __init__(self):
        """Initialize the Target Research Agent."""
        super().__init__(
            name="target_research",
            description="Research protein targets and find known compounds",
            checkpoint_dir=config.CHECKPOINT_DIR / "target_research",
            max_retries=3
        )
        
        # Initialize LLM for query parsing
        llm_config = config.get_llm_config()
        self.llm = init_chat_model(**llm_config)
        
        # Initialize API clients
        self.uniprot_client = get_uniprot_client()
        self.chembl_client = get_chembl_client()
        self.rcsb_client = get_rcsb_client()
        
    def execute(self, state: DrugDiscoveryState) -> DrugDiscoveryState:
        """
        Execute target research based on the query.
        
        Args:
            state: Current pipeline state
            
        Returns:
            Updated state with target information and known compounds
        """
        query = state["query"]
        logger.info(f"Researching target from query: {query}")
        
        # Step 1: Parse query to extract target information
        target_info = self._parse_query(query)
        
        # Step 2: Search UniProt for protein information
        if target_info.get("gene_name") or target_info.get("uniprot_id"):
            protein_data = self._search_uniprot(target_info)
            
            if protein_data:
                # Update target info with UniProt data
                target_info.update(protein_data)
                state["target_info"] = target_info
                
                # Step 3: Search for known compounds
                if target_info.get("uniprot_id"):
                    known_compounds = self._search_known_compounds(target_info["uniprot_id"])
                    state["known_compounds"] = known_compounds
                    
                    # Step 4: Check for existing structures
                    existing_structures = self.rcsb_client.search_by_uniprot(
                        target_info["uniprot_id"],
                        limit=5
                    )
                    
                    if existing_structures:
                        logger.info(f"Found {len(existing_structures)} existing PDB structures")
                        # Store structure info for later use
                        if not state.get("structure_validation"):
                            state["structure_validation"] = {}
                        state["structure_validation"]["existing_structures"] = existing_structures
            else:
                error_msg = f"Could not find protein information for target: {target_info}"
                state = self._add_error(state, error_msg)
                state["should_continue"] = False
        else:
            error_msg = "Could not extract target information from query"
            state = self._add_error(state, error_msg)
            state["should_continue"] = False
        
        # Update current step
        state["current_step"] = "target_research_complete"
        
        # Log summary
        if state.get("target_info"):
            logger.info(f"Target identified: {state['target_info'].get('protein_name')} ({state['target_info'].get('uniprot_id')})")
            if state.get("known_compounds"):
                logger.info(f"Found {len(state['known_compounds'])} known compounds")
        
        return state
    
    def validate_input(self, state: DrugDiscoveryState) -> Dict[str, Any]:
        """
        Validate that the state has required input.
        
        Args:
            state: Current pipeline state
            
        Returns:
            Validation result
        """
        errors = []
        warnings = []
        
        # Check for query
        if not state.get("query"):
            errors.append("No query provided")
        
        # Check query length
        if state.get("query") and len(state["query"]) < 10:
            warnings.append("Query seems very short - might not contain enough information")
        
        return {
            "is_valid": len(errors) == 0,
            "errors": errors,
            "warnings": warnings
        }
    
    def _parse_query(self, query: str) -> ProteinTarget:
        """
        Parse natural language query to extract target information.
        
        Args:
            query: Natural language query
            
        Returns:
            Extracted target information
        """
        # Use LLM to parse the query
        system_prompt = """You are a bioinformatics expert. Extract protein target information from the user's drug discovery query.
        
        Return a JSON object with these fields:
        - gene_name: The gene symbol (e.g., KRAS, EGFR)
        - protein_name: The protein name if mentioned
        - uniprot_id: UniProt ID if mentioned (e.g., P01112)
        - organism: The organism (default to "human" if not specified)
        - mutations: Array of mutations if mentioned (e.g., ["G12C", "G12D"])
        - disease_context: The disease context if mentioned
        
        Be precise and only include information explicitly mentioned or clearly implied."""
        
        user_prompt = f"Extract target information from this query: {query}"
        
        try:
            response = self.llm.invoke([
                SystemMessage(content=system_prompt),
                HumanMessage(content=user_prompt)
            ])
            
            # Parse the LLM response
            import json
            content = response.content
            
            # Try to extract JSON from the response
            json_match = re.search(r'\{.*\}', content, re.DOTALL)
            if json_match:
                target_data = json.loads(json_match.group())
            else:
                # Fallback to regex parsing
                target_data = self._fallback_parse_query(query)
            
            # Convert to ProteinTarget format
            target_info: ProteinTarget = {
                "uniprot_id": target_data.get("uniprot_id"),
                "protein_name": target_data.get("protein_name"),
                "gene_name": target_data.get("gene_name"),
                "organism": target_data.get("organism", "human"),
                "sequence": None,
                "sequence_length": None,
                "function_description": None,
                "disease_associations": [target_data.get("disease_context")] if target_data.get("disease_context") else None,
                "mutations": target_data.get("mutations", [])
            }
            
            return target_info
            
        except Exception as e:
            logger.warning(f"LLM parsing failed, using fallback: {e}")
            return self._fallback_parse_query(query)
    
    def _fallback_parse_query(self, query: str) -> ProteinTarget:
        """
        Fallback query parser using regex patterns.
        
        Args:
            query: Natural language query
            
        Returns:
            Extracted target information
        """
        query_upper = query.upper()
        
        # Common gene name patterns
        gene_patterns = [
            r'\b(KRAS|BRAF|EGFR|HER2|ALK|MET|RET|ROS1|NTRK|PIK3CA|PTEN|TP53|CDK4|CDK6|VEGFR|FGFR|PDGFR)\b',
            r'\b([A-Z][A-Z0-9]{2,6})\b'  # Generic pattern for gene names
        ]
        
        # Mutation patterns
        mutation_pattern = r'([A-Z]\d+[A-Z])'  # e.g., G12C, T790M
        
        # Extract gene name
        gene_name = None
        for pattern in gene_patterns:
            match = re.search(pattern, query_upper)
            if match:
                gene_name = match.group(1)
                break
        
        # Extract mutations
        mutations = re.findall(mutation_pattern, query_upper)
        
        # Extract UniProt ID if present
        uniprot_match = re.search(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', query)
        uniprot_id = uniprot_match.group(1) if uniprot_match else None
        
        target_info: ProteinTarget = {
            "uniprot_id": uniprot_id,
            "protein_name": None,
            "gene_name": gene_name,
            "organism": "human",  # Default to human
            "sequence": None,
            "sequence_length": None,
            "function_description": None,
            "disease_associations": None,
            "mutations": mutations if mutations else None
        }
        
        return target_info
    
    def _search_uniprot(self, target_info: ProteinTarget) -> Optional[Dict[str, Any]]:
        """
        Search UniProt for protein information.
        
        Args:
            target_info: Initial target information
            
        Returns:
            Updated protein information or None
        """
        # Search by UniProt ID if available
        if target_info.get("uniprot_id"):
            protein_data = self.uniprot_client.get_protein_by_id(target_info["uniprot_id"])
            if protein_data:
                return self.uniprot_client.parse_protein_info(protein_data)
        
        # Otherwise search by gene name
        if target_info.get("gene_name"):
            results = self.uniprot_client.search_by_gene_name(
                target_info["gene_name"],
                target_info.get("organism", "human")
            )
            
            if results:
                # Take the first result (usually the canonical sequence)
                protein_data = results[0]
                parsed_info = self.uniprot_client.parse_protein_info(protein_data)
                
                # Get the full sequence
                if parsed_info.get("uniprot_id"):
                    sequence = self.uniprot_client.get_protein_sequence(parsed_info["uniprot_id"])
                    if sequence:
                        parsed_info["sequence"] = sequence
                        parsed_info["sequence_length"] = len(sequence)
                
                return parsed_info
        
        return None
    
    def _search_known_compounds(self, uniprot_id: str) -> List[KnownCompound]:
        """
        Search for known compounds targeting the protein.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            List of known compounds
        """
        compounds = []
        
        # Search ChEMBL
        chembl_compounds = self.chembl_client.search_by_target(uniprot_id, limit=50)
        
        for comp in chembl_compounds:
            if comp.get("smiles"):  # Only include if we have structure
                known_compound: KnownCompound = {
                    "compound_id": comp["compound_id"],
                    "name": comp.get("name"),
                    "smiles": comp["smiles"],
                    "binding_affinity": comp.get("binding_affinity"),
                    "affinity_unit": comp.get("affinity_unit"),
                    "source": "ChEMBL",
                    "clinical_phase": comp.get("clinical_phase"),
                    "mechanism": comp.get("activity_type")
                }
                compounds.append(known_compound)
        
        # Also search for approved drugs
        approved_drugs = self.chembl_client.search_approved_drugs(uniprot_id)
        for drug in approved_drugs:
            if drug.get("smiles") and drug["compound_id"] not in [c["compound_id"] for c in compounds]:
                known_compound: KnownCompound = {
                    "compound_id": drug["compound_id"],
                    "name": drug.get("name"),
                    "smiles": drug["smiles"],
                    "binding_affinity": drug.get("binding_affinity"),
                    "affinity_unit": drug.get("affinity_unit"),
                    "source": "ChEMBL",
                    "clinical_phase": "Approved",
                    "mechanism": drug.get("activity_type")
                }
                compounds.append(known_compound)
        
        logger.info(f"Found {len(compounds)} known compounds from ChEMBL")
        
        return compounds
