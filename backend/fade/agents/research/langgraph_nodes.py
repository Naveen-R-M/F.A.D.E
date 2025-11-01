"""
LangGraph nodes for the Research module.
NO FALLBACKS - Fails immediately on any error.

These nodes handle target identification using UniProt-first approach.
(Note: Consider using langgraph_nodes_simplified.py for direct RCSB search)
"""

import re
import json
from typing import Dict, Any, List, Optional

from langchain_core.messages import HumanMessage, SystemMessage
from langchain_core.output_parsers import PydanticOutputParser
from langchain_core.prompts import PromptTemplate, ChatPromptTemplate
from langgraph.prebuilt import ToolNode

from fade.state.langgraph_state import DrugDiscoveryState
from fade.models.structured_outputs import TargetExtractionOutput
from fade.tools import get_uniprot_client, get_rcsb_client
from fade.tools.llm_client import get_llm_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.research")


def target_research_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for researching protein targets.
    NO FALLBACKS - Fails immediately on any error.
    
    This node:
    1. Parses the query to extract target information
    2. Searches UniProt for protein details
    3. Checks for existing structures in PDB
    
    Args:
        state: Current graph state
        
    Returns:
        Updated state dict with target information
    """
    query = state["query"]
    logger.info(f"[Research Node] Processing query: {query}")
    
    # Initialize LLM for parsing - Will raise if fails
    llm = get_llm_client()
    
    # Parse query to extract target information - NO FALLBACK
    target_info = _parse_query_with_llm(llm, query)
    
    if not target_info.get("gene_name") and not target_info.get("uniprot_id"):
        return {
            "error": "Could not identify protein target from query",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Search UniProt - NO ERROR HANDLING
    uniprot_client = get_uniprot_client()
    protein_data = _search_uniprot(uniprot_client, target_info)
    
    if not protein_data:
        return {
            "error": f"Could not find protein information for {target_info.get('gene_name') or target_info.get('uniprot_id')}",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Merge protein data with target info
    target_info.update(protein_data)
    
    # Known compounds removed - ChEMBL integration deprecated
    known_compounds = []
    
    # Check for existing structures with enhanced ligand information
    existing_structures = []
    if target_info.get("uniprot_id"):
        rcsb_client = get_rcsb_client()
        # Use the enhanced search that includes ligand information
        existing_structures = rcsb_client.search_by_uniprot(
            target_info["uniprot_id"],
            limit=10,  # Get more initially for better filtering
            with_ligands=True  # Enable ligand information retrieval
        )
        logger.info(f"Found {len(existing_structures)} existing PDB structures")
        
        # Log structures with drug-like ligands
        drug_like_count = sum(1 for s in existing_structures if s.get("has_drug_like_ligand"))
        if drug_like_count > 0:
            logger.info(f"  {drug_like_count} structures have drug-like ligands")
            for s in existing_structures:
                if s.get("has_drug_like_ligand"):
                    ligands = s.get("ligands", [])
                    ligand_names = [l.get("id") for l in ligands if l.get("molecular_weight", 0) > 150]
                    logger.debug(f"  {s['pdb_id']}: {', '.join(ligand_names)}")
    
    # Add message to history
    message = f"Identified target: {target_info.get('protein_name')} ({target_info.get('uniprot_id')})"
    if existing_structures:
        message += f"\nFound {len(existing_structures)} existing structures"
    
    # Store existing structures with ligand info for structure module to use
    if existing_structures:
        target_info["existing_structures"] = existing_structures
    
    return {
        "target_info": target_info,
        "known_compounds": known_compounds,
        "messages": [HumanMessage(content=message)],
        "current_step": "research_complete",
        "should_continue": True
    }


def validate_target_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for validating target information.
    NO FALLBACKS - Strict validation.
    
    This node checks if we have enough information to proceed.
    
    Args:
        state: Current graph state
        
    Returns:
        Updated state dict
    """
    target_info = state.get("target_info", {})
    
    # Check minimum requirements - NO DEFAULTS
    if not target_info.get("sequence"):
        return {
            "error": "No protein sequence found for target",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    if len(target_info["sequence"]) > config.BOLTZ2_MAX_SEQUENCE_LENGTH:
        return {
            "error": f"Protein sequence too long ({len(target_info['sequence'])} > {config.BOLTZ2_MAX_SEQUENCE_LENGTH})",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # Add validation message
    message = f"✓ Target validated: {target_info.get('protein_name')} - {len(target_info['sequence'])} amino acids"
    
    return {
        "messages": [HumanMessage(content=message)],
        "current_step": "target_validated",
        "should_continue": True
    }


def _parse_query_with_llm(llm, query: str) -> Dict[str, Any]:
    """
    Parse natural language query using LLM with structured outputs.
    
    Args:
        llm: Language model instance
        query: Natural language query
        
    Returns:
        Extracted target information as dict
        
    Raises:
        ValueError: If parsing fails
    """
    # Create output parser
    parser = PydanticOutputParser(pydantic_object=TargetExtractionOutput)
    
    # Create prompt template
    prompt = PromptTemplate(
        template="""You are a bioinformatics expert. Extract protein target information from this drug discovery query.

Query: {query}

Focus on:
- Gene symbols (KRAS, EGFR, BTK, JAK2)
- Protein names (kinases, proteases, receptors)
- UniProt IDs if mentioned (like P01112)
- Specific mutations (G12C, T790M, L858R)
- Disease context and requirements

{format_instructions}

Extracted Information:""",
        input_variables=["query"],
        partial_variables={"format_instructions": parser.get_format_instructions()}
    )
    
    # Format the prompt
    formatted_prompt = prompt.format(query=query)
    
    # Invoke the LLM
    response = llm.invoke(formatted_prompt)
    
    # Extract content
    content = response.content if hasattr(response, 'content') else str(response)
    
    # Parse with Pydantic
    try:
        result = parser.parse(content)
        logger.info(f"Successfully parsed target: {result.gene_name or result.protein_name}")
        # Convert to dict for backward compatibility
        target_info = result.dict(exclude_none=True)
        
        # Add organism default if not present (for backward compatibility)
        if 'organism' not in target_info:
            target_info['organism'] = 'Homo sapiens'
            
        return target_info
        
    except Exception as e:
        # Fallback: Try JSON extraction
        logger.warning(f"Structured parsing failed, trying JSON extraction: {e}")
        json_match = re.search(r'\{.*\}', content, re.DOTALL)
        if json_match:
            try:
                data = json.loads(json_match.group())
                # Validate with model
                result = TargetExtractionOutput(**data)
                target_info = result.dict(exclude_none=True)
                if 'organism' not in target_info:
                    target_info['organism'] = 'Homo sapiens'
                return target_info
            except Exception as json_error:
                logger.error(f"JSON parsing also failed: {json_error}")
                raise ValueError(f"Failed to parse query: {query}")
        else:
            raise ValueError(f"No valid output found in LLM response")


def _search_uniprot(uniprot_client, target_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Search UniProt for protein information.
    NO FALLBACK - Returns None if not found.
    
    Args:
        uniprot_client: UniProt API client
        target_info: Target information
        
    Returns:
        Protein data or None
    """
    # Search by UniProt ID first
    if target_info.get("uniprot_id"):
        protein_data = uniprot_client.get_protein_by_id(target_info["uniprot_id"])
        if protein_data:
            parsed = uniprot_client.parse_protein_info(protein_data)
            sequence = uniprot_client.get_protein_sequence(parsed["uniprot_id"])
            if sequence:
                parsed["sequence"] = sequence
                parsed["sequence_length"] = len(sequence)
            return parsed
    
    # Search by gene name
    if target_info.get("gene_name"):
        results = uniprot_client.search_by_gene_name(
            target_info["gene_name"],
            target_info.get("organism", "human")
        )
        
        if results:
            protein_data = results[0]
            parsed = uniprot_client.parse_protein_info(protein_data)
            sequence = uniprot_client.get_protein_sequence(parsed["uniprot_id"])
            if sequence:
                parsed["sequence"] = sequence
                parsed["sequence_length"] = len(sequence)
            return parsed
    
    # NO FALLBACK - Return None if not found
    return None


# ChEMBL compound search removed - deprecated functionality
