"""
LangGraph nodes for the Research module.

These nodes handle target identification and known compound research.
"""

import re
import json
from typing import Dict, Any, List, Optional

from langchain.chat_models import init_chat_model
from langchain_core.messages import HumanMessage, SystemMessage
from langgraph.prebuilt import ToolNode

from fade.state.langgraph_state import DrugDiscoveryState
from fade.tools import get_uniprot_client, get_chembl_client, get_rcsb_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.research")


def target_research_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for researching protein targets.
    
    This node:
    1. Parses the query to extract target information
    2. Searches UniProt for protein details
    3. Finds known compounds from ChEMBL
    4. Checks for existing structures in PDB
    
    Args:
        state: Current graph state
        
    Returns:
        Updated state dict with target information
    """
    query = state["query"]
    logger.info(f"[Research Node] Processing query: {query}")
    
    # Initialize LLM for parsing
    llm_config = config.get_llm_config()
    llm = init_chat_model(**llm_config)
    
    # Parse query to extract target information
    target_info = _parse_query_with_llm(llm, query)
    
    if not target_info.get("gene_name") and not target_info.get("uniprot_id"):
        return {
            "error": "Could not identify protein target from query",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Search UniProt
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
    
    # Search for known compounds
    known_compounds = []
    if target_info.get("uniprot_id"):
        chembl_client = get_chembl_client()
        known_compounds = _search_known_compounds(chembl_client, target_info["uniprot_id"])
        logger.info(f"Found {len(known_compounds)} known compounds")
    
    # Check for existing structures
    existing_structures = []
    if target_info.get("uniprot_id"):
        rcsb_client = get_rcsb_client()
        existing_structures = rcsb_client.search_by_uniprot(
            target_info["uniprot_id"],
            limit=5
        )
        logger.info(f"Found {len(existing_structures)} existing PDB structures")
    
    # Add message to history
    message = f"Identified target: {target_info.get('protein_name')} ({target_info.get('uniprot_id')})"
    if known_compounds:
        message += f"\nFound {len(known_compounds)} known compounds"
    if existing_structures:
        message += f"\nFound {len(existing_structures)} existing structures"
    
    # Store existing structures for structure module to use
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
    
    This node checks if we have enough information to proceed.
    
    Args:
        state: Current graph state
        
    Returns:
        Updated state dict
    """
    target_info = state.get("target_info", {})
    
    # Check minimum requirements
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
    Parse natural language query using LLM.
    
    Args:
        llm: Language model instance
        query: Natural language query
        
    Returns:
        Extracted target information
    """
    system_prompt = """You are a bioinformatics expert. Extract protein target information from the drug discovery query.
    
    Return ONLY a JSON object with these fields:
    {
        "gene_name": "gene symbol like KRAS, EGFR",
        "protein_name": "protein name if mentioned",
        "uniprot_id": "UniProt ID if mentioned like P01112",
        "organism": "organism, default to human",
        "mutations": ["list", "of", "mutations", "like", "G12C"],
        "disease_context": "disease context if mentioned"
    }
    
    Return ONLY valid JSON, no other text."""
    
    try:
        response = llm.invoke([
            SystemMessage(content=system_prompt),
            HumanMessage(content=f"Extract target from: {query}")
        ])
        
        # Parse JSON from response
        content = response.content.strip()
        
        # Try to extract JSON
        if content.startswith('{'):
            return json.loads(content)
        else:
            # Try to find JSON in the response
            json_match = re.search(r'\{.*\}', content, re.DOTALL)
            if json_match:
                return json.loads(json_match.group())
    except Exception as e:
        logger.warning(f"LLM parsing failed: {e}")
    
    # Fallback to regex parsing
    return _fallback_parse_query(query)


def _fallback_parse_query(query: str) -> Dict[str, Any]:
    """
    Fallback query parser using regex.
    
    Args:
        query: Natural language query
        
    Returns:
        Extracted target information
    """
    query_upper = query.upper()
    
    # Common gene names
    gene_pattern = r'\b(KRAS|BRAF|EGFR|HER2|ALK|MET|RET|ROS1|NTRK|PIK3CA|PTEN|TP53|CDK4|CDK6|VEGFR|FGFR|PDGFR)\b'
    gene_match = re.search(gene_pattern, query_upper)
    
    # Mutations
    mutations = re.findall(r'([A-Z]\d+[A-Z])', query_upper)
    
    return {
        "gene_name": gene_match.group(1) if gene_match else None,
        "mutations": mutations if mutations else [],
        "organism": "human"
    }


def _search_uniprot(uniprot_client, target_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Search UniProt for protein information.
    
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
    
    return None


def _search_known_compounds(chembl_client, uniprot_id: str) -> List[Dict[str, Any]]:
    """
    Search for known compounds targeting the protein.
    
    Args:
        chembl_client: ChEMBL API client
        uniprot_id: UniProt accession ID
        
    Returns:
        List of known compounds
    """
    compounds = []
    
    # Search ChEMBL
    chembl_compounds = chembl_client.search_by_target(uniprot_id, limit=50)
    
    for comp in chembl_compounds:
        if comp.get("smiles"):
            compound = {
                "compound_id": comp["compound_id"],
                "name": comp.get("name"),
                "smiles": comp["smiles"],
                "binding_affinity": comp.get("binding_affinity"),
                "affinity_unit": comp.get("affinity_unit"),
                "source": "ChEMBL",
                "clinical_phase": comp.get("clinical_phase"),
                "mechanism": comp.get("activity_type")
            }
            compounds.append(compound)
    
    # Also get approved drugs
    approved_drugs = chembl_client.search_approved_drugs(uniprot_id)
    for drug in approved_drugs:
        if drug.get("smiles") and drug["compound_id"] not in [c["compound_id"] for c in compounds]:
            compounds.append({
                "compound_id": drug["compound_id"],
                "name": drug.get("name"),
                "smiles": drug["smiles"],
                "binding_affinity": drug.get("binding_affinity"),
                "affinity_unit": drug.get("affinity_unit"),
                "source": "ChEMBL",
                "clinical_phase": "Approved"
            })
    
    return compounds
