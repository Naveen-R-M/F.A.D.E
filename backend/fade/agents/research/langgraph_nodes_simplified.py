"""
Simplified LangGraph nodes for the Research module using direct RCSB search.
ENGINEERED TO FIND SMALL MOLECULE COMPLEXES ONLY.
NO FALLBACKS - Fails fast on any error.

These nodes handle target identification directly from RCSB PDB.
"""

import re
import json
from typing import Dict, Any, List, Optional

from langchain_core.messages import HumanMessage, SystemMessage
from langchain_core.output_parsers import PydanticOutputParser
from langchain_core.prompts import PromptTemplate
from langgraph.prebuilt import ToolNode

from fade.state.langgraph_state import DrugDiscoveryState
from fade.models.structured_outputs import TargetExtractionOutput
from fade.tools.rcsb_api_enhanced import get_rcsb_enhanced_client
from fade.tools.llm_client import get_llm_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.research_simplified")


def target_research_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Simplified LangGraph node for researching protein targets.
    ENGINEERED TO FIND SMALL MOLECULE COMPLEXES.
    NO FALLBACKS - Throws errors immediately on failure.
    
    This node:
    1. Parses the query to extract target information
    2. REFINES query to specifically search for small molecule complexes
    3. Searches RCSB PDB directly for structures with small molecule inhibitors
    4. FAILS if no structures with small molecules found
    
    Args:
        state: Current graph state
        
    Returns:
        Updated state dict with target information
    """
    query = state["query"]
    logger.info(f"[Research Node] Processing query: {query}")
    
    # Initialize LLM for parsing - NO FALLBACK
    llm = get_llm_client()
    
    # Parse query to extract target information - NO FALLBACK PARSING
    target_info = _parse_query_with_llm(llm, query)
    
    # Get enhanced RCSB client
    rcsb_client = get_rcsb_enhanced_client()
    
    # BUILD SMALL MOLECULE-FOCUSED SEARCH QUERY
    search_query = _build_small_molecule_search_query(target_info, query)
    
    logger.info(f"Searching RCSB for small molecule complexes with: {search_query}")
    
    # Search RCSB and get both structures and target info
    existing_structures, extracted_target_info = rcsb_client.search_by_query(
        search_query,
        limit=20  # Get more structures to find ones with small molecules
    )
    
    # NO FALLBACK - Fail if no structures found
    if not existing_structures:
        # Try alternative search with "inhibitor" keyword
        alt_query = _build_alternative_search_query(target_info)
        logger.info(f"No structures found, trying alternative query: {alt_query}")
        
        existing_structures, extracted_target_info = rcsb_client.search_by_query(
            alt_query,
            limit=20
        )
        
        if not existing_structures:
            return {
                "error": f"No PDB structures found for '{search_query}' or '{alt_query}'. Target may not have crystallized small molecule complexes.",
                "should_continue": False,
                "current_step": "research_failed"
            }
    
    # NO FALLBACK - Fail if no target info extracted
    if not extracted_target_info:
        return {
            "error": f"Could not extract target information from RCSB for '{search_query}'",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Merge extracted target info with parsed info
    target_info.update(extracted_target_info)
    logger.info(f"Extracted target info from RCSB: {target_info.get('protein_name')}")
    
    # NO FALLBACK - Require sequence for structure prediction
    if not target_info.get("sequence"):
        return {
            "error": f"No protein sequence found in PDB data for '{search_query}'",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Filter for structures with small molecules ONLY
    structures_with_small_molecules = []
    for struct in existing_structures:
        if struct.get("has_drug_like_ligand"):
            structures_with_small_molecules.append(struct)
            ligands = struct.get("ligands", [])
            drug_like = [l for l in ligands if _is_small_molecule_ligand(l)]
            if drug_like:
                logger.info(f"  {struct['pdb_id']}: Found small molecules: {[l.get('id') for l in drug_like]}")
    
    if not structures_with_small_molecules:
        return {
            "error": f"Found {len(existing_structures)} structures but none contain small molecule inhibitors. Try searching for '{target_info.get('gene_name', 'target')} kinase inhibitor complex'",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Add message to history
    message = f"Identified target: {target_info.get('protein_name', 'Unknown')} "
    
    if target_info.get("uniprot_id"):
        message += f"({target_info['uniprot_id']})"
    
    message += f"\nFound {len(structures_with_small_molecules)} PDB structures with small molecule inhibitors"
    
    # Store only structures with small molecules
    target_info["existing_structures"] = structures_with_small_molecules
    
    return {
        "target_info": target_info,
        "known_compounds": [],  # Removed ChEMBL search
        "messages": [HumanMessage(content=message)],
        "current_step": "research_complete",
        "should_continue": True
    }


def validate_target_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for validating target information.
    NO FALLBACKS - Strict validation with immediate failure.
    
    This node checks if we have enough information to proceed.
    
    Args:
        state: Current graph state
        
    Returns:
        Updated state dict
    """
    target_info = state.get("target_info")
    
    # NO FALLBACK - Must have target info
    if not target_info:
        return {
            "error": "No target information available",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # NO FALLBACK - Must have structures with small molecules
    structures = target_info.get("existing_structures", [])
    if not structures:
        return {
            "error": "No PDB structures with small molecule inhibitors available",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # Verify at least one structure has drug-like ligand
    has_small_molecule = any(s.get("has_drug_like_ligand") for s in structures)
    if not has_small_molecule:
        return {
            "error": "No structures contain small molecule drug-like compounds",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # NO FALLBACK - Must have sequence
    if not target_info.get("sequence"):
        return {
            "error": "No protein sequence available",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # Check sequence length constraint for Boltz2 (if needed later)
    sequence = target_info["sequence"]
    if len(sequence) > config.BOLTZ2_MAX_SEQUENCE_LENGTH:
        return {
            "error": f"Protein sequence too long ({len(sequence)} > {config.BOLTZ2_MAX_SEQUENCE_LENGTH})",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # Add validation message
    message = f"✓ Target validated: {target_info.get('protein_name', 'Unknown')}"
    message += f" - {len(structures)} PDB structures with small molecule inhibitors"
    message += f" - {len(sequence)} amino acids"
    
    return {
        "messages": [HumanMessage(content=message)],
        "current_step": "validation_complete",
        "should_continue": True
    }


def _parse_query_with_llm(llm, query: str) -> Dict[str, Any]:
    """
    Use LLM with structured outputs to parse the natural language query.
    ENGINEERED TO IDENTIFY SMALL MOLECULE DRUG DISCOVERY TARGETS.
    
    Args:
        llm: LLM client
        query: Natural language query
        
    Returns:
        Parsed target information as dict
        
    Raises:
        Exception: If parsing fails
    """
    # Create output parser
    parser = PydanticOutputParser(pydantic_object=TargetExtractionOutput)
    
    # Create prompt template with format instructions
    prompt = PromptTemplate(
        template="""Extract protein target information from this SMALL MOLECULE drug discovery query.
Focus on identifying targets for SMALL MOLECULE INHIBITORS, not antibodies or biologics.

Query: "{query}"

Instructions:
- For gene_name: Extract gene symbols like KRAS, EGFR, BTK, JAK2
- For protein_name: Focus on KINASES, PROTEASES, RECEPTORS, ENZYMES
- For mutations: Use standard notation like G12C, T790M, L858R
- For drug_keywords: Extract small molecule drug names (erlotinib, gefitinib, etc.)
- For target_type: Classify as kinase, protease, receptor, enzyme, or unknown
- If a generic protein name is given (like "EGFR"), assume they want the KINASE DOMAIN with INHIBITORS

{format_instructions}

Extracted Information:""",
        input_variables=["query"],
        partial_variables={"format_instructions": parser.get_format_instructions()}
    )
    
    # Format the prompt
    formatted_prompt = prompt.format(query=query)
    
    # Invoke the LLM
    response = llm.invoke(formatted_prompt)
    
    # Extract content from response
    content = response.content if hasattr(response, 'content') else str(response)
    
    # Parse the response using the Pydantic parser
    try:
        result = parser.parse(content)
        logger.info(f"Successfully parsed target using structured output: {result.gene_name or result.protein_name}")
        # Convert Pydantic model to dict for backward compatibility
        target_info = result.dict(exclude_none=True)
    except Exception as e:
        # Fallback: Try to extract JSON manually
        logger.warning(f"Pydantic parsing failed, attempting manual extraction: {e}")
        
        json_match = re.search(r'\{.*\}', content, re.DOTALL)
        if json_match:
            try:
                data = json.loads(json_match.group())
                # Validate with Pydantic model
                result = TargetExtractionOutput(**data)
                target_info = result.dict(exclude_none=True)
                logger.info(f"Successfully parsed via JSON extraction: {result.gene_name or result.protein_name}")
            except Exception as json_error:
                logger.error(f"JSON parsing also failed: {json_error}")
                raise ValueError(f"Could not parse LLM response for query: {query}. Error: {e}")
        else:
            raise ValueError(f"Could not parse LLM response for query: {query}. Error: {e}")
    
    # Validate we got something useful
    if not target_info.get("gene_name") and not target_info.get("protein_name"):
        raise ValueError(f"Could not identify any target from query: {query}")
    
    return target_info


def _build_small_molecule_search_query(target_info: Dict[str, Any], original_query: str) -> str:
    """
    Build search query SPECIFICALLY for small molecule complexes.
    
    Args:
        target_info: Parsed target information
        original_query: Original user query
        
    Returns:
        Search query string optimized for finding small molecule complexes
    """
    # Start with gene/protein name
    base_name = target_info.get("gene_name") or target_info.get("protein_name", "")
    
    # If we have specific drug keywords, use them
    drug_keywords = target_info.get("drug_keywords", [])
    if drug_keywords:
        # Search for specific drug complex
        return f"{base_name} {' '.join(drug_keywords)}"
    
    # Add mutations if specified
    mutations = target_info.get("mutations", [])
    mutation_str = " ".join(mutations) if mutations else ""
    
    # Determine target type for better search
    target_type = target_info.get("target_type", "")
    
    # Build query based on target type
    if target_type == "kinase" or "kinase" in base_name.lower():
        # For kinases, specifically look for inhibitor complexes
        return f"{base_name} {mutation_str} kinase inhibitor complex".strip()
    
    elif target_type == "protease":
        return f"{base_name} {mutation_str} protease inhibitor".strip()
    
    else:
        # Generic search for small molecule complexes
        # Add multiple keywords to increase chances of finding small molecule structures
        return f"{base_name} {mutation_str} inhibitor antagonist ligand small molecule".strip()


def _build_alternative_search_query(target_info: Dict[str, Any]) -> str:
    """
    Build alternative search query if first search fails.
    
    Args:
        target_info: Parsed target information
        
    Returns:
        Alternative search query
    """
    base_name = target_info.get("gene_name") or target_info.get("protein_name", "")
    
    # Try with common drug names for the target
    common_drugs = _get_common_drugs_for_target(base_name)
    
    if common_drugs:
        return f"{base_name} {common_drugs[0]}"
    
    # Try with just "drug" or "compound"
    return f"{base_name} drug compound complex"


def _get_common_drugs_for_target(target: str) -> List[str]:
    """
    Get common small molecule drugs for known targets.
    
    Args:
        target: Target name
        
    Returns:
        List of common drug names
    """
    # Common target-drug mappings for fallback searches
    drug_map = {
        "EGFR": ["erlotinib", "gefitinib", "osimertinib", "afatinib"],
        "BTK": ["ibrutinib", "acalabrutinib", "zanubrutinib"],
        "JAK": ["ruxolitinib", "tofacitinib", "baricitinib"],
        "JAK2": ["ruxolitinib", "fedratinib"],
        "BCR-ABL": ["imatinib", "dasatinib", "nilotinib"],
        "ABL": ["imatinib", "dasatinib", "nilotinib"],
        "BRAF": ["vemurafenib", "dabrafenib", "encorafenib"],
        "MEK": ["trametinib", "cobimetinib", "binimetinib"],
        "CDK4": ["palbociclib", "ribociclib", "abemaciclib"],
        "CDK6": ["palbociclib", "ribociclib", "abemaciclib"],
        "HER2": ["lapatinib", "neratinib", "tucatinib"],
        "KRAS": ["sotorasib", "adagrasib"],
        "ALK": ["crizotinib", "alectinib", "lorlatinib"],
        "ROS1": ["crizotinib", "entrectinib"],
        "MET": ["crizotinib", "cabozantinib", "capmatinib"],
        "VEGFR": ["sorafenib", "sunitinib", "pazopanib"],
        "FGFR": ["erdafitinib", "pemigatinib"],
        "PI3K": ["idelalisib", "alpelisib", "copanlisib"],
        "mTOR": ["sirolimus", "everolimus", "temsirolimus"],
        "PARP": ["olaparib", "niraparib", "rucaparib"],
    }
    
    target_upper = target.upper()
    for key, drugs in drug_map.items():
        if key in target_upper:
            return drugs
    
    return []


def _is_small_molecule_ligand(ligand: Dict[str, Any]) -> bool:
    """
    Check if a ligand is a small molecule drug-like compound.
    More specific than general drug-like check.
    
    Args:
        ligand: Ligand information
        
    Returns:
        True if ligand is a small molecule
    """
    mw = ligand.get("molecular_weight", 0)
    if isinstance(mw, list):
        mw = mw[0] if mw else 0
    
    # Small molecule criteria: MW between 150-800 Da
    if not (150 < mw < 800):
        return False
    
    # Check it's not a known non-drug artifact
    ligand_id = ligand.get("id", "").upper()
    if ligand_id in config.PDB_CRYSTALLIZATION_ARTIFACTS:
        return False
    if ligand_id in config.PDB_NUCLEOTIDES:
        return False
    
    return True
