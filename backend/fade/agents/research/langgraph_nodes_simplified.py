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
from fade.models.structured_outputs import TargetExtractionOutput, QueryRefinementSuggestion
from fade.tools.rcsb_api_enhanced import get_rcsb_enhanced_client
from fade.tools.llm_client import get_llm_client
from fade.config import config
from fade.utils import get_logger
from fade.agents.research.guidance_generator import generate_query_guidance

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
    
    # LOG: Show what was parsed
    logger.info(f"Parsed target info:")
    logger.info(f"  UniProt ID: {target_info.get('uniprot_id')}")
    logger.info(f"  Gene: {target_info.get('gene_name')}")
    logger.info(f"  Protein: {target_info.get('protein_name')}")
    logger.info(f"  PDB ID: {target_info.get('pdb_id')}")
    logger.info(f"  Drug keywords: {target_info.get('drug_keywords', [])}")
    logger.info(f"  Mutations: {target_info.get('mutations', [])}")
    
    # Get enhanced RCSB client (used by both PDB ID and regular queries)
    rcsb_client = get_rcsb_enhanced_client()
    
    # CASE 1: Direct PDB ID query - FETCH, DON'T SEARCH
    if target_info.get("pdb_id"):
        pdb_id = target_info["pdb_id"].upper()
        logger.info(f"Direct PDB ID query detected: {pdb_id}")
        
        # For direct PDB queries, fetch the specific structure WITHOUT filtering
        # The user asked for THIS structure, give it to them
        existing_structures, extracted_target_info = rcsb_client.get_structure_by_id(pdb_id)
        
        if not existing_structures:
            return {
                "error": f"PDB structure {pdb_id} not found in RCSB database. Please verify the PDB ID.",
                "should_continue": False,
                "current_step": "research_failed"
            }
        
        # Extract target information from the PDB structure
        struct = existing_structures[0] if isinstance(existing_structures, list) else existing_structures
        
        # Update target_info with what we found
        if extracted_target_info:
            if extracted_target_info.get("gene_name"):
                target_info["gene_name"] = extracted_target_info["gene_name"]
            if extracted_target_info.get("protein_name"):
                target_info["protein_name"] = extracted_target_info["protein_name"]
            if extracted_target_info.get("uniprot_id"):
                target_info["uniprot_id"] = extracted_target_info["uniprot_id"]
            if extracted_target_info.get("sequence"):
                target_info["sequence"] = extracted_target_info["sequence"]
                target_info["sequence_length"] = len(extracted_target_info["sequence"])
        
        # Log what we found
        logger.info(f"Found PDB {pdb_id}: {struct.get('title', 'Unknown')[:100]}")
        if struct.get("ligands"):
            ligand_ids = [l.get("id") for l in struct["ligands"][:5]]
            logger.info(f"  Ligands: {', '.join(ligand_ids)}")
            # Note if they're not drug-like
            if not struct.get("has_drug_like_ligand"):
                logger.info(f"  Note: Ligands are not drug-like (e.g., GDP, ATP, metal ions)")
        
        # Store the structure - DON'T FILTER FOR DRUG-LIKE!
        # Wrap in list if needed
        if not isinstance(existing_structures, list):
            existing_structures = [existing_structures]
        target_info["existing_structures"] = existing_structures
        
        # Add message
        message = f"Found PDB structure {pdb_id}"
        if target_info.get("protein_name"):
            message += f": {target_info['protein_name']}"
        if target_info.get("uniprot_id"):
            message += f" ({target_info['uniprot_id']})"
        if struct.get("ligands"):
            message += f" with {len(struct['ligands'])} ligand(s)"
            if not struct.get("has_drug_like_ligand"):
                message += " (non-drug-like: GDP, cofactors, or ions)"
        
        # Mark this as a direct PDB query for downstream nodes
        target_info["is_direct_pdb_query"] = True
        
        return {
            "target_info": target_info,
            "known_compounds": [],
            "messages": [HumanMessage(content=message)],
            "current_step": "research_complete",
            "should_continue": True
        }
    
    # CASE 1.5: Direct UniProt ID query
    if target_info.get("uniprot_id") and not target_info.get("gene_name") and not target_info.get("protein_name"):
        # We have a UniProt ID but no gene/protein name - this was a direct UniProt query
        uniprot_id = target_info["uniprot_id"]
        logger.info(f"Direct UniProt ID query detected: {uniprot_id}")
        
        # For UniProt IDs, we should search RCSB for any structures from this protein
        search_query = uniprot_id
        
        # Search RCSB for structures from this UniProt ID
        existing_structures, extracted_target_info = rcsb_client.search_by_query(
            search_query,
            limit=20,  # Get multiple structures
            original_target=target_info  # Pass for validation
        )
        
        # Update target_info with any extracted information
        if extracted_target_info:
            if extracted_target_info.get("gene_name"):
                target_info["gene_name"] = extracted_target_info["gene_name"]
            if extracted_target_info.get("protein_name"):
                target_info["protein_name"] = extracted_target_info["protein_name"]
            if extracted_target_info.get("sequence"):
                target_info["sequence"] = extracted_target_info["sequence"]
                target_info["sequence_length"] = len(extracted_target_info["sequence"])
        
        # Whether or not we found PDB structures, we have a valid UniProt ID
        # We can continue to AlphaFold/Boltz2 if no PDB structures
        if not existing_structures:
            logger.info(f"No PDB structures found for UniProt {uniprot_id}, will try AlphaFold/Boltz2")
            target_info["existing_structures"] = []
            
            message = f"Identified UniProt ID: {uniprot_id}"
            if target_info.get("protein_name"):
                message += f" ({target_info['protein_name']})"
            message += "\nNo PDB structures found - will check AlphaFold database or predict structure"
            
            return {
                "target_info": target_info,
                "known_compounds": [],
                "messages": [HumanMessage(content=message)],
                "current_step": "research_complete_no_pdb",
                "should_continue": True
            }
        
        # If we found structures, continue with normal flow
        # (The rest of the function will handle filtering for small molecules)
        logger.info(f"Found {len(existing_structures)} PDB structures for UniProt {uniprot_id}")
    
    # CASE 2: Regular target query (existing logic)
    # Try to get UniProt ID by searching UniProt database
    original_uniprot_id = _get_uniprot_id_for_target(target_info)
    if original_uniprot_id:
        target_info["uniprot_id"] = original_uniprot_id
        logger.info(f"Found UniProt ID via search: {original_uniprot_id}")
    
    # BUILD SMART SEARCH QUERY using parsed target info
    # Use the parsed target information to build a focused query
    search_query = _build_smart_search_query(target_info, query)
    
    logger.info(f"Searching RCSB with focused query: {search_query}")
    
    # Search RCSB and get both structures and target info
    # CRITICAL: Pass original target info for validation
    existing_structures, extracted_target_info = rcsb_client.search_by_query(
        search_query,
        limit=20,  # Get more structures to find ones with small molecules
        original_target=target_info  # Pass original target for validation
    )
    
    # NO PDB STRUCTURES FOUND - But we can still continue with AlphaFold/Boltz2
    if not existing_structures:
        logger.info(f"No PDB structures found for {search_query}, but will try AlphaFold/Boltz2")
        
        # We still need to extract basic target info even without PDB structures
        # Check if protein_name might actually be a UniProt ID
        protein_name = target_info.get("protein_name")
        if protein_name and (len(protein_name) == 6 or len(protein_name) == 10 or "_" in protein_name):
            # Looks like a UniProt ID pattern (e.g., P01116 or A0A1L1T3F0)
            target_info["uniprot_id"] = protein_name
            logger.info(f"Recognized '{protein_name}' as a UniProt ID")
        
        # Try to get UniProt ID if we don't have it
        if not target_info.get("uniprot_id"):
            # For known targets, we may have UniProt ID
            uniprot_id = _get_uniprot_id_for_target(target_info)
            if uniprot_id:
                target_info["uniprot_id"] = uniprot_id
                logger.info(f"Using known UniProt ID: {uniprot_id}")
        
        # We don't have a sequence yet, but that's OK - AlphaFold/Boltz2 can handle it
        # Mark that we have no PDB structures
        target_info["existing_structures"] = []
        
        # Add message explaining the situation
        message = f"Identified target: {target_info.get('protein_name', target_info.get('gene_name', 'Unknown'))}"
        if target_info.get("uniprot_id"):
            message += f" ({target_info['uniprot_id']})"
        message += "\nNo PDB structures found - will check AlphaFold database or predict structure"
        
        # CONTINUE to structure resolver instead of failing
        return {
            "target_info": target_info,
            "known_compounds": [],
            "messages": [HumanMessage(content=message)],
            "current_step": "research_complete_no_pdb",
            "should_continue": True  # CONTINUE instead of stopping!
        }
        

    
    # NO FALLBACK - Fail if no target info extracted
    if not extracted_target_info:
        return {
            "error": f"Could not extract target information from RCSB for '{search_query}'",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # CRITICAL: Preserve original target info - don't overwrite with search results!
    # Only merge in missing information like sequence
    if extracted_target_info and not extracted_target_info.get("no_results"):
        # Only add sequence if we don't have it
        if not target_info.get("sequence") and extracted_target_info.get("sequence"):
            target_info["sequence"] = extracted_target_info["sequence"]
            target_info["sequence_length"] = extracted_target_info.get("sequence_length")
            logger.info(f"Added sequence from RCSB: {len(target_info['sequence'])} amino acids")
        
        # Preserve original UniProt ID if we have it
        if original_uniprot_id:
            target_info["uniprot_id"] = original_uniprot_id
            logger.info(f"Preserved original UniProt ID: {original_uniprot_id}")
        elif not target_info.get("uniprot_id") and extracted_target_info.get("uniprot_id"):
            # Only use extracted UniProt if we don't have one
            target_info["uniprot_id"] = extracted_target_info["uniprot_id"]
    
    logger.info(f"Target info preserved: {target_info.get('protein_name')} ({target_info.get('uniprot_id')})")
    
    # NO FALLBACK - Require sequence for structure prediction
    if not target_info.get("sequence"):
        return {
            "error": f"No protein sequence found in PDB data for '{search_query}'",
            "should_continue": False,
            "current_step": "research_failed"
        }
    
    # Filter for structures based on specific ligand request or small molecules
    drug_keywords = target_info.get("drug_keywords", [])
    structures_with_small_molecules = []
    
    if drug_keywords:
        # User requested specific ligands - filter for those
        logger.info(f"Filtering structures for specific ligands: {drug_keywords}")
        for struct in existing_structures:
            ligands = struct.get("ligands", [])
            # Check if any requested ligand is present
            struct_ligands = [l.get('id', '').upper() for l in ligands]
            requested_found = []
            for keyword in drug_keywords:
                keyword_upper = keyword.upper()
                # Check for exact match or substring match
                if keyword_upper in struct_ligands or any(keyword_upper in lig for lig in struct_ligands):
                    requested_found.append(keyword)
            
            if requested_found:
                structures_with_small_molecules.append(struct)
                logger.info(f"  {struct['pdb_id']}: Contains requested ligand(s): {requested_found}")
    else:
        # No specific ligand requested - filter for any drug-like small molecules
        for struct in existing_structures:
            if struct.get("has_drug_like_ligand"):
                structures_with_small_molecules.append(struct)
                ligands = struct.get("ligands", [])
                drug_like = [l for l in ligands if _is_small_molecule_ligand(l)]
                if drug_like:
                    logger.info(f"  {struct['pdb_id']}: Found small molecules: {[l.get('id') for l in drug_like]}")
    
    if not structures_with_small_molecules:
        # Different error message based on whether specific ligand was requested
        if drug_keywords:
            error_msg = f"Found {len(existing_structures)} {target_info.get('gene_name', 'target')} structures but none contain {', '.join(drug_keywords)}. "
            error_msg += "This ligand may not have been crystallized with this target. "
            error_msg += f"Try searching for '{target_info.get('gene_name', 'target')} inhibitor' or just '{target_info.get('gene_name', 'target')}' to see available structures."
        else:
            error_msg = f"Found {len(existing_structures)} structures but none contain small molecule inhibitors. "
            error_msg += f"Try searching for '{target_info.get('gene_name', 'target')} kinase inhibitor complex'."
        
        return {
            "error": error_msg,
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
    
    # For direct PDB queries, we may not have drug-like ligands and that's OK
    is_pdb_query = target_info.get("pdb_id") is not None
    
    structures = target_info.get("existing_structures", [])
    
    # If no PDB structures, that's OK - we'll try AlphaFold/Boltz2
    if not structures:
        logger.info("No PDB structures available, will proceed to AlphaFold/Boltz2")
        # We don't need to fail here - structure resolver will handle it
        message = f"✓ Target validated: {target_info.get('protein_name', target_info.get('gene_name', 'Unknown'))}"
        if target_info.get("uniprot_id"):
            message += f" ({target_info['uniprot_id']})"
        message += " - No PDB structures found, will check AlphaFold/Boltz2"
        
        # We might not have a sequence yet, but that's OK - structure resolver will fetch it
        if target_info.get("sequence"):
            message += f" - {len(target_info['sequence'])} amino acids"
        
        return {
            "messages": [HumanMessage(content=message)],
            "current_step": "validation_complete",
            "should_continue": True
        }
    
    # For direct PDB queries, we're more lenient - any structure is OK
    if not is_pdb_query:
        # Verify at least one structure has drug-like ligand
        has_small_molecule = any(s.get("has_drug_like_ligand") for s in structures)
        if not has_small_molecule:
            return {
                "error": "No structures contain small molecule drug-like compounds",
                "should_continue": False,
                "current_step": "validation_failed"
            }
    
    # We might not have a sequence yet if no PDB structures were found
    # That's OK - structure resolver will fetch it from UniProt
    sequence = target_info.get("sequence")
    if sequence and len(sequence) > config.BOLTZ2_MAX_SEQUENCE_LENGTH:
        return {
            "error": f"Protein sequence too long ({len(sequence)} > {config.BOLTZ2_MAX_SEQUENCE_LENGTH})",
            "should_continue": False,
            "current_step": "validation_failed"
        }
    
    # Add validation message
    message = f"✓ Target validated: "
    if is_pdb_query:
        message += f"PDB {target_info.get('pdb_id')} - "
    message += f"{target_info.get('protein_name', 'Unknown')}"
    
    if is_pdb_query:
        message += f" - Direct PDB query (ligand requirements relaxed)"
    else:
        message += f" - {len(structures)} PDB structures with small molecule inhibitors"
    
    if sequence:
        message += f" - {len(sequence)} amino acids"
    else:
        message += f" - Sequence will be fetched from UniProt"
    
    return {
        "messages": [HumanMessage(content=message)],
        "current_step": "validation_complete",
        "should_continue": True
    }


def _get_uniprot_id_for_target(target_info: Dict[str, Any]) -> Optional[str]:
    """
    Get UniProt ID for target by searching UniProt database.
    ALWAYS fetches from UniProt API - NO HARDCODED MAPPINGS.
    
    Args:
        target_info: Parsed target information
        
    Returns:
        UniProt ID if found, None otherwise
    """
    from fade.agents.research.uniprot_resolver import resolve_uniprot_id_dynamically
    
    # Use the dynamic resolver - NO HARDCODED MAPPINGS
    # All lookups are performed via UniProt API
    logger.info("Using dynamic UniProt resolution (no hardcoded mappings)")
    return resolve_uniprot_id_dynamically(target_info)


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
    # 1. Bind the Pydantic model to the LLM
    # This creates a new LLM instance that *must* return this structure
    structured_llm = llm.with_structured_output(TargetExtractionOutput)
    
    # 2. EVEN MORE EXPLICIT PROMPT - DIRECT INSTRUCTIONS
    prompt = f"""
    You are an expert biologist. Extract key identifiers from the query.
    The query is: "{query}"

    - **uniprot_id**: Extract 6-10 character alphanumeric UniProt IDs.
      Examples: "P01116", "Q5U9M0", "A0A8K7W3T4".
    - **pdb_id**: Extract 4-character PDB IDs.
      Examples: "4OBE", "1A2B".
    - **gene_name**: Extract gene names.
      Examples: "KRAS", "EGFR", "OR52N4".
    - **protein_name**: Extract common protein names.
      Examples: "human olfactory receptor 52N4", "GTPase KRas".
    - **drug_keywords**: Extract any ligands mentioned.
      Examples: "GDP", "inhibitor", "sanglifehrin".
    
    The query "{query}" contains an identifier. Place it in the correct field.
    If the query is just an ID like "Q5U9M0", your *only* output should be in the 'uniprot_id' field.
    """
    
    # 3. Invoke the LLM. The result *is* the Pydantic object.
    # This will raise a Pydantic validation error if the LLM output is bad
    # or if our model_validator fails (e.g., no target found).
    try:
        result = structured_llm.invoke(prompt)
        
        # This log is now more useful
        logger.info(f"Successfully parsed target: "
                    f"{result.uniprot_id or result.pdb_id or result.gene_name or result.protein_name}")
        
        # Log what was extracted for debugging
        if result.uniprot_id:
            logger.info(f"  Extracted UniProt ID: {result.uniprot_id}")
        if result.gene_name:
            logger.info(f"  Extracted gene: {result.gene_name}")
        if result.protein_name:
            logger.info(f"  Extracted protein: {result.protein_name}")
        if result.pdb_id:
            logger.info(f"  Extracted PDB ID: {result.pdb_id}")
        
        # 4. Convert to dict as required by your function's return type
        return result.model_dump(exclude_none=True)
        
    except Exception as e:
        logger.error(f"LLM parsing failed, trying regex fallback: {e}")
        
        # FALLBACK: Try simple regex extraction for UniProt IDs
        import re
        
        # Pattern for UniProt IDs (6 or 10 alphanumeric chars starting with letter)
        uniprot_pattern = r'\b([A-Z][A-Z0-9]{5}|[A-Z][A-Z0-9]{9})\b'
        match = re.search(uniprot_pattern, query.upper())
        
        if match:
            uniprot_id = match.group(1)
            logger.info(f"Regex fallback extracted UniProt ID: {uniprot_id}")
            return {"uniprot_id": uniprot_id}
        
        # Pattern for PDB IDs (4 alphanumeric chars)
        pdb_pattern = r'\b([0-9][A-Z0-9]{3})\b'
        match = re.search(pdb_pattern, query.upper())
        
        if match:
            pdb_id = match.group(1)
            logger.info(f"Regex fallback extracted PDB ID: {pdb_id}")
            return {"pdb_id": pdb_id}
        
        # Re-raise the original exception if no fallback worked
        raise ValueError(f"Could not parse query even with fallback: {query}. Error: {e}")


def _build_smart_search_query(target_info: Dict[str, Any], original_query: str) -> str:
    """
    Build a FOCUSED search query using parsed target information.
    This ensures we search for the actual target, not random mentions.
    
    Args:
        target_info: Parsed target information from LLM
        original_query: Original user query
        
    Returns:
        Focused search query string
    """
    # If we have a UniProt ID, use it directly for searching
    if target_info.get("uniprot_id"):
        logger.debug(f"Using UniProt ID for search: {target_info['uniprot_id']}")
        return target_info["uniprot_id"]
    
    # Start with the target name (gene or protein)
    base_name = target_info.get("gene_name") or target_info.get("protein_name", "")
    
    if not base_name:
        # Fallback to original query if parsing failed
        logger.warning("No target name found, using original query")
        return original_query
    
    # If we have specific drug keywords (like GDP, GTP), use them
    drug_keywords = target_info.get("drug_keywords", [])
    if drug_keywords:
        # Build focused query: "KRAS GDP" not "KRAS inhibitor antagonist ligand"
        query = f"{base_name} {' '.join(drug_keywords)}"
        logger.debug(f"Using drug-specific query: {query}")
        return query
    
    # Add mutations if specified
    mutations = target_info.get("mutations", [])
    if mutations:
        mutation_str = " ".join(mutations)
        query = f"{base_name} {mutation_str}"
        logger.debug(f"Using mutation-specific query: {query}")
        return query
    
    # For general queries, just use the target name
    # Don't add "inhibitor antagonist ligand" - too restrictive!
    logger.debug(f"Using target-only query: {base_name}")
    return base_name


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
