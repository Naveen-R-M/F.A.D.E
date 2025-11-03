"""
Dynamic UniProt ID resolver - NO HARDCODED MAPPINGS.
All lookups are performed via UniProt API in real-time.
"""

from typing import Optional, Dict, Any
from fade.tools.uniprot_api import get_uniprot_client
from fade.utils import get_logger

logger = get_logger("research.uniprot_resolver")


def resolve_uniprot_id_dynamically(target_info: Dict[str, Any]) -> Optional[str]:
    """
    Resolve UniProt ID dynamically via API - NO HARDCODED MAPPINGS.
    
    Priority:
    1. If already have uniprot_id, return it
    2. If protein_name looks like UniProt ID, validate and return
    3. Search by gene_name if available
    4. Search by protein_name as last resort
    
    Args:
        target_info: Target information dict
        
    Returns:
        UniProt ID if found, None otherwise
    """
    uniprot_client = get_uniprot_client()
    
    # Case 1: Already have UniProt ID
    if target_info.get("uniprot_id"):
        return target_info["uniprot_id"]
    
    # Case 2: Check if protein_name is actually a UniProt ID
    protein_name = target_info.get("protein_name", "").strip() if target_info.get("protein_name") else ""
    if is_uniprot_id_pattern(protein_name):
        logger.info(f"Protein name '{protein_name}' matches UniProt ID pattern")
        # Validate it exists in UniProt
        if validate_uniprot_id(protein_name):
            logger.info(f"Validated '{protein_name}' as valid UniProt ID via API")
            return protein_name
        else:
            logger.warning(f"'{protein_name}' looks like UniProt ID but not found in database")
    
    # Case 3: Search by gene name (most reliable)
    gene_name = target_info.get("gene_name", "").strip() if target_info.get("gene_name") else ""
    if gene_name:
        logger.info(f"Dynamically searching UniProt API for gene: {gene_name}")
        try:
            results = uniprot_client.search_by_gene_name(gene_name, organism="human")
            if results:
                # Take the first reviewed entry if available
                for result in results:
                    parsed = uniprot_client.parse_protein_info(result)
                    if parsed.get("uniprot_id"):
                        logger.info(f"API returned UniProt ID for {gene_name}: {parsed['uniprot_id']}")
                        return parsed["uniprot_id"]
        except Exception as e:
            logger.error(f"Error searching by gene name {gene_name}: {e}")
    
    # Case 4: Search by protein name
    if protein_name and not is_uniprot_id_pattern(protein_name):
        logger.info(f"Dynamically searching UniProt API for protein: {protein_name}")
        try:
            results = uniprot_client.search_protein(protein_name, organism="human", limit=1)
            if results:
                parsed = uniprot_client.parse_protein_info(results[0])
                if parsed.get("uniprot_id"):
                    logger.info(f"API returned UniProt ID for {protein_name}: {parsed['uniprot_id']}")
                    return parsed["uniprot_id"]
        except Exception as e:
            logger.error(f"Error searching by protein name {protein_name}: {e}")
    
    logger.warning(f"Could not dynamically resolve UniProt ID for: {gene_name or protein_name}")
    return None


def is_uniprot_id_pattern(text: str) -> bool:
    """
    Check if text matches UniProt ID patterns.
    
    Patterns:
    - 6 characters: P01116 (Swiss-Prot)
    - 10 characters: A0A1L1T3F0 (TrEMBL)
    - Contains underscore: EGFR_HUMAN (entry name)
    
    Args:
        text: Text to check
        
    Returns:
        True if matches UniProt ID pattern
    """
    if not text:
        return False
    
    text = text.strip()
    
    # Standard patterns
    if len(text) == 6:  # Swiss-Prot format
        return text[0].isalpha() and text[1:6].replace('_', '').isalnum()
    
    if len(text) == 10:  # TrEMBL format
        return text[0:3].isalnum() and text[3:].isalnum()
    
    if '_' in text and len(text) <= 20:  # Entry name format
        parts = text.split('_')
        if len(parts) == 2:
            return parts[0].isalnum() and parts[1].isalpha()
    
    return False


def validate_uniprot_id(uniprot_id: str) -> bool:
    """
    Validate if a UniProt ID exists in the database via API call.
    
    Args:
        uniprot_id: UniProt ID to validate
        
    Returns:
        True if valid and exists
    """
    uniprot_client = get_uniprot_client()
    try:
        # Try to fetch the protein - this is an API call
        logger.debug(f"Validating {uniprot_id} via UniProt API...")
        protein = uniprot_client.get_protein_by_id(uniprot_id)
        return protein is not None
    except Exception as e:
        logger.debug(f"Error validating UniProt ID {uniprot_id}: {e}")
        return False
