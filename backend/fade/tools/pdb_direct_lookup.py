"""
PDB Direct Lookup Fix

This module adds proper handling for direct PDB ID queries.
When a user asks for a specific PDB ID, they should get that structure
regardless of whether it has drug-like ligands.
"""

from typing import Dict, Any, Optional, List, Tuple
from fade.utils import get_logger

logger = get_logger("tools.pdb_direct")


def fetch_pdb_directly(pdb_id: str, rcsb_client) -> Tuple[List[Dict], Dict]:
    """
    Fetch a specific PDB structure directly WITHOUT filtering.
    
    This is used when a user provides a PDB ID directly - they want
    THIS specific structure, not a filtered search result.
    
    Args:
        pdb_id: The PDB ID to fetch (e.g., "4G5J")
        rcsb_client: The RCSB client instance
        
    Returns:
        Tuple of (structures list, target_info dict)
    """
    logger.info(f"Direct PDB fetch (no filtering): {pdb_id}")
    
    # Use the search with the PDB ID but DON'T apply drug-like filters
    # We want the structure AS-IS
    structures, target_info = rcsb_client.search_by_query(
        pdb_id.upper(),
        limit=1,
        original_target=None,
        apply_drug_filter=False  # KEY: Don't filter for drug-like ligands
    )
    
    if structures:
        # Mark the structure as a direct fetch
        structures[0]["direct_pdb_fetch"] = True
        structures[0]["skip_drug_filter"] = True
        
        # Log what we found
        struct = structures[0]
        logger.info(f"Fetched PDB {pdb_id}: {struct.get('title', 'Unknown')[:100]}")
        
        if struct.get("ligands"):
            ligand_names = []
            for lig in struct["ligands"][:5]:
                name = f"{lig.get('id', 'Unknown')}"
                if lig.get('molecular_weight'):
                    name += f" (MW: {lig['molecular_weight']})"
                ligand_names.append(name)
            logger.info(f"  Contains ligands: {', '.join(ligand_names)}")
            
            # Note if they're not drug-like but DON'T filter them out
            non_drug_like = [l.get('id') for l in struct["ligands"] 
                           if l.get('id') in ['GDP', 'GTP', 'ATP', 'ADP', 'GNP', 'MG', 'CA', 'ZN']]
            if non_drug_like:
                logger.info(f"  Including non-drug-like ligands: {', '.join(non_drug_like)}")
                struct["has_non_drug_ligands"] = True
        else:
            logger.info(f"  No ligands in structure")
            
    return structures, target_info


def is_direct_pdb_query(target_info: Dict[str, Any]) -> bool:
    """
    Check if this is a direct PDB ID query (user wants a specific structure).
    
    Args:
        target_info: Parsed target information
        
    Returns:
        True if this is a direct PDB ID query
    """
    # It's a direct PDB query if:
    # 1. We have a PDB ID
    # 2. We don't have other identifiers (it's not a discovery query)
    has_pdb = bool(target_info.get("pdb_id"))
    has_gene = bool(target_info.get("gene_name"))
    has_protein = bool(target_info.get("protein_name"))
    has_uniprot = bool(target_info.get("uniprot_id"))
    
    # If ONLY PDB ID is provided, it's a direct lookup
    if has_pdb and not (has_gene or has_protein or has_uniprot):
        return True
    
    # If PDB ID is provided with other info, it might still be a direct lookup
    # Check if the PDB ID seems to be the main focus
    if has_pdb and target_info.get("is_pdb_focused"):
        return True
        
    return False
