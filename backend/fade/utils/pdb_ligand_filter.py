"""
Utility functions for filtering and ranking PDB structures based on ligand content.

This module provides functions to identify drug-like ligands in PDB structures
and rank structures based on their suitability for drug discovery.
"""

import logging
from typing import Dict, Any, List, Optional
from datetime import datetime

from fade.config import config
from fade.utils import get_logger

logger = get_logger("utils.pdb_ligand_filter")


def is_drug_like_ligand(ligand: Dict[str, Any]) -> bool:
    """
    Check if a ligand meets drug-like criteria.
    
    Filters out crystallization artifacts, ions, nucleotides, and non-drug-like molecules.
    
    Args:
        ligand: Dictionary containing ligand information with keys:
            - id: Ligand identifier (3-letter code)
            - molecular_weight: MW in Daltons
            - heavy_atom_count: Number of non-hydrogen atoms
            - type: Ligand type/classification
            - formula: Chemical formula
    
    Returns:
        True if the ligand appears drug-like, False otherwise
    """
    ligand_id = ligand.get("id", "").upper()
    mw = ligand.get("molecular_weight", 0)
    
    # Handle case where molecular_weight might be a list or None
    if isinstance(mw, list):
        mw = mw[0] if mw else 0
    elif mw is None:
        mw = 0
    
    heavy_atoms = ligand.get("heavy_atom_count", 0)
    ligand_type = ligand.get("type", "").upper()
    
    # Quick reject for known artifacts
    if ligand_id in config.PDB_CRYSTALLIZATION_ARTIFACTS:
        logger.debug(f"Rejecting {ligand_id}: crystallization artifact")
        return False
    
    # Reject nucleotides
    if ligand_id in config.PDB_NUCLEOTIDES:
        logger.debug(f"Rejecting {ligand_id}: nucleotide")
        return False
    
    # Check molecular weight range
    if not (config.PDB_MIN_LIGAND_MW < mw < config.PDB_MAX_LIGAND_MW):
        logger.debug(f"Rejecting {ligand_id}: MW {mw} outside range {config.PDB_MIN_LIGAND_MW}-{config.PDB_MAX_LIGAND_MW}")
        return False
    
    # Check heavy atom count
    if heavy_atoms < config.PDB_MIN_HEAVY_ATOMS:
        logger.debug(f"Rejecting {ligand_id}: only {heavy_atoms} heavy atoms")
        return False
    
    # Reject peptides and other polymers (if type information available)
    if ligand_type in ["PEPTIDE", "POLYMER", "DNA", "RNA"]:
        logger.debug(f"Rejecting {ligand_id}: type is {ligand_type}")
        return False
    
    # Additional checks based on chemical formula if available
    formula = ligand.get("formula", "")
    if formula:
        # Simple heuristic: reject if only contains C, H, O in simple ratios (likely sugar)
        if _is_likely_sugar(formula):
            logger.debug(f"Rejecting {ligand_id}: likely sugar based on formula {formula}")
            return False
    
    logger.debug(f"Accepting {ligand_id} as drug-like: MW={mw}, heavy_atoms={heavy_atoms}")
    return True


def _is_likely_sugar(formula: str) -> bool:
    """
    Simple heuristic to identify likely sugar molecules.
    
    Args:
        formula: Chemical formula string
        
    Returns:
        True if formula suggests a sugar molecule
    """
    # Common sugar formulas
    sugar_formulas = {
        "C6H12O6",  # Glucose
        "C5H10O5",  # Ribose
        "C12H22O11",  # Sucrose
        "C6H10O5",  # Polysaccharide unit
    }
    
    formula_clean = formula.replace(" ", "").upper()
    return formula_clean in sugar_formulas


def filter_structures_by_ligands(structures: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Filter PDB structures to only those containing drug-like ligands.
    
    Args:
        structures: List of PDB structure dictionaries
        
    Returns:
        Filtered list containing only structures with drug-like ligands
    """
    filtered = []
    
    for structure in structures:
        ligands = structure.get("ligands", [])
        drug_like_ligands = [lig for lig in ligands if is_drug_like_ligand(lig)]
        
        if drug_like_ligands:
            # Update structure with filtered ligands
            structure["drug_like_ligands"] = drug_like_ligands
            structure["has_drug_like_ligand"] = True
            filtered.append(structure)
            logger.info(f"PDB {structure.get('pdb_id')} has {len(drug_like_ligands)} drug-like ligands")
    
    logger.info(f"Filtered {len(structures)} structures to {len(filtered)} with drug-like ligands")
    return filtered


def score_structure_with_ligands(
    structure: Dict[str, Any],
    target_mutations: Optional[List[str]] = None,
    known_compounds: Optional[List[Dict]] = None
) -> float:
    """
    Score a PDB structure based on multiple criteria including ligand quality.
    
    Args:
        structure: PDB structure dictionary with ligand information
        target_mutations: Optional list of mutations (e.g., ["G12C"])
        known_compounds: Optional list of known compounds from ChEMBL
        
    Returns:
        Numerical score (higher is better)
    """
    score = 0.0
    pdb_id = structure.get("pdb_id", "unknown")
    
    # 1. Drug-like ligand presence (most important)
    drug_like_ligands = structure.get("drug_like_ligands", [])
    if drug_like_ligands:
        score += 100
        logger.debug(f"{pdb_id}: +100 for having drug-like ligands")
        
        # Bonus for multiple drug-like ligands (SAR information)
        if len(drug_like_ligands) > 1:
            score += 10 * min(len(drug_like_ligands) - 1, 5)  # Cap at 5 extra ligands
            logger.debug(f"{pdb_id}: +{10 * min(len(drug_like_ligands) - 1, 5)} for multiple ligands")
    
    # 2. Mutation-specific binding (if applicable)
    if target_mutations and structure.get("ligand_at_mutation_site"):
        score += 50
        logger.debug(f"{pdb_id}: +50 for ligand at mutation site")
    
    # 3. Known inhibitor present
    if known_compounds and _has_known_compound(structure, known_compounds):
        score += 25
        logger.debug(f"{pdb_id}: +25 for having known inhibitor")
    
    # 4. Resolution (lower is better, max 10 points)
    resolution = structure.get("resolution")
    
    # Handle case where resolution might be a list
    if isinstance(resolution, list):
        resolution = resolution[0] if resolution else None
    
    if resolution and resolution > 0:
        try:
            resolution = float(resolution)
            resolution_score = min(10, 10 / resolution)
            score += resolution_score
            logger.debug(f"{pdb_id}: +{resolution_score:.1f} for resolution {resolution} Å")
        except (TypeError, ValueError):
            pass
    
    # 5. Recency (newer structures often have better quality)
    release_date = structure.get("release_date")
    if release_date:
        try:
            year = int(release_date[:4])
            current_year = datetime.now().year
            if current_year - year <= config.PDB_PREFER_RECENT_YEARS:
                recency_score = (config.PDB_PREFER_RECENT_YEARS - (current_year - year)) * 2
                score += recency_score
                logger.debug(f"{pdb_id}: +{recency_score} for recent structure ({year})")
        except (ValueError, TypeError):
            pass
    
    # 6. Binding affinity data available
    if structure.get("has_binding_affinity"):
        score += 10
        logger.debug(f"{pdb_id}: +10 for having binding affinity data")
    
    # 7. Experimental method bonus
    method = structure.get("experimental_method", "").upper()
    if "X-RAY" in method:
        score += 5
        logger.debug(f"{pdb_id}: +5 for X-ray crystallography")
    elif "CRYO-EM" in method:
        score += 3
        logger.debug(f"{pdb_id}: +3 for Cryo-EM")
    
    logger.info(f"Structure {pdb_id} scored {score:.1f}")
    return score


def _has_known_compound(structure: Dict[str, Any], known_compounds: List[Dict]) -> bool:
    """
    Check if structure contains any known compounds from ChEMBL.
    
    Args:
        structure: PDB structure dictionary
        known_compounds: List of known compounds from ChEMBL
        
    Returns:
        True if any ligand matches a known compound
    """
    if not known_compounds:
        return False
        
    structure_ligands = structure.get("drug_like_ligands", [])
    if not structure_ligands:
        return False
    
    for ligand in structure_ligands:
        ligand_name = ligand.get("name", "").upper()
        ligand_id = ligand.get("id", "").upper()
        
        for compound in known_compounds:
            compound_name = compound.get("molecule_name", "")
            if compound_name:
                compound_name = compound_name.upper()
                # Check for name match (accounting for variations)
                if (
                    compound_name in ligand_name or 
                    ligand_name in compound_name or
                    ligand_id == compound.get("molecule_chembl_id", "")[:3]  # Sometimes ChEMBL ID prefix matches
                ):
                    logger.debug(f"Found known compound match: {ligand_name} ~ {compound_name}")
                    return True
    
    return False


def rank_structures_by_quality(
    structures: List[Dict[str, Any]],
    target_mutations: Optional[List[str]] = None,
    known_compounds: Optional[List[Dict]] = None,
    require_drug_like: bool = True
) -> List[Dict[str, Any]]:
    """
    Rank PDB structures by quality for drug discovery.
    
    Args:
        structures: List of PDB structure dictionaries
        target_mutations: Optional list of mutations to prioritize
        known_compounds: Optional list of known compounds from ChEMBL
        require_drug_like: If True, only consider structures with drug-like ligands
        
    Returns:
        Sorted list of structures (best first) with scores
    """
    # Filter for drug-like ligands if required
    if require_drug_like:
        structures = filter_structures_by_ligands(structures)
        if not structures:
            logger.warning("No structures with drug-like ligands found")
            return []
    
    # Score each structure
    scored_structures = []
    for structure in structures:
        score = score_structure_with_ligands(structure, target_mutations, known_compounds)
        structure["quality_score"] = score
        scored_structures.append(structure)
    
    # Sort by score (highest first)
    scored_structures.sort(key=lambda x: x["quality_score"], reverse=True)
    
    # Log ranking results
    logger.info("Structure ranking results:")
    for i, structure in enumerate(scored_structures[:5], 1):
        logger.info(f"  {i}. {structure.get('pdb_id')} - Score: {structure['quality_score']:.1f}")
    
    return scored_structures


def select_best_structure(
    structures: List[Dict[str, Any]],
    target_mutations: Optional[List[str]] = None,
    known_compounds: Optional[List[Dict]] = None,
    fallback_to_any: bool = True
) -> Optional[Dict[str, Any]]:
    """
    Select the best PDB structure for drug discovery.
    
    Args:
        structures: List of PDB structure dictionaries
        target_mutations: Optional list of mutations to prioritize
        known_compounds: Optional list of known compounds from ChEMBL
        fallback_to_any: If True and no drug-like structures found, return best overall
        
    Returns:
        Best structure or None if no suitable structure found
    """
    if not structures:
        return None
    
    # Try to get best structure with drug-like ligands
    ranked = rank_structures_by_quality(
        structures, 
        target_mutations, 
        known_compounds,
        require_drug_like=True
    )
    
    if ranked:
        best = ranked[0]
        logger.info(f"Selected {best.get('pdb_id')} with drug-like ligands (score: {best['quality_score']:.1f})")
        return best
    
    # Fallback to best structure without drug-like requirement
    if fallback_to_any:
        logger.warning("No structures with drug-like ligands, falling back to best available")
        ranked = rank_structures_by_quality(
            structures,
            target_mutations,
            known_compounds,
            require_drug_like=False
        )
        if ranked:
            best = ranked[0]
            logger.info(f"Selected {best.get('pdb_id')} without drug-like ligands (score: {best['quality_score']:.1f})")
            return best
    
    logger.error("No suitable structure found")
    return None
