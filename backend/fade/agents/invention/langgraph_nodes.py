"""
Invention module LangGraph nodes.

These nodes handle molecule generation and filtering:
1. molecule_generator_node: Generates molecules using DiffSBDD
2. medicinal_filter_node: Filters molecules by drug-likeness
"""

from typing import Dict, Any, List, Optional
from pathlib import Path

from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.tools.hpc_diffsbdd import get_hpc_diffsbdd_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.invention")


def molecule_generator_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for generating molecules using DiffSBDD.
    
    Uses the selected pocket to generate novel drug-like molecules.
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with generated molecules
    """
    structure = state.get("structure", {})
    selected_pocket = state.get("selected_pocket", {})
    
    if not structure or not structure.get("structure_path"):
        return {
            "error": "No protein structure available for molecule generation",
            "should_continue": False,
            "current_step": "molecule_generation_failed"
        }
    
    if not selected_pocket:
        return {
            "error": "No pocket selected for molecule generation",
            "should_continue": False,
            "current_step": "molecule_generation_failed"
        }
    
    structure_path = structure["structure_path"]
    
    logger.info(f"[Molecule Generator] Generating molecules for pocket {selected_pocket['pocket_id']}")
    
    # Get HPC DiffSBDD client
    diffsbdd_client = get_hpc_diffsbdd_client()
    
    # Format pocket residues for DiffSBDD
    pocket_residues = _format_pocket_residues(selected_pocket)
    from Bio.PDB import PDBParser
    protein = PDBParser(QUIET=True).get_structure("target", structure_path)[0]
    valid = {f"{chain.id}:{res.id[1]}" for chain in protein for res in chain
             if res.id[0] == " " and res.id[2] == " "}
    if not pocket_residues or not set(pocket_residues).issubset(valid):
        return {"error": "Pocket residues are missing or do not match the protein PDB",
                "should_continue": False, "current_step": "molecule_generation_failed"}
    
    # Use existing job_id if available (from fpocket), otherwise generate new one
    job_id = state.get("job_id")
    if job_id:
        logger.info(f"Using existing job_id from fpocket: {job_id}")
    
    # Determine number of samples based on config
    n_samples = min(
        config.MAX_MOLECULES_TO_GENERATE,
        100  # DiffSBDD batch size limit
    )
    
    logger.info(f"Requesting generation of {n_samples} molecules")
    
    # Generate molecules (returns molecules and job_id)
    generated_molecules, used_job_id = diffsbdd_client.generate_molecules(
        pdb_file=structure_path,
        pocket_residues=pocket_residues,
        n_samples=n_samples,
        sanitize=True,
        job_id=job_id  # Pass existing job_id if available
    )
    
    # Update job_id if it was generated
    if not job_id and used_job_id:
        job_id = used_job_id
        logger.info(f"Generated new job_id: {job_id}")
    
    if not generated_molecules:
        return {
            "error": "No molecules generated",
            "should_continue": False,
            "current_step": "molecule_generation_failed"
        }
    
    # Add basic properties to molecules
    for mol in generated_molecules:
        # Add placeholder properties (will be calculated in filter node)
        mol["needs_property_calculation"] = True
        mol["passed_generation"] = True
    
    message = f"Generated {len(generated_molecules)} molecules using DiffSBDD"
    if generated_molecules and generated_molecules[0].get("smiles"):
        message += f"\nFirst molecule: {generated_molecules[0]['smiles'][:50]}..."
    
    logger.info(f"Successfully generated {len(generated_molecules)} molecules")
    
    # Store generation parameters for reference
    generation_params = {
        "method": "DiffSBDD",
        "pocket_id": selected_pocket["pocket_id"],
        "pocket_residues": pocket_residues,
        "n_requested": n_samples,
        "n_generated": len(generated_molecules),
        "sanitized": True,
        "job_id": job_id
    }
    
    updates = {
        "generated_molecules": generated_molecules,
        "generation_parameters": generation_params,
        "messages": [HumanMessage(content=message)],
        "current_step": "molecules_generated",
        "should_continue": True
    }
    
    # Update job_id if needed
    if job_id and not state.get("job_id"):
        updates["job_id"] = job_id
    
    return updates


def medicinal_filter_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for filtering molecules by drug-likeness.
    
    Applies Lipinski's Rule of Five and other medicinal chemistry filters.
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with filtered molecules
    """
    generated_molecules = state.get("generated_molecules", [])
    
    if not generated_molecules:
        return {
            "error": "No molecules available for filtering",
            "should_continue": False,
            "current_step": "filtering_failed"
        }
    
    logger.info(f"[Medicinal Filter] Filtering {len(generated_molecules)} molecules")
    
    filtered_molecules = []
    filter_stats = {
        "total": len(generated_molecules),
        "no_smiles": 0,
        "invalid_smiles": 0,
        "lipinski_violations": 0,
        "property_violations": 0,
        "passed": 0
    }
    
    for mol in generated_molecules:
        # Check if molecule has SMILES
        if not mol.get("smiles"):
            filter_stats["no_smiles"] += 1
            continue
        
        # For now, use simple filters (without RDKit)
        # In production, you'd calculate real properties with RDKit
        filter_result = _apply_simple_filters(mol)
        
        if filter_result["passed"]:
            filtered_mol = mol.copy()
            filtered_mol.update(filter_result["properties"])
            filtered_mol["passed_filters"] = True
            filtered_mol["filter_warnings"] = filter_result.get("warnings", [])
            filtered_molecules.append(filtered_mol)
            filter_stats["passed"] += 1
        else:
            if filter_result.get("reason") == "invalid_smiles":
                filter_stats["invalid_smiles"] += 1
            elif "lipinski" in filter_result.get("reason", ""):
                filter_stats["lipinski_violations"] += 1
            else:
                filter_stats["property_violations"] += 1
    
    # Sort by synthetic accessibility or other metrics
    filtered_molecules.sort(key=lambda x: x.get("qed_score", 0.5), reverse=True)
    
    # Limit to max molecules for screening
    max_for_screening = config.MAX_MOLECULES_TO_SCREEN
    if len(filtered_molecules) > max_for_screening:
        logger.info(f"Limiting to top {max_for_screening} molecules for screening")
        filtered_molecules = filtered_molecules[:max_for_screening]
    
    message = f"Filtered {len(generated_molecules)} molecules"
    message += f"\n- Passed filters: {filter_stats['passed']}"
    message += f"\n- No SMILES: {filter_stats['no_smiles']}"
    message += f"\n- Lipinski violations: {filter_stats['lipinski_violations']}"
    message += f"\n- Other violations: {filter_stats['property_violations']}"
    message += f"\n- Selected for screening: {len(filtered_molecules)}"
    
    logger.info(f"Filtered to {len(filtered_molecules)} drug-like molecules")
    
    return {
        "filtered_molecules": filtered_molecules,
        "filter_statistics": filter_stats,
        "messages": [HumanMessage(content=message)],
        "current_step": "molecules_filtered",
        "should_continue": len(filtered_molecules) > 0
    }


def _format_pocket_residues(pocket: Dict[str, Any]) -> List[str]:
    """
    Format pocket residues for DiffSBDD input.
    
    Preserve explicit chain and residue numbering; never substitute a generic site.
    
    Args:
        pocket: Pocket information
        
    Returns:
        Unique residue identifiers in CHAIN:RESNUM format
    """
    import re
    residues = pocket.get("residues") or []
    if not residues or any(not re.fullmatch(r"[^\s:]:-?\d+", r) for r in residues):
        raise ValueError("Pocket must contain explicit CHAIN:RESNUM residue identities")
    return list(dict.fromkeys(residues))


def _apply_simple_filters(mol: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply basic RDKit physicochemical filters (not a safety or efficacy assessment).
    
    Args:
        mol: Molecule dictionary
        
    Returns:
        Filter result with pass/fail and properties
    """
    smiles = mol.get("smiles", "")
    
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, QED
    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None:
        return {"passed": False, "reason": "invalid_smiles", "properties": {}}
    mw_estimate = Descriptors.MolWt(parsed)
    hba_estimate = Lipinski.NumHAcceptors(parsed)
    hbd_estimate = Lipinski.NumHDonors(parsed)
    logp_estimate = Crippen.MolLogP(parsed)
    
    # Apply Lipinski's Rule of Five
    violations = 0
    warnings = []
    
    if mw_estimate > 500:
        violations += 1
        warnings.append(f"MW estimate ({mw_estimate}) > 500")
    
    if logp_estimate > 5:
        violations += 1
        warnings.append(f"LogP estimate ({logp_estimate:.1f}) > 5")
    
    if hba_estimate > 10:
        violations += 1
        warnings.append(f"HBA estimate ({hba_estimate}) > 10")
    
    if hbd_estimate > 5:
        violations += 1
        warnings.append(f"HBD estimate ({hbd_estimate}) > 5")
    
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    has_pains = FilterCatalog(params).HasMatch(parsed)
    if has_pains:
        warnings.append("Contains a PAINS alert")
    qed_score = QED.qed(parsed)
    passed = violations <= 1 and not has_pains and mw_estimate >= 150
    
    return {
        "passed": passed,
        "reason": "lipinski" if violations > 1 else ("pains" if has_pains else ("MW<150" if mw_estimate < 150 else "")),
        "warnings": warnings,
        "properties": {
            "mol_weight": mw_estimate,
            "logp": logp_estimate,
            "hba": hba_estimate,
            "hbd": hbd_estimate,
            "lipinski_violations": violations,
            "qed_score": qed_score,
            "sa_score": None,  # Not calculated; do not fabricate a score.
            "has_pains_alert": has_pains
        }
    }
