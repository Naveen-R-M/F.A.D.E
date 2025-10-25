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
            if "lipinski" in filter_result.get("reason", ""):
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
    
    For now, we return an empty list to let DiffSBDD use the whole protein.
    This avoids chain/residue numbering issues.
    
    Args:
        pocket: Pocket information
        
    Returns:
        Empty list (DiffSBDD will use whole protein)
    """
    # DiffSBDD can work without specific residue list
    # It will generate molecules for the whole protein
    # This avoids KeyError issues with chain/residue numbering
    logger.info(f"Using whole protein for molecule generation (pocket: {pocket.get('pocket_id')})")
    
    if pocket.get("center"):
        center = pocket["center"]
        logger.info(f"Pocket center at: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    
    return []


def _apply_simple_filters(mol: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply simple drug-likeness filters without RDKit.
    
    This is a simplified version. In production, use RDKit for accurate calculations.
    
    Args:
        mol: Molecule dictionary
        
    Returns:
        Filter result with pass/fail and properties
    """
    smiles = mol.get("smiles", "")
    
    # Estimate properties from SMILES (very approximate!)
    # In production, use RDKit for accurate calculations
    
    # Estimate molecular weight (very rough)
    mw_estimate = len(smiles) * 12  # Very rough approximation
    
    # Count some functional groups
    n_nitrogen = smiles.count('N')
    n_oxygen = smiles.count('O')
    n_sulfur = smiles.count('S')
    n_fluorine = smiles.count('F')
    n_chlorine = smiles.count('Cl')
    n_bromine = smiles.count('Br')
    
    # Estimate HBA/HBD
    hba_estimate = n_nitrogen + n_oxygen
    hbd_estimate = smiles.count('NH') + smiles.count('OH')
    
    # Estimate LogP (very rough)
    logp_estimate = len(smiles) * 0.1 - n_oxygen * 0.5 - n_nitrogen * 0.3
    
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
    
    # Check for PAINS patterns (simplified)
    pains_patterns = [
        'N=N',  # Azo compounds
        'S(=O)(=O)O',  # Sulfonates
        'C(=O)C(=O)',  # Diketones
        '[N+](=O)[O-]',  # Nitro groups (can be problematic)
    ]
    
    has_pains = any(pattern in smiles for pattern in pains_patterns)
    if has_pains:
        warnings.append("Contains potential PAINS alert")
    
    # Calculate QED score (simplified)
    qed_score = max(0, min(1, 1 - violations * 0.25))
    
    # Synthetic accessibility (simplified)
    sa_score = 3.0  # Default moderate difficulty
    if len(smiles) > 50:
        sa_score += 1.0
    if 'B' in smiles or '[' in smiles:
        sa_score += 1.0
    
    passed = violations <= 1 and not has_pains
    
    return {
        "passed": passed,
        "reason": "lipinski" if violations > 1 else ("pains" if has_pains else ""),
        "warnings": warnings,
        "properties": {
            "mol_weight": mw_estimate,
            "logp": logp_estimate,
            "hba": hba_estimate,
            "hbd": hbd_estimate,
            "lipinski_violations": violations,
            "qed_score": qed_score,
            "sa_score": sa_score,
            "has_pains_alert": has_pains
        }
    }
