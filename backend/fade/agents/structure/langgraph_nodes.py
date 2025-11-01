"""
Structure module LangGraph nodes.

These nodes handle structure resolution and preparation:
1. Try to find existing structure (PDB)
2. Try AlphaFold database
3. Fall back to Boltz-2 prediction
"""

import time
from typing import Dict, Any, Optional
from pathlib import Path

from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.tools import get_rcsb_client, get_boltz2_client
from fade.tools.alphafold_api import get_alphafold_client
from fade.config import config
from fade.utils import get_logger, select_best_structure

logger = get_logger("nodes.structure")


def structure_resolver_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for resolving protein structure.
    
    Priority order:
    1. Check for existing PDB structures (from target research)
    2. Try AlphaFold database for pre-computed structure
    3. Use Boltz-2 for de novo prediction
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with structure information
    """
    target_info = state.get("target_info", {})
    
    if not target_info:
        return {
            "error": "No target information available",
            "should_continue": False,
            "current_step": "structure_failed"
        }
    
    uniprot_id = target_info.get("uniprot_id")
    sequence = target_info.get("sequence")
    
    if not sequence:
        return {
            "error": "No protein sequence available",
            "should_continue": False,
            "current_step": "structure_failed"
        }
    
    logger.info(f"[Structure Resolver] Finding structure for {uniprot_id}")
    
    # Option 1: Check existing PDB structures
    existing_structures = target_info.get("existing_structures", [])
    if existing_structures:
        # Pass target_info for enhanced selection with mutations and known compounds
        best_structure = _select_best_pdb(existing_structures, target_info)
        if best_structure:
            pdb_content = _download_pdb(best_structure["pdb_id"])
            if pdb_content:
                logger.info(f"Using PDB structure: {best_structure['pdb_id']} (resolution: {best_structure.get('resolution')} Å)")
                
                structure_info = {
                    "structure_path": _save_structure(pdb_content, f"{best_structure['pdb_id']}.pdb"),
                    "pdb_id": best_structure["pdb_id"],
                    "resolution": best_structure.get("resolution"),
                    "source": "PDB",
                    "confidence_score": None,
                    "has_ligand": best_structure.get("has_ligand", False),
                    "ligand_names": best_structure.get("ligands", [])
                }
                
                message = f"Found experimental structure: PDB {best_structure['pdb_id']} ({best_structure.get('resolution', 'N/A')} Å resolution)"
                
                return {
                    "structure": structure_info,
                    "structure_source": "PDB",
                    "messages": [HumanMessage(content=message)],
                    "current_step": "structure_found",
                    "should_continue": True
                }
        else:
            # NO FALLBACK - Must have drug-like ligands
            return {
                "error": "No PDB structures with drug-like ligands found. Try a more specific query like 'TARGET_NAME kinase domain with inhibitor' or 'TARGET_NAME with DRUG_NAME'",
                "should_continue": False,
                "current_step": "structure_failed"
            }
    
    # Option 2: Try AlphaFold database
    if uniprot_id:
        alphafold_client = get_alphafold_client()
        
        # Check if AlphaFold structure exists with good confidence
        if alphafold_client.check_confidence(uniprot_id, min_plddt=70.0):
            pdb_content = alphafold_client.download_pdb(uniprot_id)
            
            if pdb_content:
                logger.info(f"Using AlphaFold structure for {uniprot_id}")
                
                structure_info = {
                    "structure_path": _save_structure(pdb_content, f"AF-{uniprot_id}.pdb"),
                    "pdb_id": None,
                    "resolution": None,
                    "source": "AlphaFold",
                    "confidence_score": 70.0,  # We checked for min 70
                    "has_ligand": False,
                    "ligand_names": []
                }
                
                message = f"Found AlphaFold pre-computed structure for {uniprot_id} (pLDDT > 70)"
                
                return {
                    "structure": structure_info,
                    "structure_source": "AlphaFold",
                    "messages": [HumanMessage(content=message)],
                    "current_step": "structure_found",
                    "should_continue": True
                }
    
    # Option 3: Use Boltz-2 for prediction
    logger.info("No existing structure found, using Boltz-2 for prediction")
    
    boltz2_client = get_boltz2_client()
    
    try:
        # Submit Boltz-2 job
        job_result = boltz2_client.predict_structure(
            sequence=sequence,
            job_name=f"fade_{uniprot_id or 'unknown'}_{int(time.time())}"
        )
        
        job_id = job_result.get("job_id")
        
        if not job_id:
            return {
                "error": "Failed to submit Boltz-2 job",
                "should_continue": False,
                "current_step": "structure_failed"
            }
        
        logger.info(f"Boltz-2 job submitted: {job_id}")
        message = f"No existing structure found. Submitted Boltz-2 prediction job: {job_id}"
        
        # Store job ID for next node to wait for result
        structure_info = {
            "boltz2_job_id": job_id,
            "structure_path": None,
            "source": "Boltz2",
            "status": "pending"
        }
        
        return {
            "structure": structure_info,
            "structure_source": "Boltz2",
            "messages": [HumanMessage(content=message)],
            "current_step": "structure_predicting",
            "should_continue": True
        }
        
    except Exception as e:
        logger.error(f"Boltz-2 submission failed: {e}")
        return {
            "error": f"Structure prediction failed: {str(e)}",
            "should_continue": False,
            "current_step": "structure_failed"
        }


def structure_wait_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node to wait for Boltz-2 prediction to complete.
    
    This is a separate node to allow for async/parallel processing.
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with completed structure
    """
    structure_info = state.get("structure", {})
    
    # Only wait if we have a pending Boltz-2 job
    if structure_info.get("source") != "Boltz2" or structure_info.get("status") != "pending":
        return {"current_step": "structure_ready"}
    
    job_id = structure_info.get("boltz2_job_id")
    if not job_id:
        return {
            "error": "No Boltz-2 job ID found",
            "should_continue": False,
            "current_step": "structure_failed"
        }
    
    logger.info(f"[Structure Wait] Waiting for Boltz-2 job: {job_id}")
    
    boltz2_client = get_boltz2_client()
    
    # Wait for result (with timeout)
    pdb_content = boltz2_client.wait_for_structure(
        job_id=job_id,
        max_wait=1800,  # 30 minutes
        poll_interval=30  # Check every 30 seconds
    )
    
    if not pdb_content:
        return {
            "error": "Boltz-2 prediction timed out or failed",
            "should_continue": False,
            "current_step": "structure_failed"
        }
    
    # Save structure
    target_name = state.get("target_info", {}).get("gene_name", "unknown")
    structure_path = _save_structure(pdb_content, f"boltz2_{target_name}_{job_id}.pdb")
    
    # Update structure info
    structure_info.update({
        "structure_path": structure_path,
        "status": "completed",
        "confidence_score": 80.0,  # Boltz-2 typical confidence
        "has_ligand": False,
        "ligand_names": []
    })
    
    message = f"✓ Boltz-2 structure prediction completed for job {job_id}"
    
    return {
        "structure": structure_info,
        "messages": [HumanMessage(content=message)],
        "current_step": "structure_ready",
        "should_continue": True
    }


def _select_best_pdb(structures: list, target_info: Optional[Dict] = None) -> Optional[Dict[str, Any]]:
    """
    Select the best PDB structure using enhanced ligand-aware selection.
    
    Args:
        structures: List of PDB structures with ligand information
        target_info: Optional target information including mutations and known compounds
        
    Returns:
        Best structure or None
    """
    if not structures:
        return None
    
    # Extract mutations and known compounds from target info if available
    mutations = target_info.get("mutations", []) if target_info else []
    known_compounds = target_info.get("known_compounds", []) if target_info else []
    
    # Use the enhanced selection function from utils
    # NO FALLBACK - Must have drug-like ligands
    best = select_best_structure(
        structures,
        target_mutations=mutations,
        known_compounds=known_compounds,
        fallback_to_any=False  # NO FALLBACK - strict requirement
    )
    
    if best:
        logger.info(f"Selected PDB structure: {best.get('pdb_id')} (score: {best.get('quality_score', 0):.1f})")
        if best.get('drug_like_ligands'):
            ligand_names = [l.get('name', l.get('id', 'Unknown')) for l in best['drug_like_ligands']]
            logger.info(f"  Contains drug-like ligands: {', '.join(ligand_names)}")
    
    return best


def _download_pdb(pdb_id: str) -> Optional[str]:
    """Download PDB structure."""
    rcsb_client = get_rcsb_client()
    return rcsb_client.download_structure(pdb_id)


def _save_structure(pdb_content: str, filename: str) -> str:
    """Save PDB content to file."""
    structures_dir = config.DATA_DIR / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    
    filepath = structures_dir / filename
    with open(filepath, 'w') as f:
        f.write(pdb_content)
    
    logger.info(f"Saved structure to {filepath}")
    return str(filepath)
