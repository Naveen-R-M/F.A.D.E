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
    
    # If we don't have a sequence but have a UniProt ID, fetch it
    if not sequence and uniprot_id:
        logger.info(f"No sequence available, fetching from UniProt for {uniprot_id}")
        from fade.tools.uniprot_api import get_uniprot_client
        uniprot_client = get_uniprot_client()
        sequence = uniprot_client.get_protein_sequence(uniprot_id)
        if sequence:
            target_info["sequence"] = sequence
            target_info["sequence_length"] = len(sequence)
            logger.info(f"Fetched sequence: {len(sequence)} amino acids")
    
    if not sequence:
        # Check if protein_name might be a UniProt ID (like A0A1L1T3F0)
        protein_name = target_info.get("protein_name")
        if protein_name and (len(protein_name) == 6 or len(protein_name) == 10 or "_" in protein_name):
            # Looks like a UniProt ID pattern
            logger.info(f"Protein name '{protein_name}' looks like a UniProt ID, trying to fetch sequence")
            from fade.tools.uniprot_api import get_uniprot_client
            uniprot_client = get_uniprot_client()
            sequence = uniprot_client.get_protein_sequence(protein_name)
            if sequence:
                target_info["uniprot_id"] = protein_name
                target_info["sequence"] = sequence
                target_info["sequence_length"] = len(sequence)
                logger.info(f"Successfully fetched sequence for {protein_name}: {len(sequence)} amino acids")
    
    if not sequence:
        # Try to dynamically resolve UniProt ID and fetch sequence
        gene_name = target_info.get("gene_name")
        if gene_name:
            logger.info(f"Dynamically resolving UniProt ID for gene: {gene_name}")
            from fade.agents.research.uniprot_resolver import resolve_uniprot_id_dynamically
            from fade.tools.uniprot_api import get_uniprot_client
            
            # Use dynamic resolution - NO HARDCODED MAPPINGS
            resolved_uniprot_id = resolve_uniprot_id_dynamically(target_info)
            
            if resolved_uniprot_id:
                uniprot_client = get_uniprot_client()
                sequence = uniprot_client.get_protein_sequence(resolved_uniprot_id)
                if sequence:
                    target_info["uniprot_id"] = resolved_uniprot_id
                    target_info["sequence"] = sequence
                    target_info["sequence_length"] = len(sequence)
                    logger.info(f"Dynamically resolved {gene_name} → {resolved_uniprot_id}: {len(sequence)} amino acids")
    
    if not sequence:
        # Generate intelligent guidance for the user
        logger.info("Generating guidance for UniProt lookup failure")
        
        from fade.agents.structure.uniprot_guidance import (
            generate_uniprot_guidance,
            format_uniprot_guidance_message
        )
        
        # Determine what failed
        protein_name = target_info.get("protein_name", "")
        gene_name = target_info.get("gene_name", "")
        uniprot_id = target_info.get("uniprot_id", "")
        
        # The ID that failed to resolve
        failed_id = uniprot_id or protein_name or gene_name or "unknown"
        
        # Generate guidance
        guidance = generate_uniprot_guidance(
            query=state.get("query", "Unknown query"),
            uniprot_id=failed_id,
            error_type="not_found" if not uniprot_id else "no_sequence"
        )
        
        # Format the message
        guidance_message = format_uniprot_guidance_message(guidance)
        
        # Return with guidance instead of hard failure
        return {
            "error": guidance_message,
            "error_type": "uniprot_not_found",
            "guidance": guidance,
            "suggested_queries": guidance.get("suggestions", []),
            "should_continue": False,
            "current_step": "structure_failed_with_guidance"
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
    
    # If no PDB structures or couldn't download, continue to AlphaFold
    logger.info("No PDB structures available, checking AlphaFold database")
    
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
