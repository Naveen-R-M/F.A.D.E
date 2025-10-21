"""
Targeting module LangGraph nodes.

These nodes handle pocket detection and ranking:
1. pocket_mapper_node: Finds all pockets using fpocket
2. pocket_ranker_node: Ranks and selects best pocket
"""

from typing import Dict, Any, List, Optional
from pathlib import Path

from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.tools.fpocket import get_fpocket_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.targeting")


def pocket_mapper_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for detecting binding pockets.
    
    Uses fpocket to find potential binding sites in the protein structure.
    This also sets the job_id that will be used throughout the pipeline.
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with detected pockets and job_id
    """
    structure = state.get("structure", {})
    
    if not structure or not structure.get("structure_path"):
        return {
            "error": "No protein structure available for pocket detection",
            "should_continue": False,
            "current_step": "pocket_detection_failed"
        }
    
    structure_path = structure["structure_path"]
    
    if not Path(structure_path).exists():
        return {
            "error": f"Structure file not found: {structure_path}",
            "should_continue": False,
            "current_step": "pocket_detection_failed"
        }
    
    logger.info(f"[Pocket Mapper] Detecting pockets in {structure_path}")
    
    # Use fpocket for detection
    fpocket_client = get_fpocket_client()
    
    # Detect pockets with minimum druggability score
    pockets = fpocket_client.detect_pockets(
        pdb_file=structure_path,
        min_score=0.5  # Minimum druggability score
    )
    
    # Extract job_id if it was set by HPC fpocket
    job_id = None
    if pockets and pockets[0].get("remote_path"):
        # Extract job_id from remote path
        import re
        match = re.search(r'/fpocket_inputs/([^/]+)/', pockets[0]["remote_path"])
        if match:
            job_id = match.group(1)
            logger.info(f"Setting pipeline job_id: {job_id}")
    
    if not pockets:
        return {
            "error": "No pockets detected in structure",
            "should_continue": False,
            "current_step": "pocket_detection_failed"
        }
    
    # Limit to top N pockets
    max_pockets = config.MAX_POCKETS_TO_ANALYZE
    if len(pockets) > max_pockets:
        logger.info(f"Found {len(pockets)} pockets, limiting to top {max_pockets}")
        pockets = pockets[:max_pockets]
    
    # Save pocket information
    pockets_dir = config.DATA_DIR / "pockets" / state["run_id"]
    pockets_dir.mkdir(parents=True, exist_ok=True)
    
    # Save each pocket as PDB for visualization
    for pocket in pockets:
        pocket_file = pockets_dir / f"{pocket['pocket_id']}.pdb"
        fpocket_client.save_pocket_pdb(pocket, str(pocket_file))
        pocket["pocket_file"] = str(pocket_file)
    
    message = f"Found {len(pockets)} potential binding pockets"
    if pockets:
        best_pocket = pockets[0]
        message += f"\nBest pocket: {best_pocket['pocket_id']} (druggability: {best_pocket.get('druggability_score', 0):.2f})"
    
    logger.info(message)
    
    updates = {
        "pockets": pockets,
        "messages": [HumanMessage(content=message)],
        "current_step": "pockets_detected",
        "should_continue": True
    }
    
    # Add job_id if it was extracted
    if job_id:
        updates["job_id"] = job_id
    
    return updates


def pocket_ranker_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node for ranking and selecting the best pocket.
    
    Considers:
    1. Druggability score
    2. Volume (optimal range: 300-1000 Å³)
    3. Known binding sites (if available)
    4. Proximity to mutations (if specified)
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with selected pocket
    """
    pockets = state.get("pockets", [])
    
    if not pockets:
        return {
            "error": "No pockets available for ranking",
            "should_continue": False,
            "current_step": "pocket_ranking_failed"
        }
    
    logger.info(f"[Pocket Ranker] Ranking {len(pockets)} pockets")
    
    target_info = state.get("target_info", {})
    mutations = target_info.get("mutations", [])
    
    # Score each pocket
    scored_pockets = []
    for pocket in pockets:
        score = _calculate_pocket_score(pocket, mutations)
        pocket_with_score = pocket.copy()
        pocket_with_score["final_score"] = score
        scored_pockets.append(pocket_with_score)
    
    # Sort by score
    scored_pockets.sort(key=lambda x: x["final_score"], reverse=True)
    
    # Select the best pocket
    selected_pocket = scored_pockets[0]
    
    # Generate rationale
    rationale = _generate_selection_rationale(selected_pocket, mutations)
    
    message = f"Selected pocket: {selected_pocket['pocket_id']}"
    message += f"\nDruggability score: {selected_pocket.get('druggability_score', 0):.2f}"
    message += f"\nVolume: {selected_pocket.get('volume', 0):.1f} Ų"
    
    if selected_pocket.get("description"):
        message += f"\nDescription: {selected_pocket['description']}"
    
    logger.info(f"Selected {selected_pocket['pocket_id']} as best pocket")
    
    return {
        "pockets": scored_pockets,  # Update with scores
        "selected_pocket": selected_pocket,
        "pocket_selection_rationale": rationale,
        "messages": [HumanMessage(content=message)],
        "current_step": "pocket_selected",
        "should_continue": True
    }


def _calculate_pocket_score(pocket: Dict[str, Any], mutations: List[str]) -> float:
    """
    Calculate a comprehensive score for pocket ranking.
    
    Args:
        pocket: Pocket information
        mutations: List of mutations (e.g., ["G12C"])
        
    Returns:
        Final score (0-100)
    """
    score = 0.0
    
    # Druggability score (40% weight)
    druggability = pocket.get("druggability_score", 0.5)
    score += druggability * 40
    
    # Volume score (20% weight) - optimal range 300-1000 Ų
    volume = pocket.get("volume", 500)
    if 300 <= volume <= 1000:
        volume_score = 1.0
    elif volume < 300:
        volume_score = volume / 300
    else:
        volume_score = max(0, 2 - volume / 1000)
    score += volume_score * 20
    
    # Known binding site bonus (20% weight)
    if pocket.get("is_known_site", False):
        score += 20
    
    # Mutation proximity bonus (20% weight)
    if mutations and pocket.get("residues"):
        mutation_residues = set()
        for mutation in mutations:
            # Extract residue number from mutation (e.g., G12C -> 12)
            import re
            match = re.match(r'[A-Z](\d+)[A-Z]', mutation)
            if match:
                residue_num = match.group(1)
                mutation_residues.add(residue_num)
        
        # Check if any pocket residues match mutation sites
        pocket_residues = set()
        for res in pocket.get("residues", []):
            # Extract number from residue (e.g., G12 -> 12)
            match = re.match(r'[A-Z]*(\d+)', res)
            if match:
                pocket_residues.add(match.group(1))
        
        if mutation_residues & pocket_residues:
            score += 20
            logger.info(f"Pocket {pocket['pocket_id']} contains mutation site(s)")
    
    return score


def _generate_selection_rationale(pocket: Dict[str, Any], mutations: List[str]) -> str:
    """
    Generate explanation for pocket selection.
    
    Args:
        pocket: Selected pocket
        mutations: List of mutations
        
    Returns:
        Rationale string
    """
    rationale_parts = []
    
    # Druggability
    druggability = pocket.get("druggability_score", 0)
    if druggability > 0.7:
        rationale_parts.append(f"High druggability score ({druggability:.2f})")
    elif druggability > 0.5:
        rationale_parts.append(f"Moderate druggability score ({druggability:.2f})")
    else:
        rationale_parts.append(f"Low druggability score ({druggability:.2f})")
    
    # Volume
    volume = pocket.get("volume", 0)
    if 300 <= volume <= 1000:
        rationale_parts.append(f"Optimal pocket volume ({volume:.0f} Ų)")
    elif volume < 300:
        rationale_parts.append(f"Small pocket volume ({volume:.0f} Ų)")
    else:
        rationale_parts.append(f"Large pocket volume ({volume:.0f} Ų)")
    
    # Known site
    if pocket.get("is_known_site"):
        rationale_parts.append("Known binding site from literature")
    
    # Mutation proximity
    if mutations and pocket.get("description") and "G12C" in pocket.get("description", ""):
        rationale_parts.append(f"Contains target mutation site ({', '.join(mutations)})")
    
    # Hydrophobicity
    hydrophobicity = pocket.get("hydrophobicity", 0)
    if hydrophobicity > 0.5:
        rationale_parts.append("Hydrophobic character suitable for small molecules")
    
    rationale = "Selected based on: " + "; ".join(rationale_parts)
    
    return rationale
