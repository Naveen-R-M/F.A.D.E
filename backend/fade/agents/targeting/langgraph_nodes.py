"""
Enhanced targeting module LangGraph nodes with consensus detection.

These nodes handle pocket detection and ranking:
1. pocket_mapper_node: Finds pockets using single tool or consensus
2. pocket_ranker_node: Ranks and selects best pocket
"""

from typing import Dict, Any, List, Optional
from pathlib import Path

from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.tools.fpocket import get_fpocket_client
from fade.tools.consensus import get_consensus_detector
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.targeting")


def determine_pocket_detection_strategy(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Determine optimal pocket detection strategy based on structure source.
    
    Args:
        state: Current pipeline state
        
    Returns:
        Dictionary with strategy details
    """
    structure = state.get("structure", {})
    structure_source = structure.get("source", "unknown")
    has_ligand = structure.get("has_ligand", False)
    confidence = structure.get("confidence", 100)
    
    strategy = {
        "use_consensus": False,
        "tools": [],
        "min_agreement": 2,
        "structure_source": "pdb_apo",  # Default
        "reason": ""
    }
    
    # Determine structure source type for consensus
    if structure_source == "pdb" and has_ligand:
        strategy["structure_source"] = "pdb_holo"
        strategy["tools"] = ["fpocket"]
        strategy["use_consensus"] = False
        strategy["reason"] = "Experimental structure with known binding site"
        
    elif structure_source == "pdb" and not has_ligand:
        strategy["structure_source"] = "pdb_apo"
        strategy["tools"] = ["fpocket", "p2rank"]
        strategy["use_consensus"] = True
        strategy["min_agreement"] = 2
        strategy["reason"] = "Apo experimental structure needs validation"
        
    elif structure_source == "alphafold":
        strategy["structure_source"] = "alphafold"
        strategy["tools"] = ["fpocket", "p2rank", "kalasanty"]
        strategy["use_consensus"] = True
        strategy["min_agreement"] = 2
        strategy["reason"] = "AlphaFold structure benefits from multiple predictions"
        
    elif structure_source == "boltz2":
        strategy["structure_source"] = "boltz2"
        strategy["tools"] = ["fpocket", "p2rank", "kalasanty"]
        strategy["use_consensus"] = True
        strategy["min_agreement"] = 2  # Changed from 3 to be more flexible
        strategy["reason"] = "Boltz2 structure needs maximum validation"
        
    else:
        strategy["tools"] = ["fpocket", "p2rank"]
        strategy["use_consensus"] = True
        strategy["min_agreement"] = 2
        strategy["reason"] = "Default strategy for unknown structure source"
    
    # Adjust based on confidence
    if confidence < 70 and "kalasanty" not in strategy["tools"]:
        strategy["tools"].append("kalasanty")
        strategy["use_consensus"] = True
        strategy["reason"] += f" (low confidence: {confidence}%)"
    
    return strategy


def pocket_mapper_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Enhanced LangGraph node for detecting binding pockets.
    
    Now includes consensus detection for improved accuracy.
    
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
    
    # Determine detection strategy
    strategy = determine_pocket_detection_strategy(state)
    
    logger.info(f"Detection strategy: {strategy['reason']}")
    logger.info(f"Using tools: {', '.join(strategy['tools'])}")
    
    # Get or generate job_id
    job_id = state.get("job_id")
    if not job_id:
        import uuid
        job_id = str(uuid.uuid4())[:8]
        logger.info(f"Generated new job_id: {job_id}")
    
    # Run pocket detection
    if strategy["use_consensus"]:
        # Use consensus detection
        logger.info("Running consensus pocket detection...")
        
        consensus_detector = get_consensus_detector()
        
        # Run consensus detection
        consensus_pockets = consensus_detector.detect_consensus(
            pdb_file=structure_path,
            structure_source=strategy["structure_source"],
            tools=strategy["tools"],
            min_agreement=strategy["min_agreement"]
        )
        
        # Convert to standard format
        pockets = []
        for cp in consensus_pockets:
            pocket = {
                "pocket_id": cp.get("pocket_id"),
                "center": cp.get("center"),
                "residues": cp.get("residues", []),
                "volume": cp.get("volume", 0),
                "druggability_score": cp.get("druggability_score", 0.5),
                "confidence_score": cp.get("consensus_score", 0.5),
                "agreement_level": cp.get("agreement_level", 1),
                "source_tools": cp.get("source_tools", []),
                "is_consensus": True,
                "center_variance": cp.get("center_variance", 0)
            }
            pockets.append(pocket)
        
    else:
        # Use single tool (fpocket)
        logger.info("Running single-tool pocket detection (fpocket)...")
        
        fpocket_client = get_fpocket_client()
        
        # Detect pockets
        raw_pockets = fpocket_client.detect_pockets(
            pdb_file=structure_path,
            min_score=0.5
        )
        
        # Extract job_id if it was set by HPC fpocket
        if raw_pockets and raw_pockets[0].get("remote_path"):
            import re
            match = re.search(r'/fpocket_inputs/([^/]+)/', raw_pockets[0]["remote_path"])
            if match:
                job_id = match.group(1)
                logger.info(f"Updated job_id from fpocket: {job_id}")
        
        # Convert to standard format
        pockets = []
        for rp in raw_pockets:
            pocket = {
                "pocket_id": rp.get("pocket_id", f"pocket_{len(pockets)+1}"),
                "center": rp.get("center", (0, 0, 0)),
                "residues": rp.get("residues", []),
                "volume": rp.get("volume", 0),
                "druggability_score": rp.get("druggability_score", 0.5),
                "confidence_score": rp.get("druggability_score", 0.5),
                "agreement_level": 1,
                "source_tools": ["fpocket"],
                "is_consensus": False
            }
            pockets.append(pocket)
    
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
    
    # Save pocket details for visualization
    for pocket in pockets:
        pocket_file = pockets_dir / f"{pocket['pocket_id']}.txt"
        with open(pocket_file, 'w') as f:
            f.write(f"Pocket ID: {pocket['pocket_id']}\n")
            f.write(f"Center: {pocket['center']}\n")
            f.write(f"Volume: {pocket.get('volume', 0):.1f} Ų\n")
            f.write(f"Druggability: {pocket.get('druggability_score', 0):.3f}\n")
            
            if pocket.get("is_consensus"):
                f.write(f"Consensus Score: {pocket.get('confidence_score', 0):.3f}\n")
                f.write(f"Agreement: {pocket.get('agreement_level')} tools\n")
                f.write(f"Tools: {', '.join(pocket.get('source_tools', []))}\n")
                if pocket.get("center_variance"):
                    f.write(f"Spatial Variance: {pocket['center_variance']:.2f} Ų\n")
            else:
                f.write(f"Confidence: {pocket.get('confidence_score', 0):.3f}\n")
                f.write(f"Tool: {', '.join(pocket.get('source_tools', []))}\n")
            
            if pocket.get("residues"):
                f.write(f"Residues ({len(pocket['residues'])}): ")
                f.write(', '.join(pocket['residues'][:20]))
                if len(pocket['residues']) > 20:
                    f.write(f" ... and {len(pocket['residues']) - 20} more")
                f.write("\n")
        
        pocket["pocket_file"] = str(pocket_file)
    
    # Create summary message
    message = f"Found {len(pockets)} potential binding pockets"
    
    if pockets:
        best_pocket = pockets[0]
        message += f"\nBest pocket: {best_pocket['pocket_id']}"
        
        if best_pocket.get("is_consensus"):
            message += f" (consensus from {best_pocket['agreement_level']} tools: "
            message += f"{', '.join(best_pocket.get('source_tools', []))})"
        
        message += f"\n  Druggability: {best_pocket.get('druggability_score', 0):.2f}"
        message += f"\n  Confidence: {best_pocket.get('confidence_score', 0):.2f}"
        message += f"\n  Volume: {best_pocket.get('volume', 0):.1f} Ų"
        
        if best_pocket.get("center_variance") is not None and best_pocket.get("is_consensus"):
            message += f"\n  Spatial consistency: {10 - min(10, best_pocket['center_variance']):.1f}/10"
    
    logger.info(message)
    
    # Store detection metadata
    pocket_detection_metadata = {
        "strategy": strategy,
        "consensus_used": strategy["use_consensus"],
        "tools_used": strategy["tools"],
        "total_pockets_found": len(pockets)
    }
    
    return {
        "pockets": pockets,
        "pocket_detection_metadata": pocket_detection_metadata,
        "messages": [HumanMessage(content=message)],
        "current_step": "pockets_detected",
        "should_continue": True,
        "job_id": job_id
    }


def pocket_ranker_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Enhanced LangGraph node for ranking and selecting the best pocket.
    
    Now considers consensus scores and agreement levels.
    
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
        score = _calculate_enhanced_pocket_score(pocket, mutations)
        pocket_with_score = pocket.copy()
        pocket_with_score["final_score"] = score
        scored_pockets.append(pocket_with_score)
    
    # Sort by score
    scored_pockets.sort(key=lambda x: x["final_score"], reverse=True)
    
    # Select the best pocket
    selected_pocket = scored_pockets[0]
    
    # Generate rationale
    rationale = _generate_enhanced_selection_rationale(selected_pocket, mutations)
    
    message = f"Selected pocket: {selected_pocket['pocket_id']}"
    
    if selected_pocket.get("is_consensus"):
        message += f" (consensus from {selected_pocket['agreement_level']} tools)"
    
    message += f"\nDruggability score: {selected_pocket.get('druggability_score', 0):.2f}"
    message += f"\nConfidence: {selected_pocket.get('confidence_score', 0):.2f}"
    message += f"\nVolume: {selected_pocket.get('volume', 0):.1f} Ų"
    message += f"\nFinal score: {selected_pocket['final_score']:.2f}"
    
    logger.info(f"Selected {selected_pocket['pocket_id']} as best pocket")
    
    return {
        "pockets": scored_pockets,  # Update with scores
        "selected_pocket": selected_pocket,
        "pocket_selection_rationale": rationale,
        "messages": [HumanMessage(content=message)],
        "current_step": "pocket_selected",
        "should_continue": True
    }


def _calculate_enhanced_pocket_score(pocket: Dict[str, Any], mutations: List[str]) -> float:
    """
    Calculate enhanced score including consensus metrics.
    
    Args:
        pocket: Pocket information
        mutations: List of mutations
        
    Returns:
        Final score (0-100)
    """
    score = 0.0
    
    # Base scoring (60% weight)
    # Druggability score (25% weight)
    druggability = pocket.get("druggability_score", 0.5)
    score += druggability * 25
    
    # Confidence/Consensus score (25% weight)
    confidence = pocket.get("confidence_score", 0.5)
    score += confidence * 25
    
    # Volume score (10% weight)
    volume = pocket.get("volume", 500)
    if 300 <= volume <= 1000:
        volume_score = 1.0
    elif volume < 300:
        volume_score = volume / 300
    else:
        volume_score = max(0, 2 - volume / 1000)
    score += volume_score * 10
    
    # Consensus bonus (20% weight)
    if pocket.get("is_consensus"):
        agreement_level = pocket.get("agreement_level", 1)
        consensus_bonus = (agreement_level / 3) * 20  # Max 3 tools
        score += consensus_bonus
        
        # Additional bonus for low spatial variance
        if pocket.get("center_variance") is not None:
            variance_bonus = max(0, 5 * (1 - pocket["center_variance"] / 10))
            score += variance_bonus
    else:
        # Single tool gets partial credit
        score += 5
    
    # Mutation proximity bonus (15% weight)
    if mutations and pocket.get("residues"):
        mutation_residues = set()
        for mutation in mutations:
            import re
            match = re.match(r'[A-Z](\d+)[A-Z]', mutation)
            if match:
                residue_num = match.group(1)
                mutation_residues.add(residue_num)
        
        pocket_residues = set()
        for res in pocket.get("residues", []):
            match = re.match(r'[A-Z]*(\d+)', res)
            if match:
                pocket_residues.add(match.group(1))
        
        if mutation_residues & pocket_residues:
            score += 15
            logger.info(f"Pocket {pocket['pocket_id']} contains mutation site(s)")
    
    return min(100, score)


def _generate_enhanced_selection_rationale(pocket: Dict[str, Any], mutations: List[str]) -> str:
    """
    Generate enhanced explanation including consensus information.
    
    Args:
        pocket: Selected pocket
        mutations: List of mutations
        
    Returns:
        Rationale string
    """
    rationale_parts = []
    
    # Consensus information
    if pocket.get("is_consensus"):
        agreement = pocket.get("agreement_level", 1)
        tools = pocket.get("source_tools", [])
        rationale_parts.append(f"Consensus pocket from {agreement} tools ({', '.join(tools)})")
        
        if pocket.get("center_variance") is not None:
            if pocket["center_variance"] < 3:
                rationale_parts.append("Excellent spatial agreement between tools")
            elif pocket["center_variance"] < 7:
                rationale_parts.append("Good spatial agreement between tools")
    else:
        tool = pocket.get("source_tools", ["unknown"])[0]
        rationale_parts.append(f"Detected by {tool}")
    
    # Druggability
    druggability = pocket.get("druggability_score", 0)
    if druggability > 0.7:
        rationale_parts.append(f"High druggability ({druggability:.2f})")
    elif druggability > 0.5:
        rationale_parts.append(f"Moderate druggability ({druggability:.2f})")
    else:
        rationale_parts.append(f"Low druggability ({druggability:.2f})")
    
    # Confidence
    confidence = pocket.get("confidence_score", 0)
    if confidence > 0.8:
        rationale_parts.append(f"High confidence ({confidence:.2f})")
    elif confidence > 0.6:
        rationale_parts.append(f"Good confidence ({confidence:.2f})")
    
    # Volume
    volume = pocket.get("volume", 0)
    if 300 <= volume <= 1000:
        rationale_parts.append(f"Optimal volume ({volume:.0f} Ų)")
    elif volume < 300:
        rationale_parts.append(f"Small volume ({volume:.0f} Ų)")
    else:
        rationale_parts.append(f"Large volume ({volume:.0f} Ų)")
    
    # Mutation proximity
    if mutations and pocket.get("residues"):
        for mutation in mutations:
            import re
            match = re.match(r'[A-Z](\d+)[A-Z]', mutation)
            if match:
                residue_num = match.group(1)
                if any(residue_num in str(res) for res in pocket.get("residues", [])):
                    rationale_parts.append(f"Contains mutation site {mutation}")
                    break
    
    rationale = "Selected based on: " + "; ".join(rationale_parts)
    
    return rationale


def pocket_from_ligand_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    LangGraph node to extract pocket from bound ligand in holo structures.
    
    This node runs ONLY for structures with drug-like ligands.
    It extracts the pocket definition from the ligand's binding site.
    
    Args:
        state: Current pipeline state
        
    Returns:
        State updates with selected pocket from ligand
    """
    logger.info("[Pocket from Ligand] Extracting pocket from bound ligand...")
    
    structure = state.get("structure", {})
    
    if not structure or not structure.get("structure_path"):
        return {
            "error": "No structure available for pocket extraction",
            "should_continue": False,
            "current_step": "pocket_extraction_failed"
        }
    
    if not structure.get("has_drug_like_ligand"):
        return {
            "error": "No drug-like ligand found in structure",
            "should_continue": False,
            "current_step": "pocket_extraction_failed"
        }
    
    try:
        # Import the utility function
        from fade.tools.pdb_pocket_utils import extract_ligand_pocket, get_first_drug_like_ligand
        
        # Get the PDB path and ligand ID
        pdb_path = structure["structure_path"]
        ligand_id = get_first_drug_like_ligand(structure)
        
        if not ligand_id:
            return {
                "error": "Could not identify drug-like ligand ID",
                "should_continue": False,
                "current_step": "pocket_extraction_failed"
            }
        
        logger.info(f"Extracting pocket from ligand: {ligand_id}")
        
        # Extract pocket information
        pocket = extract_ligand_pocket(pdb_path, ligand_id)
        
        # Add additional metadata
        pocket.update({
            "pocket_id": "ligand_pocket",
            "rank": 1,
            "confidence_score": 1.0,  # Highest confidence - actual binding site
            "druggability_score": 1.0,  # Known to bind drugs
            "source_tools": ["ligand_extraction"],
            "is_consensus": False,
            "rationale": f"Pocket defined by bound drug-like ligand {ligand_id}. This is the actual binding site."
        })
        
        logger.info(f"Successfully extracted pocket from ligand {ligand_id}")
        logger.info(f"  Center: {pocket['center']}")
        logger.info(f"  Residues: {len(pocket['residues'])} residues")
        logger.info(f"  Volume: {pocket.get('pocket_volume', 0):.0f} Ų")
        
        message = f"✓ Extracted pocket from bound ligand {ligand_id} - proceeding to molecule generation"
        
        return {
            "selected_pocket": pocket,
            "pockets": [pocket],  # Store as single-item list for consistency
            "messages": [HumanMessage(content=message)],
            "current_step": "pocket_extracted",
            "should_continue": True
        }
        
    except Exception as e:
        logger.error(f"Failed to extract pocket from ligand: {e}")
        return {
            "error": f"Pocket extraction failed: {str(e)}",
            "should_continue": False,
            "current_step": "pocket_extraction_failed"
        }
