"""
Consensus pocket detection system for F.A.D.E.

Combines predictions from multiple pocket detection tools (fpocket, P2Rank, Kalasanty)
to improve accuracy and reliability of binding site identification.
"""

import os
import time
import numpy as np
from typing import Dict, Any, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import logging
import tempfile

from fade.tools.hpc_fpocket import HPCFpocketClient
from fade.tools.hpc_p2rank import HPCP2RankClient
from fade.tools.hpc_kalasanty import HPCKalasantyClient
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.consensus")


@dataclass
class UnifiedPocket:
    """Standard pocket format for all tools."""
    pocket_id: str
    center: Tuple[float, float, float]  # x, y, z
    residues: List[str]
    volume: float  # Ų
    
    # Scores (0-1 normalized)
    druggability_score: Optional[float] = None
    confidence_score: Optional[float] = None
    consensus_score: Optional[float] = None
    
    # Metadata
    source_tools: List[str] = field(default_factory=list)
    tool_scores: Dict[str, float] = field(default_factory=dict)
    
    # For consensus
    agreement_level: int = 1  # Number of tools that found it
    is_consensus: bool = False
    center_variance: Optional[float] = None


class ConsensusPocketDetector:
    """Consensus-based pocket detection using multiple tools."""
    
    # Tool weights based on structure source
    TOOL_WEIGHTS = {
        'fpocket': {
            'pdb_holo': 1.0,     # Best for known binding sites
            'pdb_apo': 0.8,      # Good for experimental structures
            'alphafold': 0.6,    # Moderate for predicted
            'boltz2': 0.5        # Lower for highly predicted
        },
        'p2rank': {
            'pdb_holo': 0.7,
            'pdb_apo': 0.9,      # Excellent for apo structures
            'alphafold': 0.8,
            'boltz2': 0.7
        },
        'kalasanty': {
            'pdb_holo': 0.6,
            'pdb_apo': 0.8,
            'alphafold': 0.9,    # Optimized for AF structures
            'boltz2': 0.8
        }
    }
    
    def __init__(self, grid_size: float = 5.0):
        """
        Initialize consensus detector.
        
        Args:
            grid_size: Grid cell size in Angstroms for spatial clustering
        """
        self.fpocket_client = HPCFpocketClient()
        self.p2rank_client = HPCP2RankClient()
        self.kalasanty_client = HPCKalasantyClient()
        self.grid_size = grid_size
        
    def detect_consensus(
        self,
        pdb_file: str,
        structure_source: str = "pdb_apo",
        tools: List[str] = None,
        min_agreement: int = 2
    ) -> List[Dict[str, Any]]:
        """
        Detect pockets using consensus from multiple tools.
        
        Args:
            pdb_file: Path to PDB file
            structure_source: Source type ("pdb_holo", "pdb_apo", "alphafold", "boltz2")
            tools: List of tools to use (if None, auto-select based on source)
            min_agreement: Minimum number of tools that must agree
            
        Returns:
            List of consensus pockets sorted by score
        """
        logger.info(f"Starting consensus detection for {structure_source} structure")
        
        # Auto-select tools if not specified
        if tools is None:
            tools = self._select_tools(structure_source)
        
        logger.info(f"Using tools: {', '.join(tools)}")
        
        # Run tools in parallel
        tool_results = self._run_parallel_detection(pdb_file, tools)
        
        if len(tool_results) < min_agreement:
            logger.warning(f"Only {len(tool_results)} tools succeeded, need {min_agreement}")
            if len(tool_results) == 0:
                return []
        
        # Cluster pockets spatially
        clusters = self._cluster_pockets(tool_results)
        
        # Create consensus pockets
        consensus_pockets = []
        for cluster in clusters:
            if len(cluster['tools']) >= min_agreement:
                consensus_pocket = self._create_consensus_pocket(
                    cluster, structure_source
                )
                consensus_pockets.append(consensus_pocket)
        
        # Sort by consensus score
        consensus_pockets.sort(key=lambda p: p.get('consensus_score', 0), reverse=True)
        
        logger.info(f"Generated {len(consensus_pockets)} consensus pockets")
        
        return consensus_pockets
    
    def _select_tools(self, structure_source: str) -> List[str]:
        """Select appropriate tools based on structure source."""
        if structure_source == "pdb_holo":
            return ["fpocket"]  # Single tool sufficient for holo structures
        elif structure_source == "pdb_apo":
            return ["fpocket", "p2rank"]
        elif structure_source == "alphafold":
            return ["fpocket", "p2rank", "kalasanty"]
        elif structure_source == "boltz2":
            return ["fpocket", "p2rank", "kalasanty"]
        else:
            return ["fpocket", "p2rank"]  # Default
    
    def _run_parallel_detection(
        self,
        pdb_file: str,
        tools: List[str]
    ) -> Dict[str, List[Dict]]:
        """Run multiple tools in parallel."""
        results = {}
        
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {}
            
            for tool in tools:
                if tool == "fpocket":
                    future = executor.submit(
                        self.fpocket_client.detect_pockets,
                        pdb_file, 0.3  # Lower threshold for consensus
                    )
                elif tool == "p2rank":
                    # P2Rank expects PDB content, not file path
                    with open(pdb_file, 'r') as f:
                        pdb_content = f.read()
                    future = executor.submit(
                        self._run_p2rank,
                        pdb_content
                    )
                elif tool == "kalasanty":
                    future = executor.submit(
                        self.kalasanty_client.detect_pockets,
                        pdb_file, 0.3
                    )
                else:
                    continue
                
                futures[future] = tool
            
            # Collect results
            for future in as_completed(futures, timeout=300):
                tool = futures[future]
                try:
                    tool_results = future.result()
                    if tool_results:
                        results[tool] = tool_results
                        logger.info(f"{tool} found {len(tool_results)} pockets")
                except Exception as e:
                    logger.error(f"Error running {tool}: {e}")
        
        return results
    
    def _run_p2rank(self, pdb_content: str) -> List[Dict]:
        """Run P2Rank and return results."""
        result = self.p2rank_client.predict_pockets(pdb_content)
        return result.get("pockets", [])
    
    def _cluster_pockets(self, tool_results: Dict[str, List[Dict]]) -> List[Dict]:
        """Cluster pockets from different tools based on spatial proximity."""
        clusters = []
        all_pockets = []
        
        # Collect all pockets with tool labels
        for tool, pockets in tool_results.items():
            for pocket in pockets:
                all_pockets.append({
                    'tool': tool,
                    'pocket': pocket
                })
        
        # Simple distance-based clustering
        used = set()
        for i, pocket1 in enumerate(all_pockets):
            if i in used:
                continue
            
            cluster = {
                'pockets': [pocket1],
                'tools': {pocket1['tool']},
                'centers': [pocket1['pocket'].get('center', (0, 0, 0))]
            }
            used.add(i)
            
            # Find nearby pockets from other tools
            for j, pocket2 in enumerate(all_pockets):
                if j in used or pocket2['tool'] == pocket1['tool']:
                    continue
                
                # Calculate distance between centers
                center1 = pocket1['pocket'].get('center', (0, 0, 0))
                center2 = pocket2['pocket'].get('center', (0, 0, 0))
                
                distance = np.sqrt(sum((a - b)**2 for a, b in zip(center1, center2)))
                
                # If within threshold, add to cluster
                if distance < self.grid_size * 2:  # 10Å threshold
                    cluster['pockets'].append(pocket2)
                    cluster['tools'].add(pocket2['tool'])
                    cluster['centers'].append(center2)
                    used.add(j)
            
            clusters.append(cluster)
        
        return clusters
    
    def _create_consensus_pocket(
        self,
        cluster: Dict,
        structure_source: str
    ) -> Dict[str, Any]:
        """Create a consensus pocket from a cluster."""
        # Calculate consensus center (weighted average)
        centers = cluster['centers']
        weights = []
        for pocket_data in cluster['pockets']:
            tool = pocket_data['tool']
            weight = self.TOOL_WEIGHTS[tool].get(structure_source, 0.5)
            weights.append(weight)
        
        # Weighted average center
        weighted_centers = [
            [c[i] * w for i in range(3)]
            for c, w in zip(centers, weights)
        ]
        consensus_center = tuple(
            sum(wc[i] for wc in weighted_centers) / sum(weights)
            for i in range(3)
        )
        
        # Calculate center variance
        center_variance = np.var([
            np.linalg.norm(np.array(c) - consensus_center)
            for c in centers
        ])
        
        # Collect all residues
        all_residues = set()
        for pocket_data in cluster['pockets']:
            residues = pocket_data['pocket'].get('residues', [])
            all_residues.update(residues)
        
        # Average scores
        scores = {}
        for pocket_data in cluster['pockets']:
            pocket = pocket_data['pocket']
            tool = pocket_data['tool']
            
            # Druggability
            if 'druggability_score' in pocket:
                scores.setdefault('druggability', []).append(
                    pocket['druggability_score']
                )
            
            # Confidence
            if 'confidence_score' in pocket:
                scores.setdefault('confidence', []).append(
                    pocket['confidence_score']
                )
            elif 'probability' in pocket:
                scores.setdefault('confidence', []).append(
                    pocket['probability']
                )
        
        # Calculate consensus score
        agreement_score = len(cluster['tools']) / 3  # Assume max 3 tools
        avg_druggability = np.mean(scores.get('druggability', [0.5]))
        avg_confidence = np.mean(scores.get('confidence', [0.5]))
        
        consensus_score = (
            0.35 * agreement_score +
            0.30 * avg_confidence +
            0.25 * avg_druggability +
            0.10 * (1 - min(1, center_variance / 10))  # Spatial consistency
        )
        
        # Apply structure source modifier
        source_modifiers = {
            'pdb_holo': 1.0,
            'pdb_apo': 0.95,
            'alphafold': 0.85,
            'boltz2': 0.80
        }
        consensus_score *= source_modifiers.get(structure_source, 0.9)
        
        return {
            'pocket_id': f"consensus_{len(cluster['tools'])}_{hash(consensus_center) % 10000}",
            'center': consensus_center,
            'residues': list(all_residues),
            'druggability_score': avg_druggability,
            'confidence_score': avg_confidence,
            'consensus_score': min(1.0, consensus_score),
            'agreement_level': len(cluster['tools']),
            'source_tools': list(cluster['tools']),
            'is_consensus': True,
            'center_variance': center_variance
        }


def get_consensus_detector() -> ConsensusPocketDetector:
    """Get configured consensus detector."""
    return ConsensusPocketDetector()
