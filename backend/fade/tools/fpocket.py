"""
fpocket interface for pocket detection.

fpocket is an open-source protein pocket detection algorithm.
It uses alpha spheres to detect cavities on protein surfaces.
"""

import logging
import subprocess
import tempfile
import shutil
import os
from pathlib import Path
from typing import Dict, Any, List, Optional
import re

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.fpocket")

# Try to import HPC fpocket client
try:
    from fade.tools.hpc_fpocket import get_hpc_fpocket_client
    HPC_AVAILABLE = True
except ImportError:
    HPC_AVAILABLE = False


class FpocketClient:
    """Client for running fpocket pocket detection."""
    
    def __init__(self, fpocket_path: str = None):
        """
        Initialize fpocket client.
        
        Args:
            fpocket_path: Path to fpocket executable
        """
        self.fpocket_path = fpocket_path or self._find_fpocket()
        self.temp_dir = Path(tempfile.gettempdir()) / "fade_fpocket"
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
    def _find_fpocket(self) -> str:
        """
        Try to find fpocket executable in common locations.
        
        Returns:
            Path to fpocket executable
        """
        # Common locations to check
        possible_paths = [
            "fpocket",  # In PATH
            "/usr/local/bin/fpocket",
            "/usr/bin/fpocket",
            "C:/fpocket/fpocket.exe",
            "C:/Program Files/fpocket/fpocket.exe",
            "./fpocket/bin/fpocket",
        ]
        
        for path in possible_paths:
            if shutil.which(path):
                logger.info(f"Found fpocket at: {path}")
                return path
                
        logger.error("fpocket not found locally")
        return None
    
    def detect_pockets(self, pdb_file: str, min_score: float = 0.5) -> List[Dict[str, Any]]:
        """
        Detect pockets in a protein structure.
        
        Args:
            pdb_file: Path to PDB file
            min_score: Minimum druggability score (0-1)
            
        Returns:
            List of detected pockets
        """
        # Try HPC first if available
        if HPC_AVAILABLE and os.getenv("HPC_USER"):
            logger.info("Using HPC fpocket")
            hpc_client = get_hpc_fpocket_client()
            return hpc_client.detect_pockets(pdb_file, min_score)
        
        # Fall back to local fpocket
        if not self.fpocket_path:
            raise RuntimeError("fpocket not available locally and HPC not configured")
        
        try:
            # Create temp directory for this run
            run_dir = self.temp_dir / f"run_{Path(pdb_file).stem}"
            run_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy PDB file to temp directory
            temp_pdb = run_dir / Path(pdb_file).name
            shutil.copy(pdb_file, temp_pdb)
            
            # Run fpocket
            logger.info(f"Running fpocket locally on {pdb_file}")
            cmd = [self.fpocket_path, "-f", str(temp_pdb)]
            
            result = subprocess.run(
                cmd,
                cwd=str(run_dir),
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"fpocket failed: {result.stderr}")
            
            # Parse fpocket output
            pockets = self._parse_fpocket_output(run_dir, temp_pdb.stem)
            
            # Filter by score
            filtered_pockets = [
                p for p in pockets 
                if p.get("druggability_score", 0) >= min_score
            ]
            
            logger.info(f"Found {len(filtered_pockets)} pockets with score >= {min_score}")
            
            # Clean up
            shutil.rmtree(run_dir, ignore_errors=True)
            
            return filtered_pockets
            
        except Exception as e:
            logger.error(f"Error running fpocket: {e}")
            raise
    
    def _parse_fpocket_output(self, run_dir: Path, pdb_stem: str) -> List[Dict[str, Any]]:
        """
        Parse fpocket output files.
        
        Args:
            run_dir: Directory containing fpocket output
            pdb_stem: PDB filename stem
            
        Returns:
            List of parsed pockets
        """
        pockets = []
        
        # Look for pocket info file
        info_file = run_dir / f"{pdb_stem}_out" / f"{pdb_stem}_info.txt"
        
        if not info_file.exists():
            raise FileNotFoundError(f"fpocket info file not found: {info_file}")
        
        # Parse the info file
        with open(info_file, 'r') as f:
            content = f.read()
        
        # Extract pocket information using regex
        pocket_pattern = r"Pocket\s+(\d+)\s*:.*?Score\s*:\s*([\d.]+).*?Druggability Score\s*:\s*([\d.]+).*?Volume\s*:\s*([\d.]+)"
        
        matches = re.finditer(pocket_pattern, content, re.DOTALL)
        
        for match in matches:
            pocket_id = int(match.group(1))
            score = float(match.group(2))
            druggability = float(match.group(3))
            volume = float(match.group(4))
            
            # Get pocket PDB file if it exists
            pocket_pdb = run_dir / f"{pdb_stem}_out" / "pockets" / f"pocket{pocket_id}_atm.pdb"
            
            pocket_info = {
                "pocket_id": f"pocket_{pocket_id}",
                "score": score,
                "druggability_score": druggability,
                "volume": volume,
                "pocket_pdb": str(pocket_pdb) if pocket_pdb.exists() else None
            }
            
            # Try to get center coordinates
            if pocket_pdb.exists():
                center = self._get_pocket_center(pocket_pdb)
                pocket_info["center"] = center
            
            pockets.append(pocket_info)
        
        return sorted(pockets, key=lambda x: x["druggability_score"], reverse=True)
    
    def _get_pocket_center(self, pocket_pdb: Path) -> List[float]:
        """
        Calculate geometric center of pocket.
        
        Args:
            pocket_pdb: Path to pocket PDB file
            
        Returns:
            [x, y, z] coordinates of center
        """
        coords = []
        
        try:
            with open(pocket_pdb, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
            
            if coords:
                # Calculate average
                center = [
                    sum(c[0] for c in coords) / len(coords),
                    sum(c[1] for c in coords) / len(coords),
                    sum(c[2] for c in coords) / len(coords)
                ]
                return center
        except Exception as e:
            logger.error(f"Error calculating pocket center: {e}")
            raise
        
        return [0.0, 0.0, 0.0]
    
    def save_pocket_pdb(self, pocket: Dict[str, Any], output_path: str) -> bool:
        """
        Save pocket as PDB file for visualization.
        
        Args:
            pocket: Pocket information
            output_path: Where to save the pocket PDB
            
        Returns:
            True if saved successfully
        """
        try:
            if pocket.get("pocket_pdb") and Path(pocket["pocket_pdb"]).exists():
                shutil.copy(pocket["pocket_pdb"], output_path)
                return True
            else:
                # Create a simple PDB with pocket center as a pseudo-atom
                center = pocket.get("center", [0, 0, 0])
                pdb_content = (
                    f"REMARK   1 POCKET {pocket['pocket_id']}\n"
                    f"REMARK   2 DRUGGABILITY SCORE: {pocket.get('druggability_score', 0):.2f}\n"
                    f"REMARK   3 VOLUME: {pocket.get('volume', 0):.1f}\n"
                    f"HETATM    1  C   PCK A   1    {center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00 99.00           C\n"
                    "END\n"
                )
                
                with open(output_path, 'w') as f:
                    f.write(pdb_content)
                
                return True
                
        except Exception as e:
            logger.error(f"Error saving pocket PDB: {e}")
            raise


# Singleton instance
_fpocket_client = None

def get_fpocket_client() -> FpocketClient:
    """Get or create fpocket client singleton."""
    global _fpocket_client
    if _fpocket_client is None:
        _fpocket_client = FpocketClient()
    return _fpocket_client
