"""
HPC-based fpocket implementation.

This runs fpocket on the Northeastern University HPC cluster
via SSH connection.
"""

import logging
import uuid
import time
from pathlib import Path
from typing import Dict, Any, List, Optional
import re

from fade.tools.ssh_client import get_hpc_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_fpocket")


class HPCFpocketClient:
    """Client for running fpocket on HPC cluster."""
    
    def __init__(self):
        """Initialize HPC fpocket client."""
        self.ssh_client = get_hpc_client()
        self.scratch_base = "/scratch/rajagopalmohanraj.n/F.A.D.E/fpocket_inputs"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles"
        
    def detect_pockets(self, pdb_file: str, min_score: float = 0.5) -> List[Dict[str, Any]]:
        """
        Detect pockets using fpocket on HPC cluster.
        
        Args:
            pdb_file: Path to local PDB file
            min_score: Minimum druggability score
            
        Returns:
            List of detected pockets
        """
        # Generate unique job ID
        job_id = str(uuid.uuid4())[:8]
        remote_dir = f"{self.scratch_base}/{job_id}"
        
        logger.info(f"Running fpocket on HPC with job ID: {job_id}")
        
        try:
            # Connect to HPC
            if not self.ssh_client.connect():
                raise ConnectionError("Failed to connect to HPC cluster")
            
            # Create remote directory
            logger.info(f"Creating remote directory: {remote_dir}")
            if not self.ssh_client.mkdir_remote(remote_dir):
                raise RuntimeError("Failed to create remote directory")
            
            # Upload PDB file
            pdb_filename = Path(pdb_file).name
            remote_pdb = f"{remote_dir}/{pdb_filename}"
            
            logger.info(f"Uploading PDB to: {remote_pdb}")
            if not self.ssh_client.upload_file(pdb_file, remote_pdb):
                raise RuntimeError("Failed to upload PDB file")
            
            # Run fpocket with module loading
            fpocket_cmd = f"""
cd {remote_dir}
module use {self.module_path}
module load fpocket/latest
fpocket -f {pdb_filename}
"""
            
            logger.info("Running fpocket on HPC...")
            stdout, stderr, exit_code = self.ssh_client.execute_command(
                fpocket_cmd,
                timeout=120  # 2 minutes timeout
            )
            
            if exit_code != 0:
                raise RuntimeError(f"fpocket failed with exit code {exit_code}: {stderr}")
            
            logger.info("fpocket completed successfully")
            
            # Parse results
            pockets = self._parse_remote_results(remote_dir, pdb_filename, min_score)
            
            if not pockets:
                raise ValueError(f"No pockets found with druggability score >= {min_score}")
            
            # Download pocket files for visualization
            self._download_pocket_files(remote_dir, pdb_filename, pockets, job_id)
            
            return pockets
            
        except Exception as e:
            logger.error(f"Error running fpocket on HPC: {e}")
            raise
        
        finally:
            self.ssh_client.close()
    
    def _parse_remote_results(self, remote_dir: str, pdb_filename: str, 
                            min_score: float) -> List[Dict[str, Any]]:
        """
        Parse fpocket results from remote server.
        
        Args:
            remote_dir: Remote directory containing results
            pdb_filename: Name of PDB file
            min_score: Minimum druggability score
            
        Returns:
            List of pockets
        """
        pdb_stem = Path(pdb_filename).stem
        info_file = f"{remote_dir}/{pdb_stem}_out/{pdb_stem}_info.txt"
        
        # Download and parse info file
        logger.info("Downloading fpocket info file...")
        stdout, stderr, exit_code = self.ssh_client.execute_command(
            f"cat {info_file}"
        )
        
        if exit_code != 0:
            raise FileNotFoundError(f"Failed to read fpocket info file: {stderr}")
        
        # Parse pocket information
        pockets = []
        
        # Split by pocket sections
        pocket_sections = re.split(r'Pocket\s+\d+\s*:', stdout)
        
        for i, section in enumerate(pocket_sections[1:], 1):  # Skip first empty section
            try:
                # Extract metrics using regex
                score_match = re.search(r'Score\s*:\s*([\d.]+)', section)
                drug_match = re.search(r'Druggability\s+Score\s*:\s*([\d.]+)', section)
                volume_match = re.search(r'Volume\s*(?:\((?:A|Angstrom)\^3\))?\s*:\s*([\d.]+)', section)
                if not volume_match:
                    volume_match = re.search(r'Volume\s*:\s*([\d.]+)', section)
                
                # Extract residues
                residues = []
                residue_matches = re.findall(r'([A-Z]{3})\s+(\d+)\s+([A-Z])', section)
                for res_match in residue_matches[:10]:  # Limit to first 10 residues
                    residues.append(f"{res_match[0]}{res_match[1]}")
                
                # Get scores
                score = float(score_match.group(1)) if score_match else 0.0
                druggability = float(drug_match.group(1)) if drug_match else 0.0
                volume = float(volume_match.group(1)) if volume_match else 0.0
                
                # Apply minimum score filter
                if druggability >= min_score:
                    # Get pocket center (alpha sphere coordinates)
                    center = self._get_pocket_center(remote_dir, pdb_stem, i)
                    
                    pocket_info = {
                        "pocket_id": f"pocket_{i}",
                        "score": score,
                        "druggability_score": druggability,
                        "volume": volume,
                        "residues": residues,
                        "center": center,
                        "remote_path": f"{remote_dir}/{pdb_stem}_out/pockets/pocket{i}_atm.pdb"
                    }
                    
                    pockets.append(pocket_info)
                    
            except Exception as e:
                logger.warning(f"Error parsing pocket {i}: {e}")
                continue
        
        # Sort by druggability score
        pockets.sort(key=lambda x: x["druggability_score"], reverse=True)
        
        logger.info(f"Parsed {len(pockets)} pockets with score >= {min_score}")
        return pockets
    
    def _get_pocket_center(self, remote_dir: str, pdb_stem: str, 
                          pocket_num: int) -> List[float]:
        """
        Get geometric center of pocket from alpha spheres.
        
        Args:
            remote_dir: Remote directory
            pdb_stem: PDB filename stem
            pocket_num: Pocket number
            
        Returns:
            [x, y, z] coordinates
        """
        # Try to get alpha sphere coordinates
        pocket_file = f"{remote_dir}/{pdb_stem}_out/pockets/pocket{pocket_num}_atm.pdb"
        
        stdout, stderr, exit_code = self.ssh_client.execute_command(
            f"head -20 {pocket_file}"
        )
        
        if exit_code == 0:
            coords = []
            for line in stdout.split('\n'):
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46]) 
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except:
                        continue
            
            if coords:
                # Calculate average
                center = [
                    sum(c[0] for c in coords) / len(coords),
                    sum(c[1] for c in coords) / len(coords),
                    sum(c[2] for c in coords) / len(coords)
                ]
                return center
        
        return [0.0, 0.0, 0.0]
    
    def _download_pocket_files(self, remote_dir: str, pdb_filename: str,
                              pockets: List[Dict[str, Any]], job_id: str):
        """
        Download pocket PDB files for visualization.
        
        Args:
            remote_dir: Remote directory
            pdb_filename: PDB filename
            pockets: List of pockets
            job_id: Job ID
        """
        # Create local directory for pockets
        pockets_dir = config.DATA_DIR / "pockets" / job_id
        pockets_dir.mkdir(parents=True, exist_ok=True)
        
        pdb_stem = Path(pdb_filename).stem
        
        for pocket in pockets[:5]:  # Download top 5 pockets
            pocket_num = int(pocket["pocket_id"].split("_")[1])
            remote_pocket = f"{remote_dir}/{pdb_stem}_out/pockets/pocket{pocket_num}_atm.pdb"
            local_pocket = pockets_dir / f"{pocket['pocket_id']}.pdb"
            
            if self.ssh_client.download_file(remote_pocket, str(local_pocket)):
                pocket["local_file"] = str(local_pocket)
                logger.debug(f"Downloaded pocket file: {local_pocket}")


# Singleton instance
_hpc_fpocket_client = None

def get_hpc_fpocket_client() -> HPCFpocketClient:
    """Get or create HPC fpocket client singleton."""
    global _hpc_fpocket_client
    if _hpc_fpocket_client is None:
        _hpc_fpocket_client = HPCFpocketClient()
    return _hpc_fpocket_client
