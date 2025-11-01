"""
HPC Kalasanty client for protein pocket prediction.

Kalasanty is a deep learning-based tool for predicting ligand binding sites
from protein structure, optimized for AlphaFold structures.
"""

import os
import time
import json
import re
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import tempfile

from fade.tools.ssh_client import get_hpc_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_kalasanty")


class HPCKalasantyClient:
    """Client for running Kalasanty pocket prediction on HPC."""
    
    def __init__(self):
        """Initialize HPC Kalasanty client."""
        self.ssh_client = get_hpc_client()
        self.scratch_base = "/scratch/rajagopalmohanraj.n/F.A.D.E/kalasanty_jobs"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles/"
        self.kalasanty_module = "kalasanty/latest"
        
    def detect_pockets(self, pdb_file: str, min_confidence: float = 0.5) -> List[Dict[str, Any]]:
        """
        Detect pockets using Kalasanty on HPC cluster.
        
        Args:
            pdb_file: Path to local PDB file
            min_confidence: Minimum confidence score
            
        Returns:
            List of detected pockets
        """
        import uuid
        job_id = str(uuid.uuid4())[:8]
        remote_dir = f"{self.scratch_base}/{job_id}"
        
        logger.info(f"Running Kalasanty on HPC with job ID: {job_id}")
        
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
            
            # Run Kalasanty
            kalasanty_cmd = self._build_kalasanty_command(
                pdb_filename, remote_dir
            )
            
            logger.info("Running Kalasanty on HPC...")
            stdout, stderr, exit_code = self.ssh_client.execute_command(
                kalasanty_cmd,
                timeout=300  # 5 minutes timeout
            )
            
            if exit_code != 0:
                logger.error(f"Kalasanty stderr: {stderr}")
                raise RuntimeError(f"Kalasanty failed with exit code {exit_code}")
            
            logger.info("Kalasanty completed successfully")
            
            # Parse results
            pockets = self._parse_results(remote_dir, pdb_filename, min_confidence)
            
            if not pockets:
                logger.warning(f"No pockets found with confidence >= {min_confidence}")
            
            return pockets
            
        except Exception as e:
            logger.error(f"Error running Kalasanty on HPC: {e}")
            # Return empty list instead of raising to allow pipeline to continue
            return []
        
        finally:
            self.ssh_client.disconnect()
    
    def _build_kalasanty_command(self, pdb_filename: str, remote_dir: str) -> str:
        """
        Build the Kalasanty command.
        
        Args:
            pdb_filename: Name of PDB file
            remote_dir: Remote working directory
            
        Returns:
            Command string
        """
        # Based on testing, kalasanty command syntax
        # This will be updated based on actual test results
        return f"""
cd {remote_dir}
module use {self.module_path}
module load {self.kalasanty_module}

# Run Kalasanty
kalasanty -i {pdb_filename} -o kalasanty_output/ --format json --confidence 0.3

# If above doesn't work, try alternative
if [ $? -ne 0 ]; then
    python -m kalasanty predict --pdb {pdb_filename} --out kalasanty_output/
fi

# Check output
ls -la kalasanty_output/
"""
    
    def _parse_results(self, remote_dir: str, pdb_filename: str, 
                      min_confidence: float) -> List[Dict[str, Any]]:
        """
        Parse Kalasanty output files.
        
        Args:
            remote_dir: Remote directory containing output
            pdb_filename: Original PDB filename
            min_confidence: Minimum confidence threshold
            
        Returns:
            List of pocket dictionaries
        """
        pockets = []
        
        # Try different output locations/formats
        output_files = [
            f"{remote_dir}/kalasanty_output/pockets.json",
            f"{remote_dir}/kalasanty_output/predictions.json",
            f"{remote_dir}/kalasanty_output/{pdb_filename}.pockets.json",
            f"{remote_dir}/kalasanty_output/results.json",
        ]
        
        for output_file in output_files:
            # Check if file exists
            stdout, stderr, exit_code = self.ssh_client.execute_command(
                f"cat {output_file} 2>/dev/null"
            )
            
            if exit_code == 0 and stdout:
                try:
                    data = json.loads(stdout)
                    pockets = self._extract_pockets_from_json(data, min_confidence)
                    if pockets:
                        logger.info(f"Found {len(pockets)} pockets in {output_file}")
                        break
                except json.JSONDecodeError:
                    logger.debug(f"Failed to parse JSON from {output_file}")
        
        # If no JSON found, try CSV or text format
        if not pockets:
            pockets = self._try_alternative_formats(remote_dir, min_confidence)
        
        # Convert to standard format
        formatted_pockets = []
        for i, pocket in enumerate(pockets, 1):
            formatted_pocket = {
                "pocket_id": f"kalasanty_pocket_{i}",
                "center": pocket.get("center", (0, 0, 0)),
                "confidence_score": pocket.get("confidence", 0.5),
                "volume": pocket.get("volume", 0),
                "residues": pocket.get("residues", []),
                "source": "kalasanty",
                "remote_path": remote_dir
            }
            formatted_pockets.append(formatted_pocket)
        
        return formatted_pockets
    
    def _extract_pockets_from_json(self, data: Dict, min_confidence: float) -> List[Dict]:
        """Extract pockets from JSON data."""
        pockets = []
        
        # Handle different JSON structures
        if "pockets" in data:
            pocket_list = data["pockets"]
        elif "predictions" in data:
            pocket_list = data["predictions"]
        elif "results" in data:
            pocket_list = data["results"]
        else:
            # Assume data itself is a list
            pocket_list = data if isinstance(data, list) else []
        
        for pocket in pocket_list:
            confidence = pocket.get("confidence", pocket.get("score", 0))
            if confidence >= min_confidence:
                # Extract center coordinates
                if "center" in pocket:
                    if isinstance(pocket["center"], dict):
                        center = (
                            pocket["center"].get("x", 0),
                            pocket["center"].get("y", 0),
                            pocket["center"].get("z", 0)
                        )
                    else:
                        center = tuple(pocket["center"])
                else:
                    center = (
                        pocket.get("x", 0),
                        pocket.get("y", 0),
                        pocket.get("z", 0)
                    )
                
                pockets.append({
                    "center": center,
                    "confidence": confidence,
                    "volume": pocket.get("volume", 0),
                    "residues": pocket.get("residues", [])
                })
        
        return pockets
    
    def _try_alternative_formats(self, remote_dir: str, 
                                 min_confidence: float) -> List[Dict]:
        """Try parsing CSV or text formats if JSON not found."""
        pockets = []
        
        # Try CSV format
        csv_file = f"{remote_dir}/kalasanty_output/pockets.csv"
        stdout, stderr, exit_code = self.ssh_client.execute_command(
            f"cat {csv_file} 2>/dev/null"
        )
        
        if exit_code == 0 and stdout:
            pockets = self._parse_csv_output(stdout, min_confidence)
        
        # Try text format if CSV didn't work
        if not pockets:
            txt_file = f"{remote_dir}/kalasanty_output/pockets.txt"
            stdout, stderr, exit_code = self.ssh_client.execute_command(
                f"cat {txt_file} 2>/dev/null"
            )
            
            if exit_code == 0 and stdout:
                pockets = self._parse_text_output(stdout, min_confidence)
        
        return pockets
    
    def _parse_csv_output(self, csv_content: str, min_confidence: float) -> List[Dict]:
        """Parse CSV format output."""
        import csv
        from io import StringIO
        
        pockets = []
        try:
            reader = csv.DictReader(StringIO(csv_content))
            for row in reader:
                confidence = float(row.get("confidence", row.get("score", 0)))
                if confidence >= min_confidence:
                    pockets.append({
                        "center": (
                            float(row.get("x", 0)),
                            float(row.get("y", 0)),
                            float(row.get("z", 0))
                        ),
                        "confidence": confidence,
                        "volume": float(row.get("volume", 0)),
                        "residues": row.get("residues", "").split()
                    })
        except Exception as e:
            logger.debug(f"Failed to parse CSV: {e}")
        
        return pockets
    
    def _parse_text_output(self, txt_content: str, min_confidence: float) -> List[Dict]:
        """Parse text format output."""
        pockets = []
        
        # Simple text parsing
        lines = txt_content.strip().split('\n')
        current_pocket = {}
        
        for line in lines:
            if "pocket" in line.lower() and current_pocket:
                # Save previous pocket if confidence is high enough
                if current_pocket.get("confidence", 0) >= min_confidence:
                    pockets.append(current_pocket)
                current_pocket = {}
            
            # Extract data from line
            if "center:" in line.lower():
                match = re.search(r'([-\d.]+)[,\s]+([-\d.]+)[,\s]+([-\d.]+)', line)
                if match:
                    current_pocket["center"] = (
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3))
                    )
            elif "confidence:" in line.lower() or "score:" in line.lower():
                match = re.search(r'([\d.]+)', line)
                if match:
                    current_pocket["confidence"] = float(match.group(1))
            elif "volume:" in line.lower():
                match = re.search(r'([\d.]+)', line)
                if match:
                    current_pocket["volume"] = float(match.group(1))
        
        # Don't forget last pocket
        if current_pocket and current_pocket.get("confidence", 0) >= min_confidence:
            pockets.append(current_pocket)
        
        return pockets


# Factory function
def get_kalasanty_client() -> HPCKalasantyClient:
    """Get configured Kalasanty HPC client."""
    return HPCKalasantyClient()
