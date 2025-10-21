"""
HPC-based DiffSBDD implementation for molecule generation.

This runs DiffSBDD on the Northeastern University HPC cluster
to generate novel drug-like molecules for binding pockets.
"""

import logging
import uuid
import time
import re
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import json

from fade.tools.ssh_client import get_hpc_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_diffsbdd")


class HPCDiffSBDDClient:
    """Client for running DiffSBDD on HPC cluster."""
    
    def __init__(self):
        """Initialize HPC DiffSBDD client."""
        self.ssh_client = get_hpc_client()
        self.scratch_base = "/scratch/rajagopalmohanraj.n/F.A.D.E/diffsbdd_inputs"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles"
        
    def generate_molecules(self, 
                         pdb_file: str,
                         pocket_residues: List[str],
                         n_samples: int = 100,
                         sanitize: bool = True,
                         job_id: str = None) -> Tuple[List[Dict[str, Any]], str]:
        """
        Generate molecules using DiffSBDD on HPC.
        
        Args:
            pdb_file: Path to local PDB file
            pocket_residues: List of residues defining the pocket (e.g., ["A:12", "A:61"])
            n_samples: Number of molecules to generate
            sanitize: Whether to sanitize generated molecules
            job_id: Optional job ID (will be generated if not provided)
            
        Returns:
            Tuple of (List of generated molecules, job_id used)
        """
        # Use provided job_id or generate new one
        if not job_id:
            job_id = str(uuid.uuid4())[:8]
        
        remote_dir = f"{self.scratch_base}/{job_id}"
        
        logger.info(f"Running DiffSBDD on HPC with job ID: {job_id}")
        
        try:
            # Connect to HPC
            if not self.ssh_client.connect():
                raise ConnectionError("Failed to connect to HPC cluster")
            
            # Create remote directory structure
            logger.info(f"Creating remote directories for job {job_id}")
            self.ssh_client.execute_command(f"mkdir -p {remote_dir}/inputs")
            self.ssh_client.execute_command(f"mkdir -p {remote_dir}/outputs")
            
            # Upload PDB file
            pdb_filename = Path(pdb_file).name
            remote_pdb = f"{remote_dir}/inputs/{pdb_filename}"
            
            logger.info(f"Uploading PDB to: {remote_pdb}")
            if not self.ssh_client.upload_file(pdb_file, remote_pdb):
                raise RuntimeError("Failed to upload PDB file")
            
            # Format residue list
            resi_list_str = " ".join(pocket_residues) if pocket_residues else ""
            
            # Build DiffSBDD command using sbatch_diffsbdd_generate
            output_file = f"generated_ligands_{job_id}.sdf"
            
            diffsbdd_cmd = f"""
cd {remote_dir}
module use {self.module_path}
module load diffsbdd

sbatch_diffsbdd_generate \\
    {job_id} \\
    checkpoints/crossdocked_fullatom_cond.ckpt \\
    --pdbfile {remote_pdb} \\
    --resi_list {resi_list_str} \\
    --outfile outputs/{output_file} \\
    --n_samples {n_samples} \\
    {"--sanitize" if sanitize else ""}
"""
            
            logger.info(f"Submitting DiffSBDD job to generate {n_samples} molecules...")
            logger.debug(f"Command: {diffsbdd_cmd}")
            
            # Submit the job
            stdout, stderr, exit_code = self.ssh_client.execute_command(
                diffsbdd_cmd,
                timeout=60  # 1 minute to submit
            )
            
            if exit_code != 0:
                raise RuntimeError(f"Failed to submit DiffSBDD job: {stderr}")
            
            # Extract SLURM job ID if available
            slurm_job_id = None
            job_id_match = re.search(r'Submitted batch job (\d+)', stdout)
            if job_id_match:
                slurm_job_id = job_id_match.group(1)
                logger.info(f"SLURM job submitted with ID: {slurm_job_id}")
            
            # Wait for job completion
            logger.info("Waiting for DiffSBDD job to complete...")
            self._wait_for_job_completion(remote_dir, output_file, slurm_job_id)
            
            # Download and parse the generated SDF file
            molecules = self._download_and_parse_sdf(remote_dir, job_id, output_file)
            
            if not molecules:
                raise ValueError("No molecules were generated")
            
            logger.info(f"Generated {len(molecules)} molecules")
            
            # Return molecules and job_id
            return molecules, job_id
            
        except Exception as e:
            logger.error(f"Error running DiffSBDD on HPC: {e}")
            raise
        
        finally:
            self.ssh_client.close()
    
    def _wait_for_job_completion(self, remote_dir: str, output_file: str, 
                                slurm_job_id: str = None, max_wait: int = 1200):
        """
        Wait for DiffSBDD job to complete.
        
        Args:
            remote_dir: Remote directory
            output_file: Expected output file name
            slurm_job_id: SLURM job ID if available
            max_wait: Maximum wait time in seconds (20 minutes default)
        """
        start_time = time.time()
        check_interval = 30  # Check every 30 seconds
        
        while time.time() - start_time < max_wait:
            # Check if output file exists
            check_cmd = f"test -f {remote_dir}/outputs/{output_file}"
            stdout, stderr, exit_code = self.ssh_client.execute_command(check_cmd)
            
            if exit_code == 0:
                logger.info("DiffSBDD output file found, job completed")
                return
            
            # If we have SLURM job ID, check job status
            if slurm_job_id:
                status_cmd = f"squeue -j {slurm_job_id} -h -o %T"
                stdout, stderr, exit_code = self.ssh_client.execute_command(status_cmd)
                
                if exit_code == 0 and stdout.strip():
                    job_status = stdout.strip()
                    logger.debug(f"Job {slurm_job_id} status: {job_status}")
                    
                    if job_status in ["COMPLETED", "FAILED", "CANCELLED"]:
                        if job_status != "COMPLETED":
                            raise RuntimeError(f"DiffSBDD job {slurm_job_id} {job_status}")
                        break
            
            logger.debug(f"Waiting for DiffSBDD job... ({int(time.time() - start_time)}s elapsed)")
            time.sleep(check_interval)
        
        # Final check for output file
        check_cmd = f"test -f {remote_dir}/outputs/{output_file}"
        stdout, stderr, exit_code = self.ssh_client.execute_command(check_cmd)
        
        if exit_code != 0:
            raise TimeoutError(f"DiffSBDD job did not complete within {max_wait} seconds")
    
    def _download_and_parse_sdf(self, remote_dir: str, job_id: str, 
                               output_file: str) -> List[Dict[str, Any]]:
        """
        Download and parse the generated SDF file.
        
        Args:
            remote_dir: Remote directory containing results
            job_id: Job ID
            output_file: Output file name
            
        Returns:
            List of molecule dictionaries
        """
        molecules = []
        
        # Download SDF file
        remote_sdf = f"{remote_dir}/outputs/{output_file}"
        local_sdf = config.DATA_DIR / "molecules" / f"diffsbdd_{job_id}.sdf"
        local_sdf.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Downloading generated molecules from {remote_sdf}")
        if not self.ssh_client.download_file(remote_sdf, str(local_sdf)):
            raise RuntimeError("Failed to download SDF file")
        
        # Parse SDF file to extract molecules
        try:
            # We'll need RDKit for proper parsing, but for now let's extract basic info
            molecules = self._parse_sdf_basic(local_sdf)
            
            # Also try to convert to SMILES using OpenBabel on HPC
            smiles_list = self._convert_to_smiles_hpc(remote_sdf, remote_dir)
            
            # Combine SDF data with SMILES
            for i, mol in enumerate(molecules):
                if i < len(smiles_list):
                    mol["smiles"] = smiles_list[i]
                mol["molecule_id"] = f"diffsbdd_{job_id}_{i+1}"
                mol["generation_method"] = "DiffSBDD"
                mol["source_file"] = str(local_sdf)
                mol["job_id"] = job_id
            
        except Exception as e:
            logger.error(f"Error parsing SDF file: {e}")
        
        return molecules
    
    def _parse_sdf_basic(self, sdf_file: Path) -> List[Dict[str, Any]]:
        """
        Basic SDF parser to extract molecule information.
        
        Args:
            sdf_file: Path to SDF file
            
        Returns:
            List of molecule dictionaries
        """
        molecules = []
        
        try:
            with open(sdf_file, 'r') as f:
                content = f.read()
            
            # Split by molecule delimiter
            mol_blocks = content.split("$$$$")
            
            for i, block in enumerate(mol_blocks):
                if not block.strip():
                    continue
                
                mol_info = {
                    "molecule_number": i + 1,
                    "sdf_block": block + "$$$$\n"
                }
                
                # Try to extract any properties from the SDF
                # Look for energy or score if present
                energy_match = re.search(r'<Energy>\s*([-\d.]+)', block)
                if energy_match:
                    mol_info["energy"] = float(energy_match.group(1))
                
                score_match = re.search(r'<Score>\s*([-\d.]+)', block)
                if score_match:
                    mol_info["score"] = float(score_match.group(1))
                
                molecules.append(mol_info)
                
        except Exception as e:
            logger.error(f"Error parsing SDF: {e}")
        
        return molecules
    
    def _convert_to_smiles_hpc(self, remote_sdf: str, remote_dir: str) -> List[str]:
        """
        Convert SDF to SMILES using OpenBabel on HPC.
        
        Args:
            remote_sdf: Remote SDF file path
            remote_dir: Remote directory
            
        Returns:
            List of SMILES strings
        """
        smiles_list = []
        
        try:
            # Use OpenBabel to convert SDF to SMILES
            smiles_file = f"{remote_dir}/outputs/molecules.smi"
            
            obabel_cmd = f"""
module load openbabel
obabel {remote_sdf} -O {smiles_file} -h
"""
            
            logger.info("Converting SDF to SMILES using OpenBabel...")
            stdout, stderr, exit_code = self.ssh_client.execute_command(obabel_cmd)
            
            if exit_code == 0:
                # Read SMILES file
                stdout, stderr, exit_code = self.ssh_client.execute_command(f"cat {smiles_file}")
                
                if exit_code == 0:
                    for line in stdout.strip().split('\n'):
                        if line:
                            # OpenBabel format: SMILES\tName
                            parts = line.split('\t')
                            if parts:
                                smiles_list.append(parts[0])
            else:
                logger.warning("OpenBabel conversion failed, trying alternative method")
                
        except Exception as e:
            logger.error(f"Error converting to SMILES: {e}")
        
        return smiles_list
    
    def _format_residues_for_diffsbdd(self, pocket: Dict[str, Any]) -> List[str]:
        """
        Format pocket residues for DiffSBDD input.
        
        Args:
            pocket: Pocket information dictionary
            
        Returns:
            List of residue strings in format "CHAIN:RESNUM"
        """
        formatted_residues = []
        
        if pocket.get("residues"):
            for res in pocket["residues"]:
                # Parse residue string (e.g., "GLY12" -> "A:12")
                match = re.match(r'([A-Z]{3})(\d+)', res)
                if match:
                    res_num = match.group(2)
                    # Assume chain A if not specified
                    formatted_residues.append(f"A:{res_num}")
        
        # For KRAS G12C, ensure we include key residues
        if not formatted_residues:
            # Default KRAS G12C binding site residues
            formatted_residues = ["A:12", "A:61", "A:116"]
        
        return formatted_residues


# Singleton instance
_hpc_diffsbdd_client = None

def get_hpc_diffsbdd_client() -> HPCDiffSBDDClient:
    """Get or create HPC DiffSBDD client singleton."""
    global _hpc_diffsbdd_client
    if _hpc_diffsbdd_client is None:
        _hpc_diffsbdd_client = HPCDiffSBDDClient()
    return _hpc_diffsbdd_client
