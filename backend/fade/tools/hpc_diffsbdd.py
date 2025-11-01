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
            
            # Build DiffSBDD command - create custom sbatch script
            output_file = f"generated_ligands_{job_id}.sdf"
            
            # Create a custom SLURM batch script with proper environment activation
            batch_script = self._create_batch_script(
                job_id=job_id,
                remote_dir=remote_dir,
                remote_pdb=remote_pdb,
                pocket_residues=pocket_residues,
                output_file=output_file,
                n_samples=n_samples,
                sanitize=sanitize
            )
            
            # Upload the batch script
            batch_script_path = f"{remote_dir}/run_diffsbdd.sh"
            self._upload_batch_script(batch_script, batch_script_path)
            
            # Submit the job
            diffsbdd_cmd = f"cd {remote_dir} && sbatch run_diffsbdd.sh"
            
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
                # First try squeue (for running/pending jobs)
                status_cmd = f"squeue -j {slurm_job_id} -h -o %T"
                stdout, stderr, exit_code = self.ssh_client.execute_command(status_cmd)
                
                if exit_code == 0 and stdout.strip():
                    job_status = stdout.strip()
                    logger.info(f"Job {slurm_job_id} status: {job_status}")
                    
                    if job_status in ["PENDING", "RUNNING"]:
                        # Job is still active, continue waiting
                        logger.debug(f"Waiting for DiffSBDD job... ({int(time.time() - start_time)}s elapsed)")
                        time.sleep(check_interval)
                        continue
                    elif job_status in ["COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"]:
                        logger.warning(f"Job {slurm_job_id} finished with status: {job_status}")
                        break
                else:
                    # Job not in squeue, check sacct for completed jobs
                    logger.debug(f"Job {slurm_job_id} not in queue, checking sacct...")
                    sacct_cmd = f"sacct -j {slurm_job_id} --format=JobID,State,ExitCode -n -P"
                    stdout, stderr, exit_code = self.ssh_client.execute_command(sacct_cmd)
                    
                    if exit_code == 0 and stdout.strip():
                        # Parse sacct output (format: JobID|State|ExitCode)
                        lines = stdout.strip().split('\n')
                        for line in lines:
                            if '|' in line:
                                parts = line.split('|')
                                if len(parts) >= 2:
                                    job_id_part = parts[0]
                                    job_state = parts[1]
                                    
                                    # Only check the main job (not .batch or .extern)
                                    if job_id_part == slurm_job_id:
                                        logger.info(f"Job {slurm_job_id} completed with state: {job_state}")
                                        
                                        if job_state == "COMPLETED":
                                            break
                                        else:
                                            # Job failed, get error logs
                                            self._print_job_error_logs(remote_dir, slurm_job_id)
                                            raise RuntimeError(f"DiffSBDD job {slurm_job_id} {job_state}")
                        break
                    else:
                        # Can't find job info, continue waiting
                        logger.debug(f"No job info found, continuing to wait... ({int(time.time() - start_time)}s elapsed)")
            
            time.sleep(check_interval)
        
        # Final check for output file
        check_cmd = f"test -f {remote_dir}/outputs/{output_file}"
        stdout, stderr, exit_code = self.ssh_client.execute_command(check_cmd)
        
        if exit_code != 0:
            # Check error logs before failing
            if slurm_job_id:
                self._print_job_error_logs(remote_dir, slurm_job_id)
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
    
    def _create_batch_script(self, job_id: str, remote_dir: str, remote_pdb: str,
                           pocket_residues: List[str], output_file: str,
                           n_samples: int, sanitize: bool) -> str:
        """
        Create a SLURM batch script with proper environment activation.
        
        Args:
            job_id: Job identifier
            remote_dir: Remote working directory
            remote_pdb: Path to PDB file on HPC
            pocket_residues: List of pocket residues (can be empty)
            output_file: Output filename
            n_samples: Number of molecules to generate
            sanitize: Whether to sanitize molecules
            
        Returns:
            Batch script content as string
        """
        sanitize_flag = "--sanitize" if sanitize else ""
        
        # Path to DiffSBDD installation and checkpoint
        diffsbdd_path = "/projects/SimBioSys/share/software/DiffSBDD"
        checkpoint_path = f"{diffsbdd_path}/checkpoints/crossdocked_fullatom_cond.ckpt"
        env_path = f"{diffsbdd_path}/env/"
        
        # Build residue list argument only if residues are provided
        resi_list_arg = ""
        if pocket_residues:
            resi_list_str = " ".join(pocket_residues)
            resi_list_arg = f"--resi_list {resi_list_str} \\"
        
        batch_script = f'''#!/bin/bash
#SBATCH --job-name=diffsbdd_{job_id}
#SBATCH --output={remote_dir}/slurm-%j.out
#SBATCH --error={remote_dir}/slurm-%j.err
#SBATCH --partition=gpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:1

echo "======================================================"
echo "Starting DiffSBDD Job"
echo "Job ID: $SLURM_JOB_ID"
echo "Unique ID: {job_id}"
echo "Output Directory: {remote_dir}"
echo "Running on: $(hostname)"
echo "======================================================"

# Load module and activate environment
module use {self.module_path}
module load diffsbdd

# Activate micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate {env_path}

echo "Environment activated"
echo "Python: $(which python)"
echo "PyTorch version: $(python -c 'import torch; print(torch.__version__)')"

# Change to working directory
cd {remote_dir}

# Run DiffSBDD generation
# Note: checkpoint is a positional argument and must come first
echo "Running DiffSBDD..."
python {diffsbdd_path}/generate_ligands.py \\
    {checkpoint_path} \\
    --pdbfile {remote_pdb} \\
    {resi_list_arg}
    --outfile outputs/{output_file} \\
    --n_samples {n_samples} \\
    {sanitize_flag}

EXIT_CODE=$?

echo "======================================================"
echo "DiffSBDD Job Finished with Exit Code: $EXIT_CODE"
echo "======================================================"

exit $EXIT_CODE
'''
        
        return batch_script
    
    def _upload_batch_script(self, script_content: str, remote_path: str):
        """
        Upload batch script to HPC.
        
        Args:
            script_content: Content of the batch script
            remote_path: Destination path on HPC
        """
        import tempfile
        
        # Convert Windows line endings (\r\n) to Unix line endings (\n)
        script_content = script_content.replace('\r\n', '\n')
        
        # Create temporary file with script content
        # Use newline='' to prevent Python from converting line endings
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False, newline='') as f:
            f.write(script_content)
            temp_path = f.name
        
        try:
            # Upload the script
            if not self.ssh_client.upload_file(temp_path, remote_path):
                raise RuntimeError("Failed to upload batch script")
            
            # Make it executable
            self.ssh_client.execute_command(f"chmod +x {remote_path}")
            
            logger.debug(f"Uploaded batch script to {remote_path}")
            
        finally:
            # Clean up temp file
            import os
            try:
                os.unlink(temp_path)
            except:
                pass
    
    def _print_job_error_logs(self, remote_dir: str, slurm_job_id: str):
        """
        Print SLURM job error and output logs for debugging.
        
        Args:
            remote_dir: Remote directory
            slurm_job_id: SLURM job ID
        """
        # Check for slurm output file
        slurm_out = f"{remote_dir}/slurm-{slurm_job_id}.out"
        stdout, stderr, exit_code = self.ssh_client.execute_command(f"cat {slurm_out}")
        
        if exit_code == 0 and stdout:
            logger.error(f"SLURM output for job {slurm_job_id}:")
            logger.error(stdout[-2000:])  # Last 2000 chars
        else:
            logger.warning(f"No SLURM output file found at {slurm_out}")
        
        # Check for error file
        slurm_err = f"{remote_dir}/slurm-{slurm_job_id}.err"
        stdout, stderr, exit_code = self.ssh_client.execute_command(f"cat {slurm_err}")
        
        if exit_code == 0 and stdout:
            logger.error(f"SLURM error for job {slurm_job_id}:")
            logger.error(stdout[-2000:])
    
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
