"""
HPC P2Rank client for protein pocket prediction.

P2Rank is a machine learning-based tool for predicting ligand binding sites
from protein structure. It uses a random forest classifier trained on
protein-ligand complexes.

This module provides functions to run P2Rank on the HPC cluster
for pocket detection and ranking.
"""

import os
import time
import logging
import tempfile
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

import paramiko
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_p2rank")


class HPCP2RankClient:
    """Client for running P2Rank pocket prediction on HPC."""
    
    def __init__(self, 
                 host: str = None,
                 username: str = None,
                 password: str = None,
                 key_path: str = None):
        """
        Initialize HPC P2Rank client.
        
        Args:
            host: HPC host address
            username: SSH username
            password: SSH password (if not using key)
            key_path: Path to SSH private key
        """
        self.host = host or os.getenv("HPC_HOST", "explorer.northeastern.edu")
        self.username = username or os.getenv("HPC_USER", "rajagopalmohanraj.n")
        self.password = password or os.getenv("HPC_PASSWORD")
        self.key_path = key_path or os.getenv("HPC_KEY_PATH")
        
        # HPC paths
        self.remote_base = f"/scratch/{self.username}/F.A.D.E/p2rank_jobs"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles/"
        
        # P2Rank settings
        self.p2rank_module = "p2rank/2.5.1"
        # TODO: Add more configuration options as we discover them
        
        self.ssh_client = None
        self.sftp_client = None
        
        logger.info("Initialized P2Rank HPC client")
    
    def predict_pockets(self, 
                       pdb_content: str,
                       job_name: str = None,
                       job_id: str = None) -> Dict[str, Any]:
        """
        Predict pockets in protein structure using P2Rank.
        
        Args:
            pdb_content: PDB file content as string
            job_name: Optional job name
            job_id: Optional job ID (reuse existing)
            
        Returns:
            Dictionary with job_id, status, and results
        """
        # Generate job ID if not provided
        if not job_id:
            job_id = f"p2rank_{int(time.time())}"
        
        if not job_name:
            job_name = f"fade_p2rank_{job_id}"
        
        logger.info(f"Running P2Rank prediction with job ID: {job_id}")
        
        try:
            # Connect to HPC
            self._connect()
            
            # Create remote directory
            remote_dir = f"{self.remote_base}/{job_id}"
            self._create_remote_directory(remote_dir)
            
            # Upload PDB file
            pdb_path = f"{remote_dir}/input.pdb"
            self._upload_content(pdb_content, pdb_path)
            
            # Create and upload SLURM script
            script_content = self._create_p2rank_script(
                job_id=job_id,
                job_name=job_name,
                input_file="input.pdb",
                remote_dir=remote_dir
            )
            script_path = f"{remote_dir}/run_p2rank.sh"
            self._upload_content(script_content, script_path)
            
            # Submit job
            slurm_job_id = self._submit_job(script_path)
            
            logger.info(f"P2Rank job submitted with SLURM ID: {slurm_job_id}")
            
            return {
                "job_id": job_id,
                "slurm_job_id": slurm_job_id,
                "status": "submitted",
                "remote_dir": remote_dir
            }
            
        except Exception as e:
            logger.error(f"Failed to submit P2Rank job: {e}")
            raise
        finally:
            self._disconnect()
    
    def _connect(self):
        """Establish SSH connection to HPC."""
        if self.ssh_client and self.ssh_client.get_transport() and self.ssh_client.get_transport().is_active():
            return
        
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        connect_kwargs = {
            'hostname': self.host,
            'username': self.username,
        }
        
        if self.key_path and os.path.exists(self.key_path):
            connect_kwargs['key_filename'] = self.key_path
            logger.info(f"Connecting to {self.host} with key file")
        elif self.password:
            connect_kwargs['password'] = self.password
            logger.info(f"Connecting to {self.host} with password")
        else:
            raise ValueError("No authentication method available (no key or password)")
        
        self.ssh_client.connect(**connect_kwargs)
        self.sftp_client = self.ssh_client.open_sftp()
        logger.info(f"Successfully connected to {self.host}")
    
    def _create_p2rank_script(self,
                             job_id: str,
                             job_name: str,
                             input_file: str,
                             remote_dir: str) -> str:
        """
        Create SLURM script for P2Rank.
        
        Args:
            job_id: Job ID
            job_name: SLURM job name
            input_file: Input PDB filename
            remote_dir: Remote working directory
            
        Returns:
            SLURM script content with Unix line endings
        """
        # Build the script line by line
        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={remote_dir}/p2rank.out",
            f"#SBATCH --error={remote_dir}/p2rank.err",
            "#SBATCH --partition=short",  # P2Rank doesn't need GPU
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=4",
            "#SBATCH --mem=16G",
            "#SBATCH --time=00:30:00",  # P2Rank is relatively fast
            "",
            "# Load modules",
            f"module use {self.module_path}",
            f"module load {self.p2rank_module}",
            "",
            "# Set working directory",
            f"cd {remote_dir}",
            "",
            "# Create output directory",
            "mkdir -p output",
            "",
            "# Debug: Check input file",
            "echo \"Input PDB file:\"",
            f"ls -la {input_file}",
            "echo \"\"",
            "",
            "# Run P2Rank",
            "echo \"Starting P2Rank prediction...\"",
            f"echo \"Input: {input_file}\"",
            "echo \"Output directory: output/\"",
            "",
            "# P2Rank command with offline mode to avoid download attempts",
            "# Note: P2Rank outputs predictions.csv with pocket information",
            "# The module might provide 'prank' instead of 'p2rank'",
            "echo \"Checking for P2Rank executable...\"",
            "which prank 2>/dev/null && echo \"Found: prank\" || echo \"prank not found\"",
            "which p2rank 2>/dev/null && echo \"Found: p2rank\" || echo \"p2rank not found\"",
            "which P2Rank 2>/dev/null && echo \"Found: P2Rank\" || echo \"P2Rank not found\"",
            "echo \"\"",
            "",
            "# Set P2Rank to skip downloading (work offline)",
            "export JAVA_OPTS=\"-Djava.net.useSystemProxies=true -Dhttp.proxyHost=none\"",
            "",
            "# Try different possible command names",
            "if command -v prank &> /dev/null; then",
            "    echo \"Using 'prank' command\"",
            f"    prank predict -f {input_file} -o output/ 2>&1 | grep -v \"DownloadChemCompProvider\"",
            "elif command -v P2Rank &> /dev/null; then",
            "    echo \"Using 'P2Rank' command\"",
            f"    P2Rank predict -f {input_file} -o output/ 2>&1 | grep -v \"DownloadChemCompProvider\"",
            "elif [ -f \"$P2RANK_HOME/prank\" ]; then",
            "    echo \"Using P2RANK_HOME/prank\"",
            f"    $P2RANK_HOME/prank predict -f {input_file} -o output/ 2>&1 | grep -v \"DownloadChemCompProvider\"",
            "else",
            "    echo \"ERROR: Could not find P2Rank executable!\"",
            "    echo \"Checking environment variables...\"",
            "    echo \"P2RANK_HOME: $P2RANK_HOME\"",
            "    echo \"PATH: $PATH\"",
            "    echo \"Contents of P2RANK_HOME (if set):\"",
            "    [ -n \"$P2RANK_HOME\" ] && ls -la \"$P2RANK_HOME\" 2>/dev/null",
            "    exit 1",
            "fi",
            "",
            "# Capture exit code",
            "EXIT_CODE=$?",
            "echo \"P2Rank exit code: $EXIT_CODE\"",
            "",
            "# Check for output",
            "if [ $EXIT_CODE -eq 0 ]; then",
            "    echo \"P2Rank completed successfully\"",
            "    ",
            "    # List output files",
            "    echo \"Output files:\"",
            "    ls -la output/",
            "    ",
            "    # Check for predictions file with different possible names",
            "    if [ -f \"output/input.pdb_predictions.csv\" ]; then",
            "        echo \"Predictions CSV found: input.pdb_predictions.csv\"",
            "        head -5 output/input.pdb_predictions.csv",
            "    elif [ -f \"output/input_predictions.csv\" ]; then",
            "        echo \"Predictions CSV found: input_predictions.csv\"",
            "        head -5 output/input_predictions.csv",  
            "    elif [ -f \"output/predictions.csv\" ]; then",
            "        echo \"Predictions CSV found: predictions.csv\"",
            "        head -5 output/predictions.csv",
            "    else",
            "        echo \"Looking for any CSV files in output...\"",
            "        find output/ -name \"*.csv\" 2>/dev/null",
            "    fi",
            "else",
            "    echo \"P2Rank failed with exit code $EXIT_CODE\"",
            "fi",
            "",
            "echo \"Job completed at $(date)\"",
            "",
            "exit $EXIT_CODE",
        ]
        
        # Join lines with Unix line endings
        script = "\n".join(script_lines)
        return script
    
    def _disconnect(self):
        """Close SSH connection."""
        if self.sftp_client:
            self.sftp_client.close()
            self.sftp_client = None
        if self.ssh_client:
            self.ssh_client.close()
            self.ssh_client = None
        logger.debug("SSH connection closed")
    
    def _create_remote_directory(self, path: str):
        """Create directory on remote server."""
        stdin, stdout, stderr = self.ssh_client.exec_command(f"mkdir -p {path}")
        stdout.read()
        logger.debug(f"Created remote directory: {path}")
    
    def _upload_content(self, content: str, remote_path: str):
        """Upload content as file to remote server."""
        # Convert Windows line endings to Unix line endings
        content = content.replace('\r\n', '\n')
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt', newline='\n') as f:
            f.write(content)
            temp_path = f.name
        
        try:
            self.sftp_client.put(temp_path, remote_path)
            logger.debug(f"Uploaded content to {remote_path}")
        finally:
            os.unlink(temp_path)
    
    def _download_file(self, remote_path: str) -> str:
        """Download file from remote server."""
        with tempfile.NamedTemporaryFile(mode='r', delete=False) as f:
            temp_path = f.name
        
        try:
            self.sftp_client.get(remote_path, temp_path)
            with open(temp_path, 'r') as f:
                content = f.read()
            return content
        finally:
            os.unlink(temp_path)
    
    def _file_exists(self, remote_path: str) -> bool:
        """Check if file exists on remote server."""
        try:
            self.sftp_client.stat(remote_path)
            return True
        except FileNotFoundError:
            return False
    
    def _submit_job(self, script_path: str) -> str:
        """Submit SLURM job and return job ID."""
        # First, convert the script file to Unix format on the remote server
        convert_cmd = f"dos2unix {script_path} 2>/dev/null || sed -i 's/\r$//' {script_path}"
        stdin, stdout, stderr = self.ssh_client.exec_command(convert_cmd)
        stdout.read()  # Wait for completion
        
        # Now submit the job
        submit_cmd = f"cd {os.path.dirname(script_path)} && sbatch {os.path.basename(script_path)}"
        stdin, stdout, stderr = self.ssh_client.exec_command(submit_cmd)
        output = stdout.read().decode()
        error = stderr.read().decode()
        
        if error and "Submitted batch job" not in error:
            logger.warning(f"SLURM submission warning: {error}")
        
        # Parse SLURM job ID from output or error
        import re
        combined_output = output + error
        match = re.search(r'Submitted batch job (\d+)', combined_output)
        if match:
            logger.info(f"SLURM job submitted successfully: {match.group(1)}")
            return match.group(1)
        else:
            raise RuntimeError(f"Failed to submit job. Output: {output} Error: {error}")
    
    def wait_for_results(self,
                        job_id: str,
                        max_wait: int = 1800,
                        poll_interval: int = 30) -> Optional[List[Dict[str, Any]]]:
        """
        Wait for P2Rank job to complete and parse results.
        
        Args:
            job_id: Job ID to wait for
            max_wait: Maximum wait time in seconds (default 30 minutes)
            poll_interval: Polling interval in seconds
            
        Returns:
            List of predicted pockets or None if failed
        """
        logger.info(f"Waiting for P2Rank job {job_id} (max {max_wait}s)")
        
        start_time = time.time()
        
        try:
            self._connect()
            remote_dir = f"{self.remote_base}/{job_id}"
            
            while time.time() - start_time < max_wait:
                # Check for different possible output file names
                # P2Rank might use different naming conventions
                possible_csv_files = [
                    f"{remote_dir}/output/input.pdb_predictions.csv",
                    f"{remote_dir}/output/input_predictions.csv",
                    f"{remote_dir}/output/predictions.csv",
                    f"{remote_dir}/output/prank_predictions.csv",
                ]
                
                csv_content = None
                csv_file_found = None
                
                for csv_file in possible_csv_files:
                    if self._file_exists(csv_file):
                        logger.info(f"P2Rank results found at {csv_file}")
                        csv_content = self._download_file(csv_file)
                        csv_file_found = csv_file
                        break
                
                # Also check if output directory exists and list contents
                if not csv_content:
                    output_dir = f"{remote_dir}/output"
                    try:
                        stdin, stdout, stderr = self.ssh_client.exec_command(f"ls -la {output_dir}/ 2>/dev/null")
                        output_contents = stdout.read().decode()
                        if output_contents and "predictions" in output_contents.lower():
                            logger.debug(f"Output directory contents:\n{output_contents}")
                            # Try to find any CSV file
                            stdin, stdout, stderr = self.ssh_client.exec_command(f"find {output_dir} -name '*.csv' 2>/dev/null")
                            csv_files = stdout.read().decode().strip().split('\n')
                            if csv_files and csv_files[0]:
                                csv_file = csv_files[0]
                                logger.info(f"Found CSV file: {csv_file}")
                                csv_content = self._download_file(csv_file)
                                csv_file_found = csv_file
                    except Exception as e:
                        logger.debug(f"Could not list output directory: {e}")
                
                if csv_content:
                    logger.info(f"P2Rank results ready for job {job_id}")
                    pockets = self._parse_p2rank_csv(csv_content)
                    return pockets
                
                # Check for error
                error_file = f"{remote_dir}/p2rank.err"
                if self._file_exists(error_file):
                    error_content = self._download_file(error_file)
                    if error_content and len(error_content.strip()) > 0:
                        logger.error(f"P2Rank error: {error_content[:500]}")
                        return None
                
                logger.debug(f"Waiting for P2Rank... ({int(time.time() - start_time)}s elapsed)")
                time.sleep(poll_interval)
            
            logger.error(f"P2Rank job {job_id} timed out after {max_wait}s")
            return None
            
        except Exception as e:
            logger.error(f"Error waiting for P2Rank results: {e}")
            return None
        finally:
            self._disconnect()
    
    def _parse_p2rank_csv(self, csv_content: str) -> List[Dict[str, Any]]:
        """
        Parse P2Rank CSV output.
        
        Expected CSV format:
        rank,name,score,probability,sas_points,surf_atoms,center_x,center_y,center_z,residue_ids,residue_names,residue_chains,surf_atom_ids
        
        Args:
            csv_content: CSV file content as string
            
        Returns:
            List of pocket dictionaries
        """
        pockets = []
        lines = csv_content.strip().split('\n')
        
        if len(lines) < 2:
            logger.warning("P2Rank CSV has no data")
            return pockets
        
        # Skip header
        header = lines[0].split(',')
        
        for line in lines[1:]:
            if not line.strip():
                continue
            
            parts = line.split(',')
            if len(parts) < 9:
                logger.warning(f"Skipping malformed P2Rank line: {line}")
                continue
            
            try:
                pocket = {
                    "rank": int(parts[0]),
                    "name": parts[1],
                    "score": float(parts[2]),
                    "probability": float(parts[3]),
                    "sas_points": int(parts[4]) if parts[4] else 0,
                    "surf_atoms": int(parts[5]) if parts[5] else 0,
                    "center": [
                        float(parts[6]),  # center_x
                        float(parts[7]),  # center_y
                        float(parts[8])   # center_z
                    ],
                    "residues": []
                }
                
                # Parse residues if present
                if len(parts) > 10 and parts[10]:
                    residue_names = parts[10].strip().split()
                    pocket["residues"] = residue_names
                
                pockets.append(pocket)
                
            except (ValueError, IndexError) as e:
                logger.warning(f"Error parsing P2Rank line: {e}")
                continue
        
        logger.info(f"Parsed {len(pockets)} pockets from P2Rank")
        return pockets


# Singleton instance
_hpc_p2rank_client = None

def get_hpc_p2rank_client() -> HPCP2RankClient:
    """Get or create HPC P2Rank client singleton."""
    global _hpc_p2rank_client
    if _hpc_p2rank_client is None:
        _hpc_p2rank_client = HPCP2RankClient()
    return _hpc_p2rank_client
