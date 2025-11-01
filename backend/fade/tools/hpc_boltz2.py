"""
HPC Boltz2 client for protein structure prediction.

This module provides functions to run Boltz2 on the HPC cluster
for de novo protein structure prediction.
"""

import os
import time
import logging
import tempfile
from typing import Dict, Any, Optional, Tuple
from pathlib import Path

import paramiko
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_boltz2")


class HPCBoltz2Client:
    """Client for running Boltz2 structure prediction on HPC."""
    
    def __init__(self, 
                 host: str = None,
                 username: str = None,
                 password: str = None,
                 key_path: str = None):
        """
        Initialize HPC Boltz2 client.
        
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
        self.remote_base = f"/scratch/{self.username}/F.A.D.E/boltz2_jobs"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles/"
        
        # Boltz2 settings
        self.boltz2_module = "boltz/2"
        self.max_wait_time = config.BOLTZ2_MAX_WAIT_TIME if hasattr(config, 'BOLTZ2_MAX_WAIT_TIME') else 1800
        self.poll_interval = config.BOLTZ2_POLL_INTERVAL if hasattr(config, 'BOLTZ2_POLL_INTERVAL') else 30
        self.partition = config.BOLTZ2_GPU_PARTITION if hasattr(config, 'BOLTZ2_GPU_PARTITION') else "gpu"
        self.slurm_account = config.BOLTZ2_SLURM_ACCOUNT if hasattr(config, 'BOLTZ2_SLURM_ACCOUNT') else None
        self.use_cpu = config.BOLTZ2_USE_CPU if hasattr(config, 'BOLTZ2_USE_CPU') else False
        self.BOLTZ_CACHE_DIR = config.BOLTZ_CACHE_DIR if hasattr(config, 'BOLTZ_CACHE_DIR') else "/scratch/rajagopalmohanraj.n/boltz_cache"
        
        self.ssh_client = None
        self.sftp_client = None
    
    def predict_structure(self, 
                         sequence: str,
                         job_name: str = None,
                         job_id: str = None,
                         use_msa: bool = True,
                         num_models: int = 1) -> Dict[str, Any]:
        """
        Predict protein structure using Boltz2.
        
        Args:
            sequence: Protein sequence (amino acids)
            job_name: Optional job name
            job_id: Optional job ID (reuse existing)
            use_msa: Whether to use MSA (Multiple Sequence Alignment)
            num_models: Number of models to generate
            
        Returns:
            Dictionary with job_id and status
        """
        # Generate job ID if not provided
        if not job_id:
            job_id = f"boltz2_{int(time.time())}"
        
        if not job_name:
            job_name = f"fade_boltz2_{job_id}"
        
        logger.info(f"Running Boltz2 prediction with job ID: {job_id}")
        
        try:
            # Connect to HPC
            self._connect()
            
            # Create remote directory
            remote_dir = f"{self.remote_base}/{job_id}"
            self._create_remote_directory(remote_dir)
            
            # Create and upload YAML file (Boltz2 likely expects YAML like AlphaFold3)
            yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {sequence}
"""
            yaml_path = f"{remote_dir}/input.yaml"
            self._upload_content(yaml_content, yaml_path)
            
            # Create and upload SLURM script
            script_content = self._create_boltz2_script(
                job_id=job_id,
                job_name=job_name,
                input_file="input.yaml",
                remote_dir=remote_dir,
                use_msa=use_msa,
                num_models=num_models
            )
            script_path = f"{remote_dir}/run_boltz2.sh"
            self._upload_content(script_content, script_path)
            
            # Submit job
            slurm_job_id = self._submit_job(script_path)
            
            logger.info(f"Boltz2 job submitted with SLURM ID: {slurm_job_id}")
            
            return {
                "job_id": job_id,
                "slurm_job_id": slurm_job_id,
                "status": "submitted",
                "remote_dir": remote_dir
            }
            
        except Exception as e:
            logger.error(f"Failed to submit Boltz2 job: {e}")
            raise
        finally:
            self._disconnect()
    
    def wait_for_structure(self, 
                          job_id: str,
                          max_wait: int = None,
                          poll_interval: int = None) -> Optional[str]:
        """
        Wait for Boltz2 job to complete and return structure.
        
        Args:
            job_id: Job ID to wait for
            max_wait: Maximum wait time in seconds
            poll_interval: Polling interval in seconds
            
        Returns:
            PDB content as string or None if failed
        """
        max_wait = max_wait or self.max_wait_time
        poll_interval = poll_interval or self.poll_interval
        
        logger.info(f"Waiting for Boltz2 job {job_id} (max {max_wait}s)")
        
        start_time = time.time()
        
        try:
            self._connect()
            remote_dir = f"{self.remote_base}/{job_id}"
            
            while time.time() - start_time < max_wait:
                # Check for different possible output file locations
                # Boltz2 creates files in output/boltz_results_input/predictions/input/
                possible_outputs = [
                    f"{remote_dir}/output/predicted_structure.cif",  # Our copied file
                    f"{remote_dir}/output/boltz_results_input/predictions/input/input_model_0.cif",  # Boltz2 default
                    f"{remote_dir}/output/protein.cif",
                    f"{remote_dir}/output/protein.pdb",
                ]
                
                for output_file in possible_outputs:
                    if self._file_exists(output_file):
                        logger.info(f"Boltz2 structure ready for job {job_id} at {output_file}")
                        
                        # Download structure file
                        structure_content = self._download_file(output_file)
                        
                        # Get confidence scores if available
                        confidence_files = [
                            f"{remote_dir}/output/confidence_scores.json",
                            f"{remote_dir}/output/boltz_results_input/predictions/input/confidence_input_model_0.json",
                        ]
                        confidence_data = None
                        for conf_file in confidence_files:
                            if self._file_exists(conf_file):
                                confidence_data = self._download_file(conf_file)
                                logger.info(f"Downloaded confidence scores from {conf_file}")
                                break
                        
                        return structure_content
                
                # Check for error
                error_file = f"{remote_dir}/boltz2.err"
                if self._file_exists(error_file):
                    error_content = self._download_file(error_file)
                    if error_content and len(error_content.strip()) > 0:
                        logger.error(f"Boltz2 error: {error_content[:500]}")
                        return None
                
                logger.debug(f"Waiting for Boltz2... ({int(time.time() - start_time)}s elapsed)")
                time.sleep(poll_interval)
            
            logger.error(f"Boltz2 job {job_id} timed out after {max_wait}s")
            return None
            
        except Exception as e:
            logger.error(f"Error waiting for Boltz2: {e}")
            return None
        finally:
            self._disconnect()
    
    def check_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check status of a Boltz2 job.
        
        Args:
            job_id: Job ID to check
            
        Returns:
            Status dictionary
        """
        try:
            self._connect()
            remote_dir = f"{self.remote_base}/{job_id}"
            
            status = {
                "job_id": job_id,
                "status": "unknown",
                "has_output": False,
                "has_error": False
            }
            
            # Check for output - Boltz2 creates in boltz_results_input directory
            possible_output_files = [
                f"{remote_dir}/output/predicted_structure.cif",
                f"{remote_dir}/output/boltz_results_input/predictions/input/input_model_0.cif",
            ]
            
            for output_file in possible_output_files:
                if self._file_exists(output_file):
                    status["status"] = "completed"
                    status["has_output"] = True
                    status["output_file"] = output_file
                    break
            
            # Check for error
            error_file = f"{remote_dir}/boltz2.err"
            if self._file_exists(error_file):
                error_content = self._download_file(error_file)
                if error_content and len(error_content.strip()) > 0:
                    status["status"] = "failed"
                    status["has_error"] = True
                    status["error_message"] = error_content[:500]
            
            # Check if still running
            log_file = f"{remote_dir}/boltz2.out"
            if self._file_exists(log_file) and status["status"] == "unknown":
                status["status"] = "running"
            
            return status
            
        except Exception as e:
            logger.error(f"Error checking job status: {e}")
            return {"job_id": job_id, "status": "error", "error": str(e)}
        finally:
            self._disconnect()
    
    def _create_boltz2_script(self,
                            job_id: str,
                            job_name: str,
                            input_file: str,
                            remote_dir: str,
                            use_msa: bool,
                            num_models: int) -> str:
        """
        Create SLURM script for Boltz2.
        
        Args:
            job_id: Job ID
            job_name: SLURM job name
            input_file: Input filename (YAML format)
            remote_dir: Remote working directory
            use_msa: Whether to use MSA
            num_models: Number of models to generate
            
        Returns:
            SLURM script content with Unix line endings
        """
        # Determine DEVICE based on the flag
        device_val = "cpu" if self.use_cpu else "gpu"

        # Determine MSA argument string based on the flag
        msa_arg_val = "--use_msa_server" if use_msa else "" # If false, becomes empty string

        # --- Start building the script ---
        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={remote_dir}/boltz2.out",
            f"#SBATCH --error={remote_dir}/boltz2.err",
            f"#SBATCH --partition={self.partition}",
        ]

        # Add account if specified
        if self.slurm_account:
            script_lines.append(f"#SBATCH --account={self.slurm_account}")

        # Add resource requirements based on CPU/GPU mode
        if self.use_cpu:
            script_lines.extend([
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks=4",
                "#SBATCH --cpus-per-task=4",
                "#SBATCH --mem=32G", # Adjust memory as needed for CPU runs
            ])
        else:
            script_lines.extend([
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks-per-node=4", # Adjust if needed
                "#SBATCH --mem=32G", # Adjust memory based on model size/GPU needs
                "#SBATCH --gres=gpu:1", # Request 1 GPU
            ])

        script_lines.extend([
            "#SBATCH --time=02:00:00", # Adjust time as needed
            "",
            "# Load modules",
            f"module use {self.module_path}",
            f"module load {self.boltz2_module}",
            "",
            "# Set working directory",
            f"cd {remote_dir}",
            "",
            "# Create output directory",
            "mkdir -p output",
            "",
            "# Debug: Check input file content",
            "echo \"Input file content:\"",
            f"head -10 {input_file}",
            "echo \"\"",
            "# Run Boltz2",
            "echo \"Starting Boltz2 prediction...\"",
            f"echo \"Input: {input_file}\"",
            "echo \"Output directory: output/\"",
            f"echo \"Number of models: {num_models}\"",
            f"echo \"Use MSA: {use_msa}\"",
            f"echo \"Device: {device_val}\"",
            f"echo \"Cache: {self.BOLTZ_CACHE_DIR}\"",
            "",
            "# Construct the boltz command carefully",
            "boltz predict \\", # Start command
            f"    --cache \"{self.BOLTZ_CACHE_DIR}\" \\", # Cache dir
            f"    --out_dir output/ \\", # Output relative dir
            f"    --diffusion_samples {num_models} \\", # Number of samples
        ])

        # Only add the MSA flag if msa_arg_val is not empty
        if msa_arg_val:
            script_lines.append(f"    {msa_arg_val} \\")

        script_lines.extend([
            f"    --accelerator {device_val} \\", # Device type
            f"    --seed 42 \\", # Seed
            f"    \"{input_file}\"", # Input YAML file (LAST positional argument)
            "",
            "# Capture exit code",
            "EXIT_CODE=$?",
            "echo \"Boltz2 exit code: $EXIT_CODE\"",
            "",
            "# Check for output files regardless of exit code",
            "# The error might be non-fatal since we see 100% progress",
            "echo \"Checking for all generated files...\"",
            "find . -type f -newer $0 2>/dev/null | head -20",
            "echo \"\"",
            "",
            "# Check if successful",
            "if [ $EXIT_CODE -eq 0 ]; then",
            "    echo \"Boltz2 completed successfully\"",
            "    ",
            "    # Debug: List output directory contents",
            "    echo \"Contents of output directory:\"",
            "    ls -la output/ 2>/dev/null || echo \"No output directory\"",
            "    echo \"\"",
            "    ",
            "    # Debug: Find any generated files",
            "    echo \"Looking for any structure files:\"",
            "    find . -name '*.cif' -o -name '*.pdb' 2>/dev/null | head -10",
            "    echo \"\"",
            "    ",
            "    # Check for Boltz2 output structure",
            "    # Boltz2 creates output in boltz_results_<input_name> directory",
            "    if [ -d \"output/boltz_results_input\" ]; then",
            "        echo \"Found Boltz2 results directory\"",
            "        # Look for the CIF file in the predictions subdirectory",
            "        CIF_FILE=\"output/boltz_results_input/predictions/input/input_model_0.cif\"",
            "        if [ -f \"$CIF_FILE\" ]; then",
            "            cp \"$CIF_FILE\" output/predicted_structure.cif",
            "            echo \"Structure copied from $CIF_FILE to output/predicted_structure.cif\"",
            "        else",
            "            # Try to find any CIF file in the boltz results",
            "            STRUCTURE_FILE=$(find output/boltz_results_input -name '*.cif' | head -1)",
            "            if [ -f \"$STRUCTURE_FILE\" ]; then",
            "                cp \"$STRUCTURE_FILE\" output/predicted_structure.cif",
            "                echo \"Structure copied from $STRUCTURE_FILE\"",
            "            fi",
            "        fi",
            "        ",
            "        # Also copy confidence data if available",
            "        CONF_FILE=\"output/boltz_results_input/predictions/input/confidence_input_model_0.json\"",
            "        if [ -f \"$CONF_FILE\" ]; then",
            "            cp \"$CONF_FILE\" output/confidence_scores.json",
            "            echo \"Confidence scores copied to output/confidence_scores.json\"",
            "        fi",
            "    else",
            "        echo \"Warning: Expected Boltz2 output directory not found\"",
            "    fi",
            "else",
            "    echo \"Boltz2 failed with exit code $EXIT_CODE\"",
            "fi",
            "",
            "echo \"Job completed at $(date)\"",
            "",
            "# Ensure Slurm gets the correct exit code",
            "exit $EXIT_CODE",
        ])

        # Join lines with Unix line endings
        script = "\n".join(script_lines)
        return script
    
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
        
        # Check for line ending errors specifically
        if "DOS line breaks" in error:
            logger.error("Script has DOS line breaks. Attempting to fix...")
            # Try to fix and resubmit
            fix_cmd = f"tr -d '\r' < {script_path} > {script_path}.fixed && mv {script_path}.fixed {script_path}"
            stdin, stdout, stderr = self.ssh_client.exec_command(fix_cmd)
            stdout.read()
            
            # Retry submission
            stdin, stdout, stderr = self.ssh_client.exec_command(submit_cmd)
            output = stdout.read().decode()
            error = stderr.read().decode()
        
        if error and "Submitted batch job" not in error:
            logger.warning(f"SLURM submission warning: {error}")
        
        # Parse SLURM job ID from output or error (sometimes it's in stderr)
        import re
        combined_output = output + error
        match = re.search(r'Submitted batch job (\d+)', combined_output)
        if match:
            logger.info(f"SLURM job submitted successfully: {match.group(1)}")
            return match.group(1)
        else:
            raise RuntimeError(f"Failed to submit job. Output: {output} Error: {error}")


# Singleton instance
_hpc_boltz2_client = None

def get_hpc_boltz2_client() -> HPCBoltz2Client:
    """Get or create HPC Boltz2 client singleton."""
    global _hpc_boltz2_client
    if _hpc_boltz2_client is None:
        _hpc_boltz2_client = HPCBoltz2Client()
    return _hpc_boltz2_client
