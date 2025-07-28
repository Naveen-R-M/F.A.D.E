"""
SLURM Client for F.A.D.E

This module provides an interface to interact with the SLURM job scheduler
for submitting and monitoring jobs on HPC clusters.
"""

import os
import subprocess
import re
import time
from typing import Any, Dict, List, Optional, Tuple, Union

from utils.logging import get_logger


class SlurmClient:
    """
    Client for interacting with the SLURM job scheduler.
    """
    
    def __init__(self) -> None:
        """Initialize the SLURM client."""
        self.logger = get_logger("fade.utils.slurm_client")
    
    def submit_job(self, job_script: str) -> str:
        """
        Submit a job to SLURM.
        
        Args:
            job_script: Path to the job script.
            
        Returns:
            Job ID as a string.
            
        Raises:
            ValueError: If the job submission fails.
        """
        # Check if job script exists
        if not os.path.exists(job_script):
            raise ValueError(f"Job script not found: {job_script}")
            
        # Submit job using sbatch
        try:
            result = subprocess.run(
                ["sbatch", job_script],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse job ID from output
            # Expected output format: "Submitted batch job 12345"
            match = re.search(r"Submitted batch job (\d+)", result.stdout)
            
            if match:
                job_id = match.group(1)
                self.logger.info(f"Job submitted with ID: {job_id}")
                return job_id
            else:
                raise ValueError(f"Failed to parse job ID from output: {result.stdout}")
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Job submission failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise ValueError(f"Failed to submit job: {e}")
    
    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a SLURM job.
        
        Args:
            job_id: Job ID to cancel.
            
        Returns:
            True if the job was canceled successfully, False otherwise.
        """
        try:
            subprocess.run(
                ["scancel", job_id],
                capture_output=True,
                text=True,
                check=True
            )
            
            self.logger.info(f"Job {job_id} canceled")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to cancel job {job_id}: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return False
    
    def get_job_status(self, job_id: str) -> str:
        """
        Get the status of a SLURM job.
        
        Args:
            job_id: Job ID to check.
            
        Returns:
            Job status as a string: PENDING, RUNNING, COMPLETED, FAILED, etc.
            Returns "UNKNOWN" if the job status cannot be determined.
            
        Note:
            This function uses the sacct command to query job status.
        """
        try:
            # Use sacct to get job status
            # Format: JobID,State
            result = subprocess.run(
                [
                    "sacct",
                    "-j", job_id,
                    "--format=JobID,State",
                    "--parsable2",  # Use | as delimiter
                    "--noheader"    # Skip header row
                ],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse output
            lines = result.stdout.strip().split("\n")
            
            if not lines:
                self.logger.warning(f"No output from sacct for job {job_id}")
                return "UNKNOWN"
                
            # Get status from the first line
            # Format: JobID|State
            parts = lines[0].split("|")
            
            if len(parts) < 2:
                self.logger.warning(f"Unexpected sacct output format for job {job_id}: {lines[0]}")
                return "UNKNOWN"
                
            status = parts[1].strip()
            self.logger.info(f"Job {job_id} status: {status}")
            
            return status
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to get status for job {job_id}: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return "UNKNOWN"
    
    def monitor_job(self, job_id: str, check_interval: int = 60, timeout: Optional[int] = None) -> str:
        """
        Monitor a SLURM job until it completes or fails.
        
        Args:
            job_id: Job ID to monitor.
            check_interval: Interval in seconds between status checks.
            timeout: Optional timeout in seconds. If None, wait indefinitely.
            
        Returns:
            Final job status as a string.
        """
        self.logger.info(f"Monitoring job {job_id}")
        
        start_time = time.time()
        
        while True:
            # Check if timeout has been reached
            if timeout and (time.time() - start_time) > timeout:
                self.logger.warning(f"Timeout reached for job {job_id}")
                return "TIMEOUT"
                
            # Get job status
            status = self.get_job_status(job_id)
            
            # Check if job has completed or failed
            if status in ["COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"]:
                self.logger.info(f"Job {job_id} finished with status: {status}")
                return status
                
            # Wait before next check
            time.sleep(check_interval)
    
    def get_job_output(self, job_id: str, output_file: Optional[str] = None) -> str:
        """
        Get the output of a SLURM job.
        
        Args:
            job_id: Job ID to get output for.
            output_file: Optional path to the output file. If not provided,
                        the function will try to find the default SLURM output file.
            
        Returns:
            Job output as a string.
            
        Raises:
            ValueError: If the output file cannot be found or read.
        """
        if output_file:
            # Use provided output file
            if not os.path.exists(output_file):
                raise ValueError(f"Output file not found: {output_file}")
                
            with open(output_file, "r") as f:
                return f.read()
        else:
            # Try to find default SLURM output file
            # Default format: slurm-{job_id}.out
            default_output = f"slurm-{job_id}.out"
            
            if os.path.exists(default_output):
                with open(default_output, "r") as f:
                    return f.read()
            else:
                raise ValueError(f"Default output file not found: {default_output}")
    
    def get_job_error(self, job_id: str, error_file: Optional[str] = None) -> str:
        """
        Get the error output of a SLURM job.
        
        Args:
            job_id: Job ID to get error output for.
            error_file: Optional path to the error file. If not provided,
                       the function will try to find the default SLURM error file.
            
        Returns:
            Job error output as a string.
            
        Raises:
            ValueError: If the error file cannot be found or read.
        """
        if error_file:
            # Use provided error file
            if not os.path.exists(error_file):
                raise ValueError(f"Error file not found: {error_file}")
                
            with open(error_file, "r") as f:
                return f.read()
        else:
            # Try to find default SLURM error file
            # Default format: slurm-{job_id}.err
            default_error = f"slurm-{job_id}.err"
            
            if os.path.exists(default_error):
                with open(default_error, "r") as f:
                    return f.read()
            else:
                # Some SLURM configurations combine stdout and stderr
                # Try the output file as a fallback
                try:
                    return self.get_job_output(job_id)
                except ValueError:
                    raise ValueError(f"Default error file not found: {default_error}")
