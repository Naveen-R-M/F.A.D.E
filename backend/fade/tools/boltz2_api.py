"""
Boltz2 API client wrapper.

This module provides a wrapper that redirects to the HPC Boltz2 implementation.
It maintains backward compatibility with existing code references.
"""

from fade.tools.hpc_boltz2 import get_hpc_boltz2_client, HPCBoltz2Client
from fade.utils import get_logger

logger = get_logger("tools.boltz2")


class Boltz2Client:
    """
    Wrapper client for Boltz2 that uses HPC implementation.
    
    This maintains compatibility with existing code that calls get_boltz2_client().
    """
    
    def __init__(self):
        """Initialize Boltz2 client by getting HPC client."""
        self.hpc_client = get_hpc_boltz2_client()
    
    def predict_structure(self, sequence: str, job_name: str = None, **kwargs):
        """
        Predict protein structure using Boltz2.
        
        Args:
            sequence: Protein sequence
            job_name: Optional job name
            **kwargs: Additional arguments passed to HPC client
            
        Returns:
            Dictionary with job_id and status
        """
        return self.hpc_client.predict_structure(
            sequence=sequence,
            job_name=job_name,
            **kwargs
        )
    
    def wait_for_structure(self, job_id: str, max_wait: int = 1800, poll_interval: int = 30):
        """
        Wait for structure prediction to complete.
        
        Args:
            job_id: Job ID to wait for
            max_wait: Maximum wait time in seconds
            poll_interval: Polling interval in seconds
            
        Returns:
            PDB content as string or None if failed
        """
        return self.hpc_client.wait_for_structure(
            job_id=job_id,
            max_wait=max_wait,
            poll_interval=poll_interval
        )
    
    def check_status(self, job_id: str):
        """
        Check status of a Boltz2 job.
        
        Args:
            job_id: Job ID to check
            
        Returns:
            Status dictionary
        """
        return self.hpc_client.check_job_status(job_id)


# Singleton instance
_boltz2_client = None

def get_boltz2_client() -> Boltz2Client:
    """
    Get or create Boltz2 client singleton.
    
    This function is referenced in existing code and maintains compatibility.
    """
    global _boltz2_client
    if _boltz2_client is None:
        _boltz2_client = Boltz2Client()
    return _boltz2_client
