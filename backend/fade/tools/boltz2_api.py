"""
Boltz-2 API client for protein structure prediction and screening.

This module provides interface to Boltz-2 for:
1. Structure prediction from sequence
2. Binding affinity prediction (screening)
"""

import logging
from typing import Dict, Any, Optional, List
import httpx
import time
import json
from pathlib import Path

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.boltz2")


class Boltz2Client:
    """Client for interacting with Boltz-2 server."""
    
    def __init__(self, server_url: str = None, api_key: str = None, timeout: int = 300):
        """
        Initialize Boltz-2 client.
        
        Args:
            server_url: Boltz-2 server URL
            api_key: API key if required
            timeout: Request timeout in seconds
        """
        self.server_url = server_url or config.BOLTZ2_SERVER_URL
        self.api_key = api_key or config.BOLTZ2_API_KEY
        self.timeout = timeout
        self.session = httpx.Client(timeout=timeout)
        
    def predict_structure(self, sequence: str, job_name: str = None) -> Dict[str, Any]:
        """
        Predict protein structure from amino acid sequence.
        
        Args:
            sequence: Amino acid sequence (single letter code)
            job_name: Optional name for the job
            
        Returns:
            Dictionary with job_id and status
        """
        # Validate sequence
        if not sequence:
            raise ValueError("Sequence cannot be empty")
            
        if len(sequence) > config.BOLTZ2_MAX_SEQUENCE_LENGTH:
            raise ValueError(f"Sequence too long: {len(sequence)} > {config.BOLTZ2_MAX_SEQUENCE_LENGTH}")
        
        # Prepare request
        headers = self._get_headers()
        payload = {
            "sequence": sequence,
            "job_name": job_name or f"structure_{int(time.time())}",
            "mode": "structure_prediction"
        }
        
        try:
            logger.info(f"Submitting structure prediction job for sequence of length {len(sequence)}")
            response = self.session.post(
                f"{self.server_url}/predict/structure",
                json=payload,
                headers=headers
            )
            response.raise_for_status()
            
            result = response.json()
            logger.info(f"Job submitted: {result.get('job_id')}")
            return result
            
        except Exception as e:
            logger.error(f"Error submitting structure prediction: {e}")
            raise
    
    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check status of a Boltz-2 job.
        
        Args:
            job_id: Job identifier
            
        Returns:
            Job status information
        """
        headers = self._get_headers()
        
        try:
            response = self.session.get(
                f"{self.server_url}/job/{job_id}/status",
                headers=headers
            )
            response.raise_for_status()
            
            return response.json()
            
        except Exception as e:
            logger.error(f"Error checking job status: {e}")
            raise
    
    def get_structure_result(self, job_id: str) -> Optional[str]:
        """
        Retrieve predicted structure in PDB format.
        
        Args:
            job_id: Job identifier
            
        Returns:
            PDB structure as string or None if not ready
        """
        headers = self._get_headers()
        
        try:
            response = self.session.get(
                f"{self.server_url}/job/{job_id}/result",
                headers=headers
            )
            
            if response.status_code == 202:  # Job still running
                logger.info(f"Job {job_id} still running")
                return None
                
            response.raise_for_status()
            
            result = response.json()
            return result.get("structure_pdb")
            
        except Exception as e:
            logger.error(f"Error retrieving structure: {e}")
            raise
    
    def wait_for_structure(self, job_id: str, max_wait: int = 3600, poll_interval: int = 30) -> Optional[str]:
        """
        Wait for structure prediction to complete.
        
        Args:
            job_id: Job identifier
            max_wait: Maximum wait time in seconds
            poll_interval: Polling interval in seconds
            
        Returns:
            PDB structure or None if timeout
        """
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            status = self.get_job_status(job_id)
            
            if status.get("status") == "completed":
                logger.info(f"Job {job_id} completed")
                return self.get_structure_result(job_id)
                
            elif status.get("status") == "failed":
                logger.error(f"Job {job_id} failed: {status.get('error')}")
                raise RuntimeError(f"Structure prediction failed: {status.get('error')}")
            
            logger.debug(f"Job {job_id} status: {status.get('status')} - waiting {poll_interval}s")
            time.sleep(poll_interval)
        
        logger.warning(f"Job {job_id} timed out after {max_wait}s")
        return None
    
    def predict_binding_affinity(self, protein_pdb: str, ligand_smiles: str) -> Dict[str, Any]:
        """
        Predict binding affinity between protein and ligand.
        
        Args:
            protein_pdb: Protein structure in PDB format
            ligand_smiles: Ligand SMILES string
            
        Returns:
            Binding affinity prediction results
        """
        headers = self._get_headers()
        payload = {
            "protein_pdb": protein_pdb,
            "ligand_smiles": ligand_smiles,
            "mode": "binding_prediction"
        }
        
        try:
            logger.info(f"Predicting binding affinity for ligand: {ligand_smiles[:30]}...")
            response = self.session.post(
                f"{self.server_url}/predict/binding",
                json=payload,
                headers=headers
            )
            response.raise_for_status()
            
            result = response.json()
            return result
            
        except Exception as e:
            logger.error(f"Error predicting binding affinity: {e}")
            raise
    
    def batch_predict_binding(self, protein_pdb: str, ligand_smiles_list: List[str], 
                            batch_size: int = 10) -> List[Dict[str, Any]]:
        """
        Predict binding affinities for multiple ligands.
        
        Args:
            protein_pdb: Protein structure in PDB format
            ligand_smiles_list: List of SMILES strings
            batch_size: Number of ligands per batch
            
        Returns:
            List of binding predictions
        """
        results = []
        
        for i in range(0, len(ligand_smiles_list), batch_size):
            batch = ligand_smiles_list[i:i + batch_size]
            
            headers = self._get_headers()
            payload = {
                "protein_pdb": protein_pdb,
                "ligands": batch,
                "mode": "batch_binding_prediction"
            }
            
            try:
                logger.info(f"Batch prediction: {i+1}-{min(i+batch_size, len(ligand_smiles_list))} of {len(ligand_smiles_list)}")
                response = self.session.post(
                    f"{self.server_url}/predict/batch_binding",
                    json=payload,
                    headers=headers
                )
                response.raise_for_status()
                
                batch_results = response.json().get("predictions", [])
                results.extend(batch_results)
                
            except Exception as e:
                logger.error(f"Error in batch prediction: {e}")
                # Continue with next batch
                for smiles in batch:
                    results.append({
                        "smiles": smiles,
                        "affinity": None,
                        "error": str(e)
                    })
        
        return results
    
    def _get_headers(self) -> Dict[str, str]:
        """Get headers for API requests."""
        headers = {
            "Content-Type": "application/json"
        }
        
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        
        return headers
    
    def __del__(self):
        """Clean up session on deletion."""
        if hasattr(self, 'session'):
            self.session.close()


# Singleton instance
_boltz2_client = None

def get_boltz2_client() -> Boltz2Client:
    """Get or create Boltz-2 client singleton."""
    global _boltz2_client
    if _boltz2_client is None:
        _boltz2_client = Boltz2Client()
    return _boltz2_client
