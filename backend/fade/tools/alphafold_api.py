"""
AlphaFold Database API client for retrieving pre-computed structures.

This provides fast access to pre-computed AlphaFold structures
before falling back to Boltz-2 prediction.
"""

import logging
from typing import Dict, Any, Optional
import httpx

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.alphafold")


class AlphaFoldDBClient:
    """Client for AlphaFold Database API."""
    
    def __init__(self, base_url: str = None, timeout: int = 30):
        """
        Initialize AlphaFold DB client.
        
        Args:
            base_url: Base URL for AlphaFold API
            timeout: Request timeout in seconds
        """
        self.base_url = base_url or config.ALPHAFOLD_API_URL
        self.timeout = timeout
        self.session = httpx.Client(timeout=timeout)
    
    def get_structure_by_uniprot(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get AlphaFold structure by UniProt ID.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            Structure information or None if not found
        """
        try:
            logger.debug(f"Checking AlphaFold DB for {uniprot_id}")
            
            # Check if structure exists
            response = self.session.get(
                f"{self.base_url}/prediction/{uniprot_id}"
            )
            
            if response.status_code == 404:
                logger.info(f"No AlphaFold structure for {uniprot_id}")
                return None
            
            response.raise_for_status()
            data = response.json()
            
            # Get the latest version
            latest = data[0] if isinstance(data, list) else data
            
            logger.info(f"Found AlphaFold structure for {uniprot_id}, version {latest.get('latestVersion')}")
            
            return {
                "uniprot_id": uniprot_id,
                "pdb_url": latest.get("pdbUrl"),
                "pae_url": latest.get("paeImageUrl"),  # Predicted Aligned Error
                "confidence_url": latest.get("confidenceImageUrl"),
                "model_date": latest.get("modelCreatedDate"),
                "mean_plddt": latest.get("globalMetricValue"),  # Mean confidence score
                "version": latest.get("latestVersion")
            }
            
        except Exception as e:
            logger.error(f"Error fetching AlphaFold structure: {e}")
            return None
    
    def download_pdb(self, uniprot_id: str) -> Optional[str]:
        """
        Download AlphaFold structure in PDB format.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            PDB content as string or None
        """
        try:
            # Get structure info first
            info = self.get_structure_by_uniprot(uniprot_id)
            if not info or not info.get("pdb_url"):
                return None
            
            # Download PDB file
            logger.info(f"Downloading AlphaFold PDB for {uniprot_id}")
            response = self.session.get(info["pdb_url"])
            response.raise_for_status()
            
            pdb_content = response.text
            
            # Add confidence score as a remark
            if info.get("mean_plddt"):
                pdb_lines = pdb_content.split('\n')
                pdb_lines.insert(1, f"REMARK   1 AlphaFold Mean pLDDT: {info['mean_plddt']:.1f}")
                pdb_content = '\n'.join(pdb_lines)
            
            return pdb_content
            
        except Exception as e:
            logger.error(f"Error downloading AlphaFold PDB: {e}")
            return None
    
    def check_confidence(self, uniprot_id: str, min_plddt: float = 70.0) -> bool:
        """
        Check if AlphaFold structure meets confidence threshold.
        
        Args:
            uniprot_id: UniProt accession ID
            min_plddt: Minimum mean pLDDT score (0-100)
            
        Returns:
            True if structure meets threshold
        """
        info = self.get_structure_by_uniprot(uniprot_id)
        
        if not info:
            return False
        
        mean_plddt = info.get("mean_plddt", 0)
        meets_threshold = mean_plddt >= min_plddt
        
        if meets_threshold:
            logger.info(f"AlphaFold structure confidence: {mean_plddt:.1f} (threshold: {min_plddt})")
        else:
            logger.warning(f"AlphaFold structure confidence too low: {mean_plddt:.1f} < {min_plddt}")
        
        return meets_threshold
    
    def __del__(self):
        """Clean up session on deletion."""
        if hasattr(self, 'session'):
            self.session.close()


# Singleton instance
_alphafold_client = None

def get_alphafold_client() -> AlphaFoldDBClient:
    """Get or create AlphaFold DB client singleton."""
    global _alphafold_client
    if _alphafold_client is None:
        _alphafold_client = AlphaFoldDBClient()
    return _alphafold_client
