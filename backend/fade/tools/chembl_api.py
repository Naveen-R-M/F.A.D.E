"""
ChEMBL API client for retrieving known compounds and drug information.

This module provides functions to search for known inhibitors and drugs.
"""

import logging
from typing import Dict, Any, List, Optional
import httpx
import time

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.chembl")


class ChEMBLClient:
    """Client for interacting with ChEMBL REST API."""
    
    def __init__(self, base_url: str = None, timeout: int = 30):
        """
        Initialize ChEMBL client.
        
        Args:
            base_url: Base URL for ChEMBL API
            timeout: Request timeout in seconds
        """
        self.base_url = base_url or config.CHEMBL_API_URL
        self.timeout = timeout
        self.session = httpx.Client(timeout=timeout)
        
    def search_by_target(self, uniprot_id: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search for compounds that target a specific protein.
        
        Args:
            uniprot_id: UniProt accession ID of target protein
            limit: Maximum number of compounds to retrieve
            
        Returns:
            List of compound information
        """
        try:
            # First, find the ChEMBL target ID for this UniProt ID
            logger.debug(f"Searching ChEMBL target for UniProt ID: {uniprot_id}")
            
            # Add delay to avoid rate limiting
            time.sleep(1)
            
            headers = {
                "Accept": "application/json",
                "User-Agent": "F.A.D.E/2.0"
            }
            
            # Try the target search endpoint with format parameter
            target_url = f"{self.base_url}/target/search.json"
            target_response = self.session.get(
                target_url,
                params={"q": uniprot_id, "limit": 1},
                headers=headers
            )
            
            # Check if we got JSON response
            content_type = target_response.headers.get("content-type", "")
            if "application/json" not in content_type:
                logger.warning(f"ChEMBL returned non-JSON response: {content_type}")
                return []
            
            target_response.raise_for_status()
            target_data = target_response.json()
            
            if not target_data.get("targets"):
                logger.warning(f"No ChEMBL target found for {uniprot_id}")
                return []
            
            target_chembl_id = target_data["targets"][0]["target_chembl_id"]
            logger.info(f"Found ChEMBL target: {target_chembl_id}")
            
            # Add delay before next request
            time.sleep(1)
            
            # Now get activities for this target
            activities_url = f"{self.base_url}/activity.json"
            activities_response = self.session.get(
                activities_url,
                params={
                    "target_chembl_id": target_chembl_id,
                    "limit": limit,
                    "assay_type": "B",  # Binding assays
                    "pchembl_value__gte": 5  # pChEMBL >= 5 (10 µM or better)
                },
                headers=headers
            )
            
            # Check content type again
            if "application/json" not in activities_response.headers.get("content-type", ""):
                logger.warning("ChEMBL activities endpoint returned non-JSON")
                return []
            
            activities_response.raise_for_status()
            activities_data = activities_response.json()
            
            compounds = self._process_activities(activities_data.get("activities", []))
            logger.info(f"Found {len(compounds)} active compounds for {uniprot_id}")
            
            return compounds
            
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error from ChEMBL: {e.response.status_code}")
            return []
        except Exception as e:
            logger.error(f"Error searching ChEMBL: {e}")
            return []
    
    def get_compound_by_id(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed compound information by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL compound ID (e.g., CHEMBL1234)
            
        Returns:
            Compound information or None
        """
        try:
            logger.debug(f"Fetching compound: {chembl_id}")
            
            headers = {
                "Accept": "application/json",
                "User-Agent": "F.A.D.E/2.0"
            }
            
            response = self.session.get(
                f"{self.base_url}/molecule/{chembl_id}.json",
                headers=headers
            )
            
            if response.status_code == 404:
                logger.warning(f"Compound {chembl_id} not found")
                return None
                
            response.raise_for_status()
            data = response.json()
            
            return self._parse_compound(data)
            
        except Exception as e:
            logger.error(f"Error fetching compound {chembl_id}: {e}")
            return None
    
    def search_approved_drugs(self, target_uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Search for approved drugs targeting a protein.
        
        Args:
            target_uniprot_id: UniProt ID of target
            
        Returns:
            List of approved drugs
        """
        compounds = self.search_by_target(target_uniprot_id)
        
        # Filter for approved drugs
        approved_drugs = []
        for compound in compounds:
            if compound.get("max_phase", 0) == 4:  # Phase 4 = Approved
                approved_drugs.append(compound)
        
        logger.info(f"Found {len(approved_drugs)} approved drugs for {target_uniprot_id}")
        return approved_drugs
    
    def _process_activities(self, activities: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Process activity data into compound information.
        
        Args:
            activities: Raw activity data from ChEMBL
            
        Returns:
            Processed compound information
        """
        compounds_dict = {}
        
        for activity in activities:
            chembl_id = activity.get("molecule_chembl_id")
            if not chembl_id:
                continue
                
            # Get the best (lowest) activity value for each compound
            if chembl_id not in compounds_dict:
                compounds_dict[chembl_id] = {
                    "compound_id": chembl_id,
                    "name": activity.get("molecule_pref_name"),
                    "smiles": activity.get("canonical_smiles"),
                    "binding_affinity": None,
                    "affinity_unit": None,
                    "pchembl_value": activity.get("pchembl_value"),
                    "assay_type": activity.get("assay_type"),
                    "max_phase": activity.get("molecule_max_phase", 0),
                    "clinical_phase": self._phase_to_string(activity.get("molecule_max_phase", 0))
                }
            
            # Update with best activity value
            if activity.get("standard_value") and activity.get("standard_type"):
                current_value = compounds_dict[chembl_id].get("binding_affinity")
                try:
                    new_value = float(activity["standard_value"])
                    
                    if current_value is None or new_value < current_value:
                        compounds_dict[chembl_id]["binding_affinity"] = new_value
                        compounds_dict[chembl_id]["affinity_unit"] = activity.get("standard_units", "nM")
                        compounds_dict[chembl_id]["activity_type"] = activity.get("standard_type")
                except (ValueError, TypeError):
                    pass  # Skip invalid values
        
        return list(compounds_dict.values())
    
    def _parse_compound(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse ChEMBL compound data into simplified format.
        
        Args:
            data: Raw compound data from ChEMBL
            
        Returns:
            Simplified compound information
        """
        molecule_properties = data.get("molecule_properties", {})
        
        return {
            "compound_id": data.get("molecule_chembl_id"),
            "name": data.get("pref_name"),
            "smiles": data.get("molecule_structures", {}).get("canonical_smiles"),
            "molecular_weight": molecule_properties.get("mw_freebase"),
            "logp": molecule_properties.get("alogp"),
            "num_ro5_violations": molecule_properties.get("num_ro5_violations"),
            "max_phase": data.get("max_phase"),
            "clinical_phase": self._phase_to_string(data.get("max_phase", 0)),
            "first_approval": data.get("first_approval"),
            "therapeutic_flags": data.get("therapeutic_flag"),
            "molecule_type": data.get("molecule_type"),
            "natural_product": data.get("natural_product"),
            "oral_bioavailability": molecule_properties.get("oral"),
            "molecule_synonyms": [s["molecule_synonym"] for s in data.get("molecule_synonyms", [])]
        }
    
    def _phase_to_string(self, phase: int) -> str:
        """Convert clinical phase number to string."""
        phase_map = {
            0: "Preclinical",
            1: "Phase I",
            2: "Phase II",
            3: "Phase III",
            4: "Approved",
            -1: "Unknown"
        }
        return phase_map.get(phase, "Unknown")
    
    def __del__(self):
        """Clean up session on deletion."""
        if hasattr(self, 'session'):
            self.session.close()


# Singleton instance
_chembl_client = None

def get_chembl_client() -> ChEMBLClient:
    """Get or create ChEMBL client singleton."""
    global _chembl_client
    if _chembl_client is None:
        _chembl_client = ChEMBLClient()
    return _chembl_client
