"""
RCSB PDB API client for retrieving protein structures.

This module provides functions to search for existing crystal structures.
"""

import logging
from typing import Dict, Any, List, Optional
import httpx
import json

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.rcsb")


class RCSBClient:
    """Client for interacting with RCSB PDB API."""
    
    def __init__(self, base_url: str = None, timeout: int = 30):
        """
        Initialize RCSB client.
        
        Args:
            base_url: Base URL for RCSB API
            timeout: Request timeout in seconds
        """
        self.base_url = base_url or config.RCSB_API_URL
        self.search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.data_url = "https://data.rcsb.org/rest/v1/core"
        self.timeout = timeout
        self.session = httpx.Client(timeout=timeout)
        
    def search_by_uniprot(self, uniprot_id: str, limit: int = 10, with_ligands: bool = False) -> List[Dict[str, Any]]:
        """
        Search for PDB structures by UniProt ID.
        
        Args:
            uniprot_id: UniProt accession ID
            limit: Maximum number of structures to return
            
        Returns:
            List of PDB entries
        """
        try:
            logger.debug(f"Searching RCSB for structures of {uniprot_id}")
            
            # Build search query
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id
                    }
                },
                "request_options": {
                    "return_all_hits": False,
                    "results_content_type": ["experimental"],
                    "sort": [
                        {
                            "sort_by": "score",
                            "direction": "desc"
                        }
                    ],
                    "scoring_strategy": "combined"
                },
                "return_type": "entry"
            }
            
            response = self.session.post(
                self.search_url,
                json=query,
                params={"rows": limit}
            )
            response.raise_for_status()
            
            data = response.json()
            pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]
            
            # Get detailed info for each PDB
            structures = []
            for pdb_id in pdb_ids[:limit]:
                structure_info = self.get_structure_info(pdb_id)
                if structure_info:
                    # Optionally add detailed ligand information
                    if with_ligands:
                        ligand_info = self.get_structure_ligands(pdb_id)
                        structure_info["ligands"] = ligand_info
                        structure_info["has_drug_like_ligand"] = any(
                            self._is_potentially_drug_like(lig) for lig in ligand_info
                        )
                    structures.append(structure_info)
            
            logger.info(f"Found {len(structures)} structures for {uniprot_id}")
            return structures
            
        except Exception as e:
            logger.error(f"Error searching RCSB: {e}")
            return []
    
    def get_structure_info(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a PDB structure using REST API.
        
        Args:
            pdb_id: 4-character PDB ID
            
        Returns:
            Structure information or None
        """
        try:
            logger.debug(f"Fetching PDB structure info: {pdb_id}")
            
            # Use REST API instead of GraphQL
            entry_url = f"{self.data_url}/entry/{pdb_id.upper()}"
            
            response = self.session.get(entry_url)
            
            if response.status_code == 404:
                logger.warning(f"PDB {pdb_id} not found")
                return None
            
            response.raise_for_status()
            entry_data = response.json()
            
            # Get experimental info
            exptl_url = f"{self.data_url}/exptl/{pdb_id.upper()}"
            exptl_response = self.session.get(exptl_url)
            exptl_data = exptl_response.json() if exptl_response.status_code == 200 else []
            
            # Parse the structure information
            structure_info = self._parse_structure_info(pdb_id, entry_data, exptl_data)
            return structure_info
            
        except Exception as e:
            logger.error(f"Error fetching PDB {pdb_id}: {e}")
            return None
    
    def search_by_sequence(self, sequence: str, identity_cutoff: float = 0.9) -> List[Dict[str, Any]]:
        """
        Search for structures by sequence similarity.
        
        Args:
            sequence: Protein sequence
            identity_cutoff: Minimum sequence identity (0-1)
            
        Returns:
            List of similar structures
        """
        try:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff": 1e-5,
                        "identity_cutoff": identity_cutoff,
                        "sequence_type": "protein",
                        "value": sequence
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "scoring_strategy": "sequence"
                }
            }
            
            response = self.session.post(
                self.search_url,
                json=query,
                params={"rows": 10}
            )
            response.raise_for_status()
            
            data = response.json()
            pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]
            
            structures = []
            for pdb_id in pdb_ids:
                structure_info = self.get_structure_info(pdb_id)
                if structure_info:
                    structures.append(structure_info)
            
            return structures
            
        except Exception as e:
            logger.error(f"Error in sequence search: {e}")
            return []
    
    def _parse_structure_info(self, pdb_id: str, entry_data: Dict[str, Any], 
                             exptl_data: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Parse RCSB entry data into simplified format.
        
        Args:
            pdb_id: PDB ID
            entry_data: Entry data from REST API
            exptl_data: Experimental data from REST API
            
        Returns:
            Simplified structure information
        """
        # Get basic info
        struct_info = entry_data.get("struct", {})
        title = struct_info.get("title", "Unknown")
        
        # Get experimental method and resolution
        method = "Unknown"
        resolution = None
        
        if exptl_data and isinstance(exptl_data, list) and len(exptl_data) > 0:
            exptl = exptl_data[0]
            method = exptl.get("method", "Unknown")
            
        # Try to get resolution from different sources
        rcsb_info = entry_data.get("rcsb_entry_info", {})
        resolution = rcsb_info.get("resolution_combined")
        
        # Check for ligands (simplified)
        has_ligand = False
        ligands = []
        
        # Get non-polymer entities (ligands)
        nonpolymer_info = entry_data.get("nonpolymer_entities", [])
        if nonpolymer_info:
            has_ligand = True
            for entity in nonpolymer_info[:3]:  # Limit to first 3 ligands
                comp = entity.get("nonpolymer_comp", {})
                if comp:
                    ligands.append({
                        "name": comp.get("chem_comp", {}).get("name", "Unknown"),
                        "formula": comp.get("chem_comp", {}).get("formula")
                    })
        
        return {
            "pdb_id": pdb_id.upper(),
            "title": title[:100] if title else "Unknown",  # Truncate long titles
            "method": method,
            "resolution": resolution,
            "has_ligand": has_ligand,
            "ligands": ligands,
            "download_url": f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        }
    
    def get_structure_ligands(self, pdb_id: str) -> List[Dict[str, Any]]:
        """
        Get detailed ligand information for a PDB structure.
        
        Args:
            pdb_id: 4-character PDB ID
            
        Returns:
            List of ligand information dictionaries
        """
        try:
            logger.debug(f"Fetching ligand information for {pdb_id}")
            
            # Get chemical component information
            chem_comp_url = f"{self.data_url}/nonpolymer_entity/{pdb_id.upper()}"
            response = self.session.get(chem_comp_url)
            
            if response.status_code == 404:
                logger.debug(f"No ligand data found for {pdb_id}")
                return []
            
            ligands = []
            if response.status_code == 200:
                nonpolymer_data = response.json()
                
                for entity in nonpolymer_data:
                    comp_info = entity.get("rcsb_nonpolymer_entity", {})
                    comp_id = comp_info.get("comp_id", "")
                    
                    if not comp_id:
                        continue
                    
                    # Get detailed chemical component info
                    chem_detail = self._get_chemical_component_details(comp_id)
                    
                    ligand_data = {
                        "id": comp_id,
                        "name": comp_info.get("pdbx_description", "Unknown"),
                        "molecular_weight": self._safe_float(comp_info.get("formula_weight")),
                        "formula": comp_info.get("formula"),
                        "heavy_atom_count": self._count_heavy_atoms(comp_info.get("formula", "")),
                        "type": chem_detail.get("type", "UNKNOWN"),
                        "smiles": chem_detail.get("smiles"),
                        "inchi": chem_detail.get("inchi"),
                    }
                    
                    ligands.append(ligand_data)
                    
            logger.debug(f"Found {len(ligands)} ligands in {pdb_id}")
            return ligands
            
        except Exception as e:
            logger.error(f"Error fetching ligands for {pdb_id}: {e}")
            return []
    
    def _safe_float(self, value) -> float:
        """
        Safely convert a value to float, handling None, lists, and strings.
        
        Args:
            value: Value to convert
            
        Returns:
            Float value or 0 if conversion fails
        """
        if value is None:
            return 0
        if isinstance(value, list):
            return float(value[0]) if value else 0
        if isinstance(value, str):
            try:
                return float(value)
            except ValueError:
                return 0
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0
    
    def _get_chemical_component_details(self, comp_id: str) -> Dict[str, Any]:
        """
        Get detailed chemical component information.
        
        Args:
            comp_id: Chemical component ID (3-letter code)
            
        Returns:
            Dictionary with chemical details
        """
        try:
            # Try to get from RCSB chemical component dictionary
            chem_comp_url = f"https://files.rcsb.org/ligands/view/{comp_id}.json"
            response = self.session.get(chem_comp_url, timeout=5)
            
            if response.status_code == 200:
                data = response.json()
                return {
                    "type": data.get("type", "UNKNOWN"),
                    "smiles": data.get("pdbx_smiles", {}).get("canonical"),
                    "inchi": data.get("pdbx_inchi", {}).get("inchi"),
                }
        except:
            pass
        
        return {"type": "UNKNOWN"}
    
    def _count_heavy_atoms(self, formula: str) -> int:
        """
        Count heavy atoms (non-hydrogen) from chemical formula.
        
        Args:
            formula: Chemical formula string
            
        Returns:
            Number of heavy atoms
        """
        if not formula:
            return 0
        
        # Simple parser for chemical formulas
        import re
        
        # Remove hydrogen atoms
        formula_no_h = re.sub(r'H\d*', '', formula)
        
        # Count remaining atoms
        atom_pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(atom_pattern, formula_no_h)
        
        count = 0
        for element, number in matches:
            if number:
                count += int(number)
            else:
                count += 1
        
        return count
    
    def _is_potentially_drug_like(self, ligand: Dict[str, Any]) -> bool:
        """
        Quick check if a ligand is potentially drug-like.
        
        Args:
            ligand: Ligand information dictionary
            
        Returns:
            True if potentially drug-like
        """
        # Quick checks without importing the full filter module
        mw = ligand.get("molecular_weight", 0)
        
        # Handle case where molecular_weight might be a list or None
        if isinstance(mw, list):
            mw = mw[0] if mw else 0
        elif mw is None:
            mw = 0
        
        ligand_id = ligand.get("id", "").upper()
        
        # Check molecular weight range
        if not (150 < mw < 800):
            return False
        
        # Exclude common artifacts
        artifacts = {'HOH', 'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'CL', 'NA', 'K', 'MG', 'CA', 'ZN'}
        if ligand_id in artifacts:
            return False
        
        # Exclude nucleotides
        nucleotides = {'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP'}
        if ligand_id in nucleotides:
            return False
        
        return True
    
    def search_by_uniprot_with_ligands(self, uniprot_id: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search for PDB structures with drug-like ligands by UniProt ID.
        
        This is a convenience method that combines search with ligand filtering.
        
        Args:
            uniprot_id: UniProt accession ID
            limit: Maximum number of structures to return
            
        Returns:
            List of PDB entries with ligand information
        """
        # Get structures with ligand information
        structures = self.search_by_uniprot(uniprot_id, limit=limit * 2, with_ligands=True)  # Get more initially
        
        # Filter to only those with drug-like ligands
        filtered = [s for s in structures if s.get("has_drug_like_ligand")]
        
        logger.info(f"Found {len(filtered)} structures with drug-like ligands for {uniprot_id}")
        return filtered[:limit]  # Return up to limit
    
    def download_structure(self, pdb_id: str, output_path: str = None) -> Optional[str]:
        """
        Download a PDB structure file.
        
        Args:
            pdb_id: 4-character PDB ID
            output_path: Path to save the file (optional)
            
        Returns:
            Path to downloaded file or PDB content as string
        """
        try:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            response = self.session.get(url)
            response.raise_for_status()
            
            pdb_content = response.text
            
            if output_path:
                with open(output_path, 'w') as f:
                    f.write(pdb_content)
                logger.info(f"Downloaded PDB {pdb_id} to {output_path}")
                return output_path
            else:
                return pdb_content
                
        except Exception as e:
            logger.error(f"Error downloading PDB {pdb_id}: {e}")
            return None
    
    def __del__(self):
        """Clean up session on deletion."""
        if hasattr(self, 'session'):
            self.session.close()


# Singleton instance
_rcsb_client = None

def get_rcsb_client() -> RCSBClient:
    """Get or create RCSB client singleton."""
    global _rcsb_client
    if _rcsb_client is None:
        _rcsb_client = RCSBClient()
    return _rcsb_client
