"""
Enhanced RCSB PDB API client with direct search capabilities.
OPTIMIZED FOR SMALL MOLECULE COMPLEX DISCOVERY.
FIXED: Correct API endpoints for ligand information.
NO FALLBACKS - Fails immediately on any error.

This module provides direct search for protein structures with small molecule inhibitors.
"""

import logging
from typing import Dict, Any, List, Optional, Tuple
import httpx
import json
import re

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.rcsb_enhanced")


class RCSBEnhancedClient:
    """Enhanced client for RCSB PDB optimized for small molecule drug discovery."""
    
    def __init__(self, timeout: int = 30):
        """
        Initialize enhanced RCSB client.
        
        Args:
            timeout: Request timeout in seconds
        """
        self.search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.data_url = "https://data.rcsb.org/rest/v1/core"
        self.timeout = timeout
        self.session = httpx.Client(timeout=timeout)
        
    def search_by_query(self, query: str, limit: int = 10) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
        """
        Search RCSB for structures with small molecule complexes.
        OPTIMIZED to find drug-like small molecule inhibitors.
        NO ERROR HANDLING - Fails on any error.
        
        Args:
            query: Search query (should include inhibitor/drug keywords)
            limit: Maximum number of structures to return
            
        Returns:
            Tuple of (list of structures, target info dict)
            
        Raises:
            Any exception from the search or data retrieval
        """
        logger.info(f"Searching RCSB for small molecule complexes: {query}")
        # Build advanced search query to prioritize small molecule complexes
        search_query = self._build_small_molecule_search_query(query)
        
        # NO ERROR HANDLING - Let errors propagate
        response = self.session.post(
            self.search_url,
            json=search_query,
            params={"rows": limit * 3}  # Get more initially to filter for small molecules
        )
        response.raise_for_status()  # Will raise on HTTP errors
        
        data = response.json()
        pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]
        
        if not pdb_ids:
            # Try broader search without strict filtering
            logger.info(f"No results with strict criteria, trying broader search")
            broader_query = self._build_broader_search_query(query)
            response = self.session.post(
                self.search_url,
                json=broader_query,
                params={"rows": limit * 2}
            )
            response.raise_for_status()
            data = response.json()
            pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]
            
            if not pdb_ids:
                raise ValueError(f"No PDB structures found for query: {query}")
        
        # Get detailed info for each PDB and extract target information
        structures = []
        target_info = {}
        
        for pdb_id in pdb_ids[:limit * 2]:  # Check more structures to find ones with small molecules
            try:
                structure_data = self.get_enhanced_structure_info(pdb_id)
                if not structure_data:
                    continue
                
                # Only keep structures with drug-like small molecules
                if structure_data["structure"].get("has_drug_like_ligand"):
                    structures.append(structure_data["structure"])
                    
                    # Extract and merge target information from first valid structure
                    if not target_info and structure_data.get("target_info"):
                        target_info = structure_data["target_info"]
                        logger.info(f"Extracted target info from {pdb_id}: {target_info.get('protein_name')}")
                    
                    # Stop if we have enough structures with small molecules
                    if len(structures) >= limit:
                        break
                        
            except Exception as e:
                logger.debug(f"Skipping {pdb_id}: {e}")
                continue
        
        if not structures:
            raise ValueError(f"No structures with small molecule inhibitors found for: {query}")
        
        if not target_info:
            raise ValueError(f"Could not extract target information from PDB structures for: {query}")
        
        logger.info(f"Found {len(structures)} structures with small molecule inhibitors")
        return structures, target_info
    
    def _build_small_molecule_search_query(self, query: str) -> Dict[str, Any]:
        """
        Build search query optimized for small molecule complexes.
        
        Args:
            query: User query string
            
        Returns:
            RCSB search query dictionary
        """
        # Extract keywords that suggest we want small molecules
        small_molecule_keywords = [
            "inhibitor", "antagonist", "agonist", "ligand", 
            "drug", "compound", "small molecule", "complex"
        ]
        
        # Check if query already has small molecule keywords
        has_sm_keyword = any(kw in query.lower() for kw in small_molecule_keywords)
        
        if not has_sm_keyword:
            # Add "inhibitor" to help find small molecule complexes
            query = f"{query} inhibitor"
        
        # Build query with resolution and method filters
        return {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "full_text",
                        "parameters": {
                            "value": query
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "exptl.method",
                            "operator": "in",
                            "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal",
                            "value": 2.5  # Good resolution for drug discovery
                        }
                    }
                ]
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
    
    def _build_broader_search_query(self, query: str) -> Dict[str, Any]:
        """
        Build a broader search query as fallback.
        
        Args:
            query: User query string
            
        Returns:
            Broader RCSB search query
        """
        return {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": query
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
    
    def get_enhanced_structure_info(self, pdb_id: str) -> Dict[str, Any]:
        """
        Get structure info with protein sequence and metadata.
        FOCUSES on identifying small molecule ligands.
        NO ERROR HANDLING - Fails on any error.
        
        Args:
            pdb_id: 4-character PDB ID
            
        Returns:
            Dictionary with structure info and extracted target data
            
        Raises:
            Any exception from API calls
        """
        logger.debug(f"Fetching enhanced info for PDB {pdb_id}")
        
        # Get basic entry data - NO ERROR HANDLING
        entry_url = f"{self.data_url}/entry/{pdb_id.upper()}"
        response = self.session.get(entry_url)
        
        if response.status_code == 404:
            raise ValueError(f"PDB {pdb_id} not found")
        
        response.raise_for_status()
        entry_data = response.json()
        
        # Get polymer entity data (contains sequence and UniProt info) - NO ERROR HANDLING
        polymer_url = f"{self.data_url}/polymer_entity/{pdb_id.upper()}/1"
        polymer_response = self.session.get(polymer_url)
        polymer_response.raise_for_status()
        polymer_data = polymer_response.json()
        
        # Extract structure information
        structure_info = self._parse_structure_info(pdb_id, entry_data)
        
        # Extract target information from polymer data
        target_info = self._extract_target_info(polymer_data, entry_data)
        
        # Validate that we have required data
        if not target_info.get("sequence"):
            raise ValueError(f"No sequence found for PDB {pdb_id}")
        
        # Get ligand information using CORRECT API endpoint
        ligands = self.get_structure_ligands_fixed(pdb_id, entry_data)

        # Filter for small molecule drug-like compounds
        small_molecule_ligands = [
            lig for lig in ligands 
            if self._is_small_molecule_drug_like(lig)
        ]
        
        structure_info["ligands"] = ligands
        structure_info["small_molecule_ligands"] = small_molecule_ligands
        structure_info["has_drug_like_ligand"] = len(small_molecule_ligands) > 0
        
        if small_molecule_ligands:
            logger.debug(f"{pdb_id} has {len(small_molecule_ligands)} small molecule ligands")
        
        return {
            "structure": structure_info,
            "target_info": target_info
        }
    
    def get_structure_ligands_fixed(self, pdb_id: str, entry_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Get detailed ligand information for a PDB structure using CORRECT API endpoints.
        
        Args:
            pdb_id: 4-character PDB ID
            entry_data: Entry data already fetched
            
        Returns:
            List of ligand information dictionaries
        """
        try:
            ligands = []
            
            # CORRECT METHOD: Get non-polymer entity IDs from entry data
            container_ids = entry_data.get("rcsb_entry_container_identifiers", {})
            non_polymer_ids = container_ids.get("non_polymer_entity_ids", [])
            
            if not non_polymer_ids:
                logger.debug(f"No non-polymer entities found for {pdb_id}")
                return []
            
            logger.debug(f"Found non-polymer entity IDs for {pdb_id}: {non_polymer_ids}")
            
            # Fetch data for each non-polymer entity
            for entity_id in non_polymer_ids:
                try:
                    # CORRECT URL with entity ID
                    nonpolymer_url = f"{self.data_url}/nonpolymer_entity/{pdb_id.upper()}/{entity_id}"
                    response = self.session.get(nonpolymer_url)
                    
                    if response.status_code == 200:
                        nonpolymer_data = response.json()
                        
                        # Extract ligand information
                        comp_info = nonpolymer_data.get("rcsb_nonpolymer_entity", {})
                        comp_id = nonpolymer_data.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("nonpolymer_comp_id")
                        
                        if not comp_id:
                            continue
                        
                        # Extract molecular weight
                        mw = self._safe_float(comp_info.get("formula_weight")) * 1000.0 # Convert kDa to g/mol
                        
                        # Get chemical name and formula
                        name = comp_info.get("pdbx_description", "Unknown")
                        formula = comp_info.get("formula", "")
                        
                        # Count heavy atoms from formula if possible
                        heavy_atoms = self._get_heavy_atoms_count(comp_id)
                        
                        ligand_data = {
                            "id": comp_id,
                            "name": name,
                            "molecular_weight": mw,
                            "formula": formula,
                            "heavy_atom_count": heavy_atoms,
                            "entity_id": entity_id
                        }
                        
                        ligands.append(ligand_data)
                        logger.debug(f"  Found ligand: {comp_id} (MW: {mw}, {name})")
                        
                except Exception as e:
                    logger.debug(f"Error fetching non-polymer entity {entity_id}: {e}")
                    continue
            
            return ligands
            
        except Exception as e:
            logger.error(f"Error fetching ligands for {pdb_id}: {e}")
            return []
    
    def get_structure_ligands(self, pdb_id: str) -> List[Dict[str, Any]]:
        """
        DEPRECATED - Use get_structure_ligands_fixed instead.
        This method uses incorrect API endpoint.
        """
        logger.warning(f"Using deprecated get_structure_ligands for {pdb_id}")
        # Try to get entry data and use fixed method
        try:
            entry_url = f"{self.data_url}/entry/{pdb_id.upper()}"
            response = self.session.get(entry_url)
            if response.status_code == 200:
                entry_data = response.json()
                return self.get_structure_ligands_fixed(pdb_id, entry_data)
        except:
            pass
        return []
    
    def _count_heavy_atoms_from_formula(self, formula: str) -> int:
        """
        Count non-hydrogen atoms from chemical formula.
        
        Args:
            formula: Chemical formula string
            
        Returns:
            Number of heavy atoms
        """
        if not formula:
            return 0
        
        import re
        
        # Remove spaces
        formula = formula.replace(" ", "")
        
        # Pattern to match element symbols with optional counts
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        count = 0
        for element, number in matches:
            if element != 'H':  # Exclude hydrogen
                if number:
                    count += int(number)
                else:
                    count += 1
        
        return count
    
    def _get_heavy_atoms_count(self, comp_id: str) -> int:
        """
        Get heavy atoms from chemical component api.
        
        Args:
            comp_id: Component Identifier string
            
        Returns:
            Number of heavy atoms
        """
        chemcomp_url = f"{self.data_url}/chemcomp/{comp_id.upper()}"
        response = self.session.get(chemcomp_url)
        
        if response.status_code == 404:
            raise ValueError(f"Chemical Component for {comp_id} not found")
        
        response.raise_for_status()
        chemcomp_data = response.json()
        heavy_atoms_count = chemcomp_data.get("rcsb_chem_comp_info", {}).get("atom_count_heavy")
        return heavy_atoms_count
    
    def _extract_target_info(self, polymer_data: Dict[str, Any], entry_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract target protein information from PDB data.
        NO FALLBACKS - Must extract required information.
        
        Args:
            polymer_data: Polymer entity data from RCSB
            entry_data: Entry data from RCSB
            
        Returns:
            Target information dictionary
            
        Raises:
            ValueError: If required data cannot be extracted
        """
        target_info = {}
        
        # Extract from polymer entity data - REQUIRED
        if not polymer_data:
            raise ValueError("No polymer data available")
        
        entity_info = polymer_data.get("rcsb_polymer_entity", {})
        
        # Get protein name - REQUIRED
        protein_name = entity_info.get("pdbx_description", "")
        if not protein_name:
            raise ValueError("No protein name found in PDB data")
        target_info["protein_name"] = protein_name
        
        # Get sequence - REQUIRED
        container = polymer_data.get("entity_poly", {})
        sequence = container.get("pdbx_seq_one_letter_code_can", "")
        if not sequence:
            raise ValueError("No sequence found in PDB data")
        
        # Clean sequence (remove newlines, spaces)
        sequence = re.sub(r'[^A-Z]', '', sequence.upper())
        if not sequence:
            raise ValueError("Empty sequence after cleaning")
        
        target_info["sequence"] = sequence
        target_info["sequence_length"] = len(sequence)
        
        # Get UniProt ID if available (optional)
        ref_seq = polymer_data.get("rcsb_polymer_entity_container_identifiers", {})
        ref_ids = ref_seq.get("reference_sequence_identifiers", [])
        for ref in ref_ids:
            if ref.get("database_name") == "UniProt":
                target_info["uniprot_id"] = ref.get("database_accession")
                break
        
        # UniProt ID is optional - don't log if missing
        
        # Get gene name from entity source (optional)
        entity_src = polymer_data.get("entity_src_gen", [])
        if entity_src and isinstance(entity_src, list):
            gene_name = entity_src[0].get("gene_src_common_name") or \
                       entity_src[0].get("pdbx_gene_src_gene")
            if gene_name:
                target_info["gene_name"] = gene_name
        
        # Extract organism if available (optional)
        struct_info = entry_data.get("struct", {})
        target_info["organism"] = struct_info.get("pdbx_descriptor", "")
        
        # Extract any mutations from title (optional)
        title = struct_info.get("title", "")
        mutations = self._extract_mutations_from_title(title)
        if mutations:
            target_info["mutations"] = mutations
        
        return target_info
    
    def _extract_mutations_from_title(self, title: str) -> List[str]:
        """
        Extract mutation information from PDB title.
        
        Args:
            title: PDB structure title
            
        Returns:
            List of mutations found
        """
        mutations = []
        
        # Common mutation patterns (e.g., G12C, T790M)
        mutation_pattern = r'\b([A-Z]\d+[A-Z])\b'
        found_mutations = re.findall(mutation_pattern, title)
        mutations.extend(found_mutations)
        
        # Also check for written out mutations (e.g., "mutant", "variant")
        if re.search(r'\bmutant\b|\bvariant\b|\bmutation\b', title, re.IGNORECASE):
            # Add a flag that this is a mutant structure
            if not mutations:
                mutations.append("MUTANT")
        
        return mutations
    
    def _parse_structure_info(self, pdb_id: str, entry_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse RCSB entry data into simplified structure format.
        
        Args:
            pdb_id: PDB ID
            entry_data: Entry data from REST API
            
        Returns:
            Simplified structure information
        """
        struct_info = entry_data.get("struct", {})
        rcsb_info = entry_data.get("rcsb_entry_info", {})
        
        return {
            "pdb_id": pdb_id.upper(),
            "title": struct_info.get("title", "Unknown")[:100],
            "method": rcsb_info.get("experimental_method", "Unknown"),
            "resolution": rcsb_info.get("resolution_combined"),
            "release_date": rcsb_info.get("initial_release_date"),
            "download_url": f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        }
    
    def _safe_float(self, value) -> float:
        """Safely convert a value to float."""
        if value is None:
            return 0
        if isinstance(value, list):
            return float(value[0]) if value else 0
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0
    
    def _is_small_molecule_drug_like(self, ligand: Dict[str, Any]) -> bool:
        """
        Check if a ligand is a small molecule drug-like compound.
        Stricter criteria for drug discovery.
        
        Args:
            ligand: Ligand information
            
        Returns:
            True if ligand is drug-like small molecule
        """
        mw = ligand.get("molecular_weight", 0)
        ligand_id = ligand.get("id", "").upper()
        
        # Small molecule drug criteria
        if not (150 < mw < 800):  # Typical drug MW range
            return False
        
        # Exclude common artifacts
        artifacts = config.PDB_CRYSTALLIZATION_ARTIFACTS
        if ligand_id in artifacts:
            return False
        
        # Exclude nucleotides
        if ligand_id in config.PDB_NUCLEOTIDES:
            return False
        
        # Exclude common cofactors that aren't drug targets
        cofactors = {"NAD", "NADH", "NADP", "FAD", "FMN", "COA", "SAM", "SAH"}
        if ligand_id in cofactors:
            return False
        
        return True
    
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
_rcsb_enhanced_client = None

def get_rcsb_enhanced_client() -> RCSBEnhancedClient:
    """Get or create enhanced RCSB client singleton."""
    global _rcsb_enhanced_client
    if _rcsb_enhanced_client is None:
        _rcsb_enhanced_client = RCSBEnhancedClient()
    return _rcsb_enhanced_client
