"""
API interfaces for protein databases and compound libraries.

This module provides unified interfaces to:
- UniProt (protein sequences and annotations)
- RCSB PDB (protein structures)
- ChEMBL (bioactivity data and known compounds)
- AlphaFold DB (predicted structures)
"""

import asyncio
import json
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime
import httpx
import aiohttp
from pathlib import Path

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.protein_apis")


class UniProtAPI:
    """Interface to UniProt REST API for protein information."""
    
    def __init__(self, base_url: str = None):
        self.base_url = base_url or config.UNIPROT_API_URL
        self.session = None
        
    async def __aenter__(self):
        self.session = aiohttp.ClientSession()
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.session:
            await self.session.close()
    
    async def search_protein(self, 
                           query: str, 
                           organism: str = "Homo sapiens",
                           limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search for proteins by gene name, protein name, or UniProt ID.
        
        Args:
            query: Search term (e.g., "KRAS", "P01116")
            organism: Organism name (default: "Homo sapiens")
            limit: Maximum number of results
            
        Returns:
            List of protein entries
        """
        # Build search query
        search_query = f"({query})"
        if organism:
            search_query += f" AND (organism_name:\"{organism}\")"
        
        params = {
            "query": search_query,
            "format": "json",
            "size": limit,
            "fields": "accession,id,gene_names,protein_name,organism_name,sequence,ec,cc_function,cc_disease,ft_variant"
        }
        
        url = f"{self.base_url}/search"
        
        try:
            if not self.session:
                self.session = aiohttp.ClientSession()
                
            async with self.session.get(url, params=params) as response:
                response.raise_for_status()
                data = await response.json()
                
                results = []
                for entry in data.get("results", []):
                    protein_info = self._parse_uniprot_entry(entry)
                    results.append(protein_info)
                
                logger.info(f"Found {len(results)} proteins for query: {query}")
                return results
                
        except Exception as e:
            logger.error(f"UniProt search failed: {e}")
            raise
    
    async def get_protein_by_id(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed protein information by UniProt ID.
        
        Args:
            uniprot_id: UniProt accession (e.g., "P01116")
            
        Returns:
            Protein information dictionary
        """
        url = f"{self.base_url}/{uniprot_id}.json"
        
        try:
            if not self.session:
                self.session = aiohttp.ClientSession()
                
            async with self.session.get(url) as response:
                if response.status == 404:
                    logger.warning(f"Protein {uniprot_id} not found")
                    return None
                    
                response.raise_for_status()
                data = await response.json()
                
                protein_info = self._parse_uniprot_entry(data)
                logger.info(f"Retrieved protein info for {uniprot_id}")
                return protein_info
                
        except Exception as e:
            logger.error(f"Failed to get protein {uniprot_id}: {e}")
            raise
    
    def _parse_uniprot_entry(self, entry: Dict[str, Any]) -> Dict[str, Any]:
        """Parse UniProt entry to extract relevant information."""
        # Extract gene names
        gene_names = entry.get("genes", [])
        primary_gene = gene_names[0].get("geneName", {}).get("value", "") if gene_names else ""
        
        # Extract protein name
        protein_names = entry.get("proteinDescription", {}).get("recommendedName", {})
        protein_name = protein_names.get("fullName", {}).get("value", "")
        
        # Extract function
        functions = entry.get("comments", [])
        function_text = ""
        disease_associations = []
        
        for comment in functions:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")
            elif comment.get("commentType") == "DISEASE":
                disease = comment.get("disease", {})
                disease_associations.append({
                    "name": disease.get("diseaseId", ""),
                    "description": disease.get("description", "")
                })
        
        # Extract sequence
        sequence_info = entry.get("sequence", {})
        sequence = sequence_info.get("value", "")
        
        # Extract variants (mutations)
        features = entry.get("features", [])
        variants = []
        for feature in features:
            if feature.get("type") == "VARIANT":
                location = feature.get("location", {})
                start = location.get("start", {}).get("value", "")
                end = location.get("end", {}).get("value", "")
                original = feature.get("alternativeSequence", {}).get("originalSequence", "")
                variation = feature.get("alternativeSequence", {}).get("alternativeSequences", [""])[0]
                
                if start and original and variation:
                    variants.append(f"{original}{start}{variation}")
        
        return {
            "uniprot_id": entry.get("primaryAccession", ""),
            "entry_name": entry.get("uniProtkbId", ""),
            "protein_name": protein_name,
            "gene_name": primary_gene,
            "organism": entry.get("organism", {}).get("scientificName", ""),
            "sequence": sequence,
            "sequence_length": len(sequence),
            "function": function_text,
            "disease_associations": disease_associations,
            "variants": variants,
            "ec_number": entry.get("proteinDescription", {}).get("ecNumbers", []),
        }


class ChEMBLAPI:
    """Interface to ChEMBL REST API for bioactivity data."""
    
    def __init__(self, base_url: str = None):
        self.base_url = base_url or config.CHEMBL_API_URL
        self.session = None
        
    async def __aenter__(self):
        self.session = aiohttp.ClientSession()
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.session:
            await self.session.close()
    
    async def search_compounds_by_target(self, 
                                        target_name: str,
                                        min_pchembl: float = 5.0,
                                        limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search for compounds active against a protein target.
        
        Args:
            target_name: Protein/gene name (e.g., "KRAS")
            min_pchembl: Minimum pChEMBL value (higher = more potent)
            limit: Maximum number of compounds
            
        Returns:
            List of compound information
        """
        compounds = []
        
        try:
            if not self.session:
                self.session = aiohttp.ClientSession()
            
            # First, find the target in ChEMBL
            target_search_url = f"{self.base_url}/target/search.json"
            params = {"q": target_name, "limit": 5}
            
            async with self.session.get(target_search_url, params=params) as response:
                response.raise_for_status()
                data = await response.json()
                targets = data.get("targets", [])
                
                if not targets:
                    logger.warning(f"No ChEMBL target found for {target_name}")
                    return compounds
                
                # Use the first matching target
                target_chembl_id = targets[0].get("target_chembl_id")
                logger.info(f"Found ChEMBL target: {target_chembl_id}")
            
            # Get bioactivity data for this target
            activity_url = f"{self.base_url}/activity.json"
            params = {
                "target_chembl_id": target_chembl_id,
                "type": "IC50",
                "limit": limit * 2  # Get more to filter later
            }
            
            async with self.session.get(activity_url, params=params) as response:
                response.raise_for_status()
                data = await response.json()
                activities = data.get("activities", [])
            
            # Process activities to extract unique compounds
            seen_compounds = set()
            
            for activity in activities:
                pchembl = activity.get("pchembl_value")
                if pchembl and float(pchembl) >= min_pchembl:
                    compound_chembl_id = activity.get("molecule_chembl_id")
                    
                    if compound_chembl_id not in seen_compounds:
                        seen_compounds.add(compound_chembl_id)
                        
                        compound_info = {
                            "compound_id": compound_chembl_id,
                            "name": activity.get("molecule_pref_name", ""),
                            "smiles": activity.get("canonical_smiles", ""),
                            "pchembl_value": float(pchembl),
                            "standard_value": activity.get("standard_value"),
                            "standard_units": activity.get("standard_units", "nM"),
                            "assay_type": activity.get("assay_type", ""),
                            "target_chembl_id": target_chembl_id
                        }
                        
                        # Convert pChEMBL to approximate IC50 in nM
                        if pchembl:
                            compound_info["ic50_nm"] = 10 ** (9 - float(pchembl))
                        
                        compounds.append(compound_info)
                        
                        if len(compounds) >= limit:
                            break
            
            logger.info(f"Found {len(compounds)} active compounds for {target_name}")
            return compounds
            
        except Exception as e:
            logger.error(f"ChEMBL search failed: {e}")
            return compounds


class RCSBAPI:
    """Interface to RCSB PDB for structure and ligand information."""
    
    def __init__(self):
        self.graphql_url = config.RCSB_API_URL
        self.session = None
        
    async def __aenter__(self):
        self.session = aiohttp.ClientSession()
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.session:
            await self.session.close()
    
    async def search_structures_by_uniprot(self, 
                                          uniprot_id: str,
                                          with_ligand: bool = True) -> List[Dict[str, Any]]:
        """
        Search for PDB structures by UniProt ID.
        
        Args:
            uniprot_id: UniProt accession
            with_ligand: Only return structures with bound ligands
            
        Returns:
            List of PDB structure information
        """
        query = """
        query searchByUniProt($uniprot_id: String!) {
          polymer_entities(
            where: {
              rcsb_polymer_entity_container_identifiers: {
                uniprot_ids: {
                  _has_key: $uniprot_id
                }
              }
            }
          ) {
            rcsb_polymer_entity_container_identifiers {
              entry_id
              entity_id
            }
            entity_poly {
              pdbx_seq_one_letter_code_can
            }
            rcsb_entity_source_organism {
              scientific_name
            }
          }
        }
        """
        
        try:
            if not self.session:
                self.session = aiohttp.ClientSession()
                
            async with self.session.post(
                self.graphql_url,
                json={"query": query, "variables": {"uniprot_id": uniprot_id}}
            ) as response:
                response.raise_for_status()
                data = await response.json()
                
                structures = []
                for entity in data.get("data", {}).get("polymer_entities", []):
                    entry_id = entity["rcsb_polymer_entity_container_identifiers"]["entry_id"]
                    
                    # Get additional structure details
                    structure_info = await self.get_structure_details(entry_id)
                    
                    if with_ligand and not structure_info.get("has_ligand"):
                        continue
                    
                    structures.append(structure_info)
                
                logger.info(f"Found {len(structures)} PDB structures for {uniprot_id}")
                return structures
                
        except Exception as e:
            logger.error(f"RCSB search failed: {e}")
            return []
    
    async def get_structure_details(self, pdb_id: str) -> Dict[str, Any]:
        """Get detailed information about a PDB structure."""
        query = """
        query getStructureDetails($entry_id: String!) {
          entry(entry_id: $entry_id) {
            rcsb_id
            struct {
              title
            }
            rcsb_entry_info {
              resolution_combined
              experimental_method
              structure_determination_methodology
            }
            polymer_entities {
              entity_poly {
                pdbx_seq_one_letter_code_can
              }
            }
            nonpolymer_entities {
              nonpolymer_comp {
                chem_comp {
                  id
                  name
                  formula
                  formula_weight
                }
              }
            }
          }
        }
        """
        
        try:
            async with self.session.post(
                self.graphql_url,
                json={"query": query, "variables": {"entry_id": pdb_id.upper()}}
            ) as response:
                response.raise_for_status()
                data = await response.json()
                
                entry = data.get("data", {}).get("entry", {})
                
                # Extract ligand information
                ligands = []
                for entity in entry.get("nonpolymer_entities", []):
                    comp = entity.get("nonpolymer_comp", {}).get("chem_comp", {})
                    if comp.get("id") and comp["id"] not in ["HOH", "WAT"]:  # Exclude water
                        ligands.append({
                            "id": comp.get("id"),
                            "name": comp.get("name"),
                            "formula": comp.get("formula"),
                            "weight": comp.get("formula_weight")
                        })
                
                return {
                    "pdb_id": pdb_id.upper(),
                    "title": entry.get("struct", {}).get("title", ""),
                    "resolution": entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                    "method": entry.get("rcsb_entry_info", {}).get("experimental_method", ""),
                    "has_ligand": len(ligands) > 0,
                    "ligands": ligands
                }
                
        except Exception as e:
            logger.error(f"Failed to get details for {pdb_id}: {e}")
            return {"pdb_id": pdb_id, "error": str(e)}


class AlphaFoldAPI:
    """Interface to AlphaFold Database."""
    
    def __init__(self):
        self.base_url = config.ALPHAFOLD_API_URL
        self.session = None
        
    async def __aenter__(self):
        self.session = aiohttp.ClientSession()
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.session:
            await self.session.close()
    
    async def check_structure_availability(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Check if AlphaFold structure is available for a UniProt ID.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Structure information if available
        """
        url = f"{self.base_url}/prediction/{uniprot_id}"
        
        try:
            if not self.session:
                self.session = aiohttp.ClientSession()
                
            async with self.session.get(url) as response:
                if response.status == 404:
                    logger.info(f"No AlphaFold structure for {uniprot_id}")
                    return None
                    
                response.raise_for_status()
                data = await response.json()
                
                return {
                    "uniprot_id": uniprot_id,
                    "pdb_url": data[0].get("pdbUrl"),
                    "pae_url": data[0].get("paeUrl"),
                    "confidence_version": data[0].get("confidenceVersion"),
                    "model_version": data[0].get("modelCreatedDate")
                }
                
        except Exception as e:
            logger.error(f"AlphaFold API error: {e}")
            return None


# Convenience functions for synchronous usage
def search_uniprot(query: str, organism: str = "Homo sapiens") -> List[Dict[str, Any]]:
    """Synchronous wrapper for UniProt search."""
    async def _search():
        async with UniProtAPI() as api:
            return await api.search_protein(query, organism)
    return asyncio.run(_search())


def get_chembl_compounds(target: str, min_pchembl: float = 5.0) -> List[Dict[str, Any]]:
    """Synchronous wrapper for ChEMBL search."""
    async def _search():
        async with ChEMBLAPI() as api:
            return await api.search_compounds_by_target(target, min_pchembl)
    return asyncio.run(_search())


def search_pdb_structures(uniprot_id: str) -> List[Dict[str, Any]]:
    """Synchronous wrapper for PDB search."""
    async def _search():
        async with RCSBAPI() as api:
            return await api.search_structures_by_uniprot(uniprot_id)
    return asyncio.run(_search())
