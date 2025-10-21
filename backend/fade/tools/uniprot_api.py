"""
UniProt API client for protein information retrieval.

This module provides functions to search and retrieve protein data from UniProt.
"""

import logging
from typing import Dict, Any, List, Optional
import httpx
from urllib.parse import quote

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.uniprot")


class UniProtClient:
    """Client for interacting with UniProt REST API."""
    
    def __init__(self, base_url: str = None, timeout: int = 30):
        """
        Initialize UniProt client.
        
        Args:
            base_url: Base URL for UniProt API
            timeout: Request timeout in seconds
        """
        self.base_url = base_url or config.UNIPROT_API_URL
        self.timeout = timeout
        self.session = httpx.Client(timeout=timeout)
        
    def search_protein(self, query: str, organism: str = "human", 
                      limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search for proteins in UniProt.
        
        Args:
            query: Search query (gene name, protein name, or UniProt ID)
            organism: Organism to filter by
            limit: Maximum number of results
            
        Returns:
            List of protein entries
        """
        # Build search query
        search_query = query
        if organism:
            organism_id = "9606" if organism.lower() == "human" else organism
            search_query = f"({query}) AND (organism_id:{organism_id})"
            
        params = {
            "query": search_query,
            "format": "json",
            "size": limit,
            "fields": "accession,id,gene_names,protein_name,organism_name,sequence,length,ft_domain,cc_function,cc_disease"
        }
        
        try:
            logger.debug(f"Searching UniProt for: {search_query}")
            response = self.session.get(
                f"{self.base_url}/search",
                params=params
            )
            response.raise_for_status()
            
            data = response.json()
            results = data.get("results", [])
            
            logger.info(f"Found {len(results)} proteins for query: {query}")
            return results
            
        except Exception as e:
            logger.error(f"Error searching UniProt: {e}")
            return []
    
    def get_protein_by_id(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed protein information by UniProt ID.
        
        Args:
            uniprot_id: UniProt accession ID (e.g., P01112)
            
        Returns:
            Protein entry data or None if not found
        """
        try:
            logger.debug(f"Fetching protein: {uniprot_id}")
            response = self.session.get(
                f"{self.base_url}/{uniprot_id}.json"
            )
            
            if response.status_code == 404:
                logger.warning(f"Protein {uniprot_id} not found")
                return None
                
            response.raise_for_status()
            data = response.json()
            
            logger.info(f"Retrieved protein {uniprot_id}: {data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')}")
            return data
            
        except Exception as e:
            logger.error(f"Error fetching protein {uniprot_id}: {e}")
            return None
    
    def get_protein_sequence(self, uniprot_id: str) -> Optional[str]:
        """
        Get protein sequence in FASTA format.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            Protein sequence string or None
        """
        try:
            response = self.session.get(
                f"{self.base_url}/{uniprot_id}.fasta"
            )
            
            if response.status_code == 404:
                return None
                
            response.raise_for_status()
            fasta_text = response.text
            
            # Extract sequence from FASTA (skip header)
            lines = fasta_text.strip().split('\n')
            sequence = ''.join(lines[1:])
            
            return sequence
            
        except Exception as e:
            logger.error(f"Error fetching sequence for {uniprot_id}: {e}")
            return None
    
    def search_by_gene_name(self, gene_name: str, organism: str = "human") -> List[Dict[str, Any]]:
        """
        Search for proteins by gene name.
        
        Args:
            gene_name: Gene symbol (e.g., KRAS, EGFR)
            organism: Organism name
            
        Returns:
            List of matching proteins
        """
        query = f"gene_exact:{gene_name}"
        return self.search_protein(query, organism)
    
    def parse_protein_info(self, entry: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse UniProt entry into simplified format.
        
        Args:
            entry: Raw UniProt entry
            
        Returns:
            Simplified protein information
        """
        # Extract gene names
        gene_names = entry.get("genes", [])
        primary_gene = gene_names[0]["geneName"]["value"] if gene_names else None
        
        # Extract protein name
        protein_desc = entry.get("proteinDescription", {})
        recommended_name = protein_desc.get("recommendedName", {})
        protein_name = recommended_name.get("fullName", {}).get("value", "Unknown")
        
        # Extract function
        comments = entry.get("comments", [])
        function_comments = [c for c in comments if c.get("commentType") == "FUNCTION"]
        function = function_comments[0]["texts"][0]["value"] if function_comments else None
        
        # Extract disease associations
        disease_comments = [c for c in comments if c.get("commentType") == "DISEASE"]
        diseases = []
        for comment in disease_comments:
            if "disease" in comment:
                diseases.append(comment["disease"]["diseaseId"])
        
        # Get sequence
        sequence_info = entry.get("sequence", {})
        
        return {
            "uniprot_id": entry.get("primaryAccession"),
            "entry_name": entry.get("uniProtkbId"),
            "gene_name": primary_gene,
            "protein_name": protein_name,
            "organism": entry.get("organism", {}).get("scientificName"),
            "sequence": sequence_info.get("value"),
            "sequence_length": sequence_info.get("length"),
            "function": function,
            "disease_associations": diseases,
            "molecular_weight": sequence_info.get("molWeight"),
            "database_cross_refs": entry.get("uniProtKBCrossReferences", [])
        }
    
    def __del__(self):
        """Clean up session on deletion."""
        if hasattr(self, 'session'):
            self.session.close()


# Singleton instance
_uniprot_client = None

def get_uniprot_client() -> UniProtClient:
    """Get or create UniProt client singleton."""
    global _uniprot_client
    if _uniprot_client is None:
        _uniprot_client = UniProtClient()
    return _uniprot_client
