"""
UniProt API client for F.A.D.E framework.
"""

import json
import logging
from typing import Any, Dict, List, Optional, Union

import requests


class UniProtClient:
    """
    Client for interacting with the UniProt API to retrieve protein information.
    """
    
    def __init__(self) -> None:
        """Initialize the UniProt client."""
        self.base_url = "https://rest.uniprot.org/uniprotkb/"
        self.logger = logging.getLogger("fade.uniprot")
    
    def search_protein(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search for proteins in UniProt.
        
        Args:
            query: Search query string.
            limit: Maximum number of results to return.
            
        Returns:
            List of protein entries matching the query.
        """
        search_url = f"{self.base_url}search"
        params = {
            "query": query,
            "format": "json",
            "size": limit
        }
        
        response = requests.get(search_url, params=params)
        
        if response.status_code != 200:
            self.logger.error(f"Error searching UniProt: {response.status_code} {response.text}")
            return []
        
        results = response.json()
        return results.get("results", [])
    
    def get_protein_by_accession(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get protein information by UniProt accession number.
        
        Args:
            accession: UniProt accession number.
            
        Returns:
            Protein information or None if not found.
        """
        entry_url = f"{self.base_url}{accession}"
        params = {"format": "json"}
        
        response = requests.get(entry_url, params=params)
        
        if response.status_code != 200:
            self.logger.error(f"Error getting UniProt entry: {response.status_code} {response.text}")
            return None
        
        return response.json()
    
    def get_sequence(self, accession: str) -> Optional[str]:
        """
        Get protein sequence by UniProt accession number.
        
        Args:
            accession: UniProt accession number.
            
        Returns:
            Protein sequence or None if not found.
        """
        protein = self.get_protein_by_accession(accession)
        
        if not protein:
            return None
        
        return protein.get("sequence", {}).get("value")
    
    def get_protein_by_gene_name(self, gene_name: str, organism: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        Get protein information by gene name and optional organism.
        
        Args:
            gene_name: Gene name.
            organism: Optional organism name or taxonomy ID.
            
        Returns:
            Protein information or None if not found.
        """
        query = f"gene:{gene_name}"
        
        if organism:
            query += f" AND organism:{organism}"
        
        results = self.search_protein(query, limit=1)
        
        if not results:
            return None
        
        return results[0]
    
    def apply_mutation(self, sequence: str, mutation_desc: str) -> str:
        """
        Apply a mutation to a protein sequence.
        
        Args:
            sequence: Original protein sequence.
            mutation_desc: Mutation description in standard format (e.g., "G12D").
            
        Returns:
            Mutated protein sequence.
        """
        # Parse mutation description
        if len(mutation_desc) < 3:
            raise ValueError(f"Invalid mutation format: {mutation_desc}")
        
        original_residue = mutation_desc[0]
        mutated_residue = mutation_desc[-1]
        
        # Extract position (handling cases where position might have multiple digits)
        position_str = mutation_desc[1:-1]
        try:
            position = int(position_str)
        except ValueError:
            raise ValueError(f"Invalid position in mutation description: {mutation_desc}")
        
        # Convert to 0-based index
        index = position - 1
        
        # Validate sequence and mutation
        if index < 0 or index >= len(sequence):
            raise ValueError(f"Mutation position {position} out of range for sequence of length {len(sequence)}")
        
        if sequence[index] != original_residue:
            self.logger.warning(
                f"Original residue mismatch at position {position}: "
                f"expected {original_residue}, found {sequence[index]}"
            )
        
        # Apply mutation
        mutated_sequence = sequence[:index] + mutated_residue + sequence[index+1:]
        
        return mutated_sequence
    
    def generate_fasta(self, accession: str, sequence: str, description: Optional[str] = None) -> str:
        """
        Generate a FASTA format string for a protein sequence.
        
        Args:
            accession: Protein accession number or identifier.
            sequence: Protein sequence.
            description: Optional description for the FASTA header.
            
        Returns:
            FASTA format string.
        """
        header = f">{accession}"
        if description:
            header += f" {description}"
        
        # Format sequence with line breaks
        formatted_sequence = ""
        for i in range(0, len(sequence), 60):
            formatted_sequence += sequence[i:i+60] + "\n"
        
        return f"{header}\n{formatted_sequence}"
    
    def save_fasta(self, file_path: str, accession: str, sequence: str, description: Optional[str] = None) -> None:
        """
        Save a protein sequence to a FASTA file.
        
        Args:
            file_path: Path to the output FASTA file.
            accession: Protein accession number or identifier.
            sequence: Protein sequence.
            description: Optional description for the FASTA header.
        """
        fasta_content = self.generate_fasta(accession, sequence, description)
        
        with open(file_path, "w") as f:
            f.write(fasta_content)
        
        self.logger.info(f"Saved FASTA file to {file_path}")
