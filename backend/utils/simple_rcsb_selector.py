"""
Simple RCSB Target Selector - No LLM Dependencies
Quick implementation to test RCSB functionality without API calls
"""

import os
import json
import requests
from typing import Any, Dict, List, Optional


class SimpleRCSBTargetSelector:
    """Simple RCSB target selector without LLM dependencies."""
    
    def __init__(self, name: str = "simple_rcsb_selector"):
        self.name = name
        
    def process(self, query: str) -> Dict[str, Any]:
        """Process query and return RCSB structure."""
        print(f"Processing query: {query}")
        
        # Simple keyword extraction (no LLM)
        protein_name = self._extract_protein_name(query)
        mutations = self._extract_mutations(query)
        
        # Search RCSB
        structure_info = self._search_rcsb_simple(protein_name)
        
        if structure_info:
            return {
                "parsed_data": {
                    "protein_targets": [{
                        "name": protein_name,
                        "mutations": mutations
                    }],
                    "drug_requirements": {
                        "binding_affinity": "< -8 kcal/mol",
                        "blood_brain_barrier": "permeable" if "BBB" in query else "not specified"
                    }
                },
                "structures": {protein_name: structure_info},
                "source": "simple_rcsb"
            }
        else:
            raise Exception(f"No RCSB structure found for {protein_name}")
    
    def _extract_protein_name(self, query: str) -> str:
        """Extract protein name using simple keyword matching."""
        query_upper = query.upper()
        
        # Common protein targets
        proteins = ["KRAS", "EGFR", "BRAF", "TP53", "PIK3CA", "PTEN", "MYC", "ACE2"]
        
        for protein in proteins:
            if protein in query_upper:
                print(f"Found protein: {protein}")
                return protein
        
        # Default fallback
        return "KRAS"
    
    def _extract_mutations(self, query: str) -> List[Dict[str, Any]]:
        """Extract mutations using simple regex."""
        import re
        
        mutations = []
        # Look for patterns like G12D, K117N, etc.
        mutation_pattern = r'([A-Z])(\d+)([A-Z])'
        matches = re.findall(mutation_pattern, query)
        
        for match in matches:
            mutations.append({
                "original_residue": match[0],
                "position": int(match[1]), 
                "mutated_residue": match[2],
                "notation": f"{match[0]}{match[1]}{match[2]}"
            })
            
        return mutations
    
    def _search_rcsb_simple(self, protein_name: str) -> Optional[Dict[str, Any]]:
        """Simple RCSB search without LLM."""
        print(f"Searching RCSB for: {protein_name}")
        
        # Use only title search since gene name attributes are not searchable
        search_strategies = [
            # Strategy 1: Title search
            {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_words",
                        "value": protein_name
                    }
                },
                "request_options": {"return_all_hits": False},
                "return_type": "entry"
            }
        ]
        
        for i, search_query in enumerate(search_strategies, 1):
            try:
                print(f"Searching RCSB with title search strategy")
                
                # Search RCSB
                search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
                response = requests.post(search_url, json=search_query, timeout=15)
                
                if response.status_code == 200:
                    data = response.json()
                    results = data.get("result_set", [])
                    
                    if results:
                        print(f"Found {len(results)} results")
                        # Take first result
                        pdb_id = results[0].get("identifier")
                        print(f"Selected PDB: {pdb_id}")
                        
                        # Download structure
                        pdb_file = self._download_pdb(pdb_id)
                        
                        return {
                            "pdb_id": pdb_id,
                            "pdb_file": pdb_file,
                            "source": "simple_rcsb",
                            "organism": "Homo sapiens",
                            "method": "experimental",
                            "resolution": "unknown",
                            "strategy": i
                        }
                    else:
                        print(f"No results found")
                else:
                    print(f"HTTP {response.status_code} - {response.text[:100]}")
                    
            except Exception as e:
                print(f"Search failed: {e}")
                continue
        
        print(f"All search strategies failed for {protein_name}")
        return None
    
    def _download_pdb(self, pdb_id: str) -> str:
        """Download PDB file."""
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        
        try:
            response = requests.get(pdb_url, timeout=15)
            response.raise_for_status()
            
            pdb_file = f"{pdb_id}.pdb"
            with open(pdb_file, "w") as f:
                f.write(response.text)
            
            print(f"Downloaded: {pdb_file}")
            return pdb_file
            
        except Exception as e:
            print(f"Download failed: {e}")
            raise
