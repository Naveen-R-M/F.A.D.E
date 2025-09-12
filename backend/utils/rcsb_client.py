"""
RCSB PDB Client for F.A.D.E Target Selector
Enhanced client for retrieving protein structures from RCSB PDB with LLM-powered search
"""

import os
import json
import requests
import sys
from typing import Dict, List, Optional, Any, Tuple
from utils.logging import get_logger


class RCSBClient:
    """
    Enhanced RCSB PDB client that integrates with F.A.D.E target selection workflow.
    Uses LLM-generated search queries to find appropriate protein structures.
    """
    
    def __init__(self):
        """Initialize the RCSB client."""
        self.base_search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.base_data_url = "https://data.rcsb.org/rest/v1/core/entry"
        self.base_download_url = "https://files.rcsb.org/download"
        self.logger = get_logger("fade.rcsb_client")
        
    def extract_pdb_id_from_query(self, llm_client, query: str) -> Optional[str]:
        """
        Use LLM to extract PDB ID or generate search terms from user query.
        
        Args:
            llm_client: LLM client (Gemini) for processing the query
            query: User's natural language query
            
        Returns:
            PDB ID if found, or None to proceed with protein name search
        """
        try:
            # Prompt to extract PDB ID from query
            extraction_prompt = f"""
            Analyze the following drug discovery query and extract any specific PDB ID mentioned:
            
            Query: "{query}"
            
            Look for:
            1. Direct PDB ID mentions (4-character codes like "1ABC", "2XYZ")
            2. Indirect references to specific structures
            
            If a specific PDB ID is mentioned, return just the PDB ID (like "1ABC").
            If no PDB ID is mentioned, return "NONE".
            
            Response format: Just the PDB ID or "NONE"
            """
            
            response = llm_client.generate_content(extraction_prompt)
            pdb_id = response.strip().upper()
            
            if pdb_id != "NONE" and len(pdb_id) == 4 and pdb_id.isalnum():
                self.logger.info(f"Extracted PDB ID from query: {pdb_id}")
                return pdb_id
            else:
                self.logger.info("No specific PDB ID found in query")
                return None
                
        except Exception as e:
            self.logger.error(f"Error extracting PDB ID from query: {e}")
            return None
    
    def generate_search_terms(self, llm_client, target_info: Dict[str, Any]) -> Dict[str, Any]:
        """
        Use LLM to generate optimal search terms for RCSB PDB search.
        
        Args:
            llm_client: LLM client for generating search terms
            target_info: Target information from query parsing
            
        Returns:
            Dictionary with search strategy and terms
        """
        try:
            # Create prompt for search term generation
            prompt = f"""
            Generate optimal search terms for finding protein structures in the RCSB Protein Data Bank.
            
            Target Information:
            - Protein targets: {target_info.get('protein_targets', [])}
            - Disease context: {target_info.get('disease_context', 'not specified')}
            - Mutations: {target_info.get('mutations', [])}
            - Requirements: {target_info.get('requirements', {})}
            
            Generate search terms for the most important protein target. Consider:
            1. Official protein name
            2. Gene symbol
            3. Common aliases
            4. Organism (prefer human/mouse structures)
            5. Specific mutations if mentioned
            
            Return a JSON object with:
            {{
                "primary_name": "official protein name",
                "gene_symbol": "gene symbol", 
                "aliases": ["alias1", "alias2"],
                "organism": "organism name",
                "search_priority": ["term1", "term2", "term3"],
                "mutations": ["mutation1", "mutation2"]
            }}
            
            Focus on the most druggable and well-studied target.
            """
            
            response = llm_client.generate_content(prompt)
            
            # Parse JSON response
            try:
                search_terms = json.loads(response)
                self.logger.info(f"Generated search terms for {search_terms.get('primary_name', 'unknown')}")
                return search_terms
            except json.JSONDecodeError:
                # Fallback to extracting from target_info directly
                targets = target_info.get('protein_targets', [])
                if targets:
                    primary_target = targets[0]
                    return {
                        "primary_name": primary_target.get('name', 'unknown'),
                        "gene_symbol": primary_target.get('name', 'unknown'),
                        "aliases": [],
                        "organism": "Homo sapiens",
                        "search_priority": [primary_target.get('name', 'unknown')],
                        "mutations": primary_target.get('mutations', [])
                    }
                else:
                    return None
                    
        except Exception as e:
            self.logger.error(f"Error generating search terms: {e}")
            return None
    
    def search_structures(self, search_terms: Dict[str, Any], limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search RCSB PDB for structures matching the search terms.
        
        Args:
            search_terms: Search terms generated by LLM
            limit: Maximum number of results to return
            
        Returns:
            List of structure entries from RCSB
        """
        results = []
        
        # Try different search strategies in order of preference
        search_strategies = self._build_search_strategies(search_terms)
        
        for strategy in search_strategies:
            try:
                self.logger.info(f"Trying search strategy: {strategy['description']}")
                
                search_query = {
                    "query": strategy["query"],
                    "request_options": {"return_all_hits": False},
                    "return_type": "entry"
                }
                
                response = requests.post(self.base_search_url, json=search_query, timeout=30)
                response.raise_for_status()
                
                data = response.json()
                search_results = data.get("result_set", [])
                
                if search_results:
                    self.logger.info(f"Found {len(search_results)} results with {strategy['description']}")
                    results.extend(search_results[:limit])
                    break  # Use first successful strategy
                    
            except Exception as e:
                self.logger.warning(f"Search strategy failed: {strategy['description']} - {e}")
                continue
        
        # Enhance results with metadata
        enhanced_results = []
        for result in results[:limit]:
            pdb_id = result.get("identifier")
            if pdb_id:
                metadata = self._get_structure_metadata(pdb_id)
                if metadata:
                    enhanced_result = {**result, **metadata}
                    enhanced_results.append(enhanced_result)
        
        return enhanced_results
    
    def _build_search_strategies(self, search_terms: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Build different search strategies to try in order.
        
        Args:
            search_terms: Search terms from LLM
            
        Returns:
            List of search strategies
        """
        strategies = []
        
        primary_name = search_terms.get("primary_name", "")
        gene_symbol = search_terms.get("gene_symbol", "")
        organism = search_terms.get("organism", "Homo sapiens")
        
        # Strategy 1: Gene symbol + organism (most specific)
        if gene_symbol:
            strategies.append({
                "description": f"Gene symbol '{gene_symbol}' in {organism}",
                "query": {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                                "operator": "exact_match",
                                "value": organism
                            }
                        },
                        {
                            "type": "terminal", 
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_polymer_entity.rcsb_gene_name.value",
                                "operator": "exact_match",
                                "value": gene_symbol.upper()
                            }
                        }
                    ]
                }
            })
        
        # Strategy 2: Protein name + organism
        if primary_name:
            strategies.append({
                "description": f"Protein name '{primary_name}' in {organism}",
                "query": {
                    "type": "group",
                    "logical_operator": "and", 
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                                "operator": "exact_match",
                                "value": organism
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text", 
                            "parameters": {
                                "attribute": "struct.title",
                                "operator": "contains_words",
                                "value": primary_name
                            }
                        }
                    ]
                }
            })
        
        # Strategy 3: Gene symbol only (broader)
        if gene_symbol:
            strategies.append({
                "description": f"Gene symbol '{gene_symbol}' (any organism)",
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity.rcsb_gene_name.value", 
                        "operator": "exact_match",
                        "value": gene_symbol.upper()
                    }
                }
            })
        
        # Strategy 4: Protein name only (broadest)
        if primary_name:
            strategies.append({
                "description": f"Protein name '{primary_name}' (any organism)",
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_words", 
                        "value": primary_name
                    }
                }
            })
        
        return strategies
    
    def _get_structure_metadata(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Get metadata for a specific PDB structure.
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            Metadata dictionary or None
        """
        try:
            entry_url = f"{self.base_data_url}/{pdb_id}"
            response = requests.get(entry_url, timeout=10)
            response.raise_for_status()
            
            entry_data = response.json()
            
            # Extract relevant metadata
            metadata = {
                "pdb_id": pdb_id,
                "title": entry_data.get("struct", {}).get("title", ""),
                "experimental_method": entry_data.get("exptl", [{}])[0].get("method", ""),
                "resolution": entry_data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                "deposition_date": entry_data.get("rcsb_accession_info", {}).get("deposit_date", ""),
                "organism": None,
                "ligands": []
            }
            
            # Extract organism info
            entity_info = entry_data.get("rcsb_entry_info", {})
            polymer_entity_count = entity_info.get("polymer_entity_count", 0)
            
            if polymer_entity_count > 0:
                # Try to get organism from entity source organism
                source_organisms = entry_data.get("rcsb_entity_source_organism", [])
                if source_organisms:
                    metadata["organism"] = source_organisms[0].get("taxonomy_lineage", [{}])[-1].get("name", "")
            
            # Extract ligand information
            nonpolymer_comp_ids = entity_info.get("nonpolymer_bound_components", [])
            if nonpolymer_comp_ids:
                # Filter out water and common buffer components
                filtered_ligands = [comp for comp in nonpolymer_comp_ids 
                                  if comp not in ["HOH", "SO4", "PO4", "GOL", "EDO", "PEG"]]
                metadata["ligands"] = filtered_ligands[:3]  # Top 3 ligands
            
            return metadata
            
        except Exception as e:
            self.logger.warning(f"Failed to get metadata for {pdb_id}: {e}")
            return None
    
    def select_best_structure(self, llm_client, structures: List[Dict[str, Any]], 
                            target_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Use LLM to select the best structure from search results.
        
        Args:
            llm_client: LLM client for decision making
            structures: List of structure candidates
            target_info: Original target information
            
        Returns:
            Best structure or None
        """
        if not structures:
            return None
            
        if len(structures) == 1:
            return structures[0]
        
        try:
            # Create structure summaries for LLM
            structure_summaries = []
            for i, struct in enumerate(structures[:5]):  # Limit to top 5 for LLM processing
                summary = {
                    "index": i,
                    "pdb_id": struct.get("pdb_id", "unknown"),
                    "title": struct.get("title", "")[:100],  # Truncate long titles
                    "method": struct.get("experimental_method", "unknown"),
                    "resolution": struct.get("resolution"),
                    "organism": struct.get("organism", "unknown"),
                    "ligands": struct.get("ligands", [])
                }
                structure_summaries.append(summary)
            
            # Create selection prompt
            prompt = f"""
            Select the best protein structure for drug discovery from these candidates:
            
            Target Information:
            - Main target: {target_info.get('protein_targets', [{}])[0].get('name', 'unknown')}
            - Disease: {target_info.get('disease_context', 'not specified')}
            - Requirements: {target_info.get('requirements', {})}
            
            Structure Candidates:
            {json.dumps(structure_summaries, indent=2)}
            
            Selection Criteria (in order of importance):
            1. Human protein (Homo sapiens) preferred
            2. High resolution (< 2.5 Å preferred)
            3. X-ray crystallography preferred
            4. Contains bound ligands (indicates druggable sites)
            5. Recent structure (newer is better)
            6. Clear, descriptive title matching target
            
            Return only the index number (0-{len(structure_summaries)-1}) of the best structure.
            """
            
            response = llm_client.generate_content(prompt)
            
            try:
                selected_index = int(response.strip())
                if 0 <= selected_index < len(structures):
                    selected = structures[selected_index]
                    self.logger.info(f"LLM selected structure: {selected.get('pdb_id', 'unknown')}")
                    return selected
            except ValueError:
                pass
            
            # Fallback: return first structure with good resolution
            for struct in structures:
                resolution = struct.get("resolution")
                if resolution and resolution < 3.0:
                    self.logger.info(f"Fallback selected structure: {struct.get('pdb_id', 'unknown')}")
                    return struct
            
            # Last fallback: return first structure
            return structures[0]
            
        except Exception as e:
            self.logger.error(f"Error in structure selection: {e}")
            return structures[0] if structures else None
    
    def download_structure(self, pdb_id: str, output_dir: str) -> str:
        """
        Download PDB structure file.
        
        Args:
            pdb_id: PDB identifier
            output_dir: Directory to save the file
            
        Returns:
            Path to downloaded PDB file
        """
        os.makedirs(output_dir, exist_ok=True)
        pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
        
        try:
            pdb_url = f"{self.base_download_url}/{pdb_id}.pdb"
            response = requests.get(pdb_url, timeout=30)
            response.raise_for_status()
            
            with open(pdb_path, "w") as f:
                f.write(response.text)
            
            self.logger.info(f"Downloaded PDB structure: {pdb_path}")
            return pdb_path
            
        except Exception as e:
            self.logger.error(f"Failed to download {pdb_id}: {e}")
            raise
    
    def process_target_query(self, llm_client, query: str, target_info: Dict[str, Any], 
                           output_dir: str) -> Optional[Dict[str, Any]]:
        """
        Complete process: extract/search for structures and download the best one.
        
        Args:
            llm_client: LLM client
            query: Original user query
            target_info: Parsed target information
            output_dir: Directory for outputs
            
        Returns:
            Structure information dictionary
        """
        try:
            # Step 1: Check if specific PDB ID mentioned in query
            specific_pdb = self.extract_pdb_id_from_query(llm_client, query)
            
            if specific_pdb:
                # Direct PDB download
                pdb_path = self.download_structure(specific_pdb, output_dir)
                metadata = self._get_structure_metadata(specific_pdb)
                
                return {
                    "pdb_id": specific_pdb,
                    "pdb_file": pdb_path,
                    "source": "rcsb_direct",
                    "title": metadata.get("title", "") if metadata else "",
                    "method": metadata.get("experimental_method", "") if metadata else "",
                    "resolution": metadata.get("resolution") if metadata else None,
                    "organism": metadata.get("organism", "") if metadata else "",
                    "ligands": metadata.get("ligands", []) if metadata else []
                }
            
            # Step 2: Generate search terms using LLM
            search_terms = self.generate_search_terms(llm_client, target_info)
            if not search_terms:
                self.logger.error("Failed to generate search terms")
                return None
            
            # Step 3: Search for structures
            structures = self.search_structures(search_terms, limit=10)
            if not structures:
                self.logger.error("No structures found in RCSB")
                return None
            
            # Step 4: Select best structure using LLM
            best_structure = self.select_best_structure(llm_client, structures, target_info)
            if not best_structure:
                self.logger.error("Failed to select best structure")
                return None
            
            # Step 5: Download the selected structure
            pdb_id = best_structure.get("pdb_id")
            pdb_path = self.download_structure(pdb_id, output_dir)
            
            # Step 6: Return comprehensive structure information
            return {
                "pdb_id": pdb_id,
                "pdb_file": pdb_path,
                "source": "rcsb_search",
                "title": best_structure.get("title", ""),
                "method": best_structure.get("experimental_method", ""),
                "resolution": best_structure.get("resolution"),
                "organism": best_structure.get("organism", ""),
                "ligands": best_structure.get("ligands", []),
                "search_terms": search_terms,
                "total_candidates": len(structures)
            }
            
        except Exception as e:
            self.logger.error(f"Error processing target query: {e}")
            return None
