"""
Enhanced RCSB PDB Fetcher for F.A.D.E
Dynamic script for retrieving protein structures from RCSB PDB with LLM integration
"""

import requests
import os
import sys
import json
import argparse
from typing import Dict, List, Optional, Any


def search_by_protein_name(protein_name: str, organism: str = "Homo sapiens") -> List[Dict[str, Any]]:
    """
    Search RCSB PDB by protein name and organism.
    
    Args:
        protein_name: Name of the protein to search for
        organism: Organism name (default: Homo sapiens)
        
    Returns:
        List of matching PDB entries
    """
    search_query = {
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
                        "value": protein_name
                    }
                }
            ]
        },
        "request_options": {"return_all_hits": True},
        "return_type": "entry"
    }
    
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    try:
        resp = requests.post(search_url, json=search_query, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        return data.get("result_set", [])
    except Exception as e:
        print(f"Search failed: {e}")
        return []


def search_by_gene_name(gene_name: str, organism: str = "Homo sapiens") -> List[Dict[str, Any]]:
    """
    Search RCSB PDB by gene name and organism.
    
    Args:
        gene_name: Gene name to search for
        organism: Organism name (default: Homo sapiens)
        
    Returns:
        List of matching PDB entries
    """
    search_query = {
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
                        "value": gene_name.upper()
                    }
                }
            ]
        },
        "request_options": {"return_all_hits": True},
        "return_type": "entry"
    }
    
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    try:
        resp = requests.post(search_url, json=search_query, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        return data.get("result_set", [])
    except Exception as e:
        print(f"Gene search failed: {e}")
        return []


def get_structure_metadata(pdb_id: str) -> Optional[Dict[str, Any]]:
    """
    Get detailed metadata for a PDB structure.
    
    Args:
        pdb_id: PDB identifier
        
    Returns:
        Metadata dictionary or None
    """
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(entry_url, timeout=10)
        response.raise_for_status()
        entry_data = response.json()
        
        # Extract key metadata
        metadata = {
            "pdb_id": pdb_id,
            "title": entry_data.get("struct", {}).get("title", ""),
            "experimental_method": entry_data.get("exptl", [{}])[0].get("method", ""),
            "resolution": entry_data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
            "deposition_date": entry_data.get("rcsb_accession_info", {}).get("deposit_date", ""),
            "ligands": []
        }
        
        # Extract ligand information
        entity_info = entry_data.get("rcsb_entry_info", {})
        nonpolymer_comp_ids = entity_info.get("nonpolymer_bound_components", [])
        if nonpolymer_comp_ids:
            # Filter out water and common components
            filtered_ligands = [comp for comp in nonpolymer_comp_ids 
                              if comp not in ["HOH", "SO4", "PO4", "GOL", "EDO", "PEG"]]
            metadata["ligands"] = filtered_ligands
        
        return metadata
        
    except Exception as e:
        print(f"Failed to get metadata for {pdb_id}: {e}")
        return None


def rank_structures(structures: List[Dict[str, Any]], metadata_list: List[Dict[str, Any]]) -> List[tuple]:
    """
    Rank structures by quality for drug discovery.
    
    Args:
        structures: List of structure entries
        metadata_list: List of corresponding metadata
        
    Returns:
        List of (structure, metadata, score) tuples sorted by score (descending)
    """
    ranked = []
    
    for struct, meta in zip(structures, metadata_list):
        if not meta:
            continue
            
        score = 0
        
        # Resolution scoring (lower is better)
        resolution = meta.get("resolution")
        if resolution:
            if resolution <= 1.5:
                score += 100
            elif resolution <= 2.0:
                score += 80
            elif resolution <= 2.5:
                score += 60
            elif resolution <= 3.0:
                score += 40
            else:
                score += 20
        
        # Method scoring
        method = meta.get("experimental_method", "").upper()
        if "X-RAY" in method:
            score += 50
        elif "NMR" in method:
            score += 30
        elif "CRYO" in method:
            score += 40
        
        # Ligand scoring (indicates druggable sites)
        ligands = meta.get("ligands", [])
        if ligands:
            score += 30 + min(len(ligands) * 10, 50)
        
        # Recent structures get bonus
        dep_date = meta.get("deposition_date", "")
        if dep_date and dep_date.startswith("20"):
            year = int(dep_date[:4])
            if year >= 2020:
                score += 20
            elif year >= 2015:
                score += 10
        
        ranked.append((struct, meta, score))
    
    return sorted(ranked, key=lambda x: x[2], reverse=True)


def download_pdb_structure(pdb_id: str, output_dir: str) -> str:
    """
    Download PDB structure file.
    
    Args:
        pdb_id: PDB identifier
        output_dir: Directory to save the file
        
    Returns:
        Path to downloaded file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    try:
        response = requests.get(pdb_url, timeout=30)
        response.raise_for_status()
        
        with open(pdb_path, "w") as f:
            f.write(response.text)
        
        print(f"Downloaded: {pdb_path}")
        return pdb_path
        
    except Exception as e:
        print(f"Download failed for {pdb_id}: {e}")
        raise


def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(description="Fetch protein structures from RCSB PDB")
    parser.add_argument("--protein-name", "-p", help="Protein name to search for")
    parser.add_argument("--gene-name", "-g", help="Gene name to search for") 
    parser.add_argument("--pdb-id", "-i", help="Specific PDB ID to download")
    parser.add_argument("--organism", "-o", default="Homo sapiens", help="Organism name")
    parser.add_argument("--output-dir", "-d", default="./pdb_files", help="Output directory")
    parser.add_argument("--limit", "-l", type=int, default=5, help="Maximum structures to consider")
    parser.add_argument("--json-output", "-j", help="Save results as JSON to this file")
    
    args = parser.parse_args()
    
    # Validate arguments
    if not any([args.protein_name, args.gene_name, args.pdb_id]):
        print("Error: Must specify --protein-name, --gene-name, or --pdb-id")
        sys.exit(1)
    
    try:
        # Direct PDB ID download
        if args.pdb_id:
            print(f"Downloading specific PDB: {args.pdb_id}")
            metadata = get_structure_metadata(args.pdb_id)
            pdb_path = download_pdb_structure(args.pdb_id, args.output_dir)
            
            result = {
                "pdb_id": args.pdb_id,
                "pdb_file": pdb_path,
                "source": "direct_id",
                "metadata": metadata
            }
            
            if args.json_output:
                with open(args.json_output, "w") as f:
                    json.dump(result, f, indent=2)
            
            print(f"Success: {args.pdb_id}")
            return
        
        # Search-based retrieval
        structures = []
        
        if args.gene_name:
            print(f"Searching by gene name: {args.gene_name} in {args.organism}")
            structures = search_by_gene_name(args.gene_name, args.organism)
            search_type = "gene_name"
            search_term = args.gene_name
            
        elif args.protein_name:
            print(f"Searching by protein name: {args.protein_name} in {args.organism}")
            structures = search_by_protein_name(args.protein_name, args.organism)
            search_type = "protein_name"
            search_term = args.protein_name
        
        if not structures:
            print(f"No structures found for {search_term}")
            sys.exit(1)
        
        print(f"Found {len(structures)} structures")
        
        # Get metadata for all structures
        print("Fetching metadata...")
        metadata_list = []
        for struct in structures[:args.limit]:
            pdb_id = struct.get("identifier")
            if pdb_id:
                meta = get_structure_metadata(pdb_id)
                metadata_list.append(meta)
            else:
                metadata_list.append(None)
        
        # Rank structures
        ranked = rank_structures(structures[:args.limit], metadata_list)
        
        if not ranked:
            print("No valid structures with metadata found")
            sys.exit(1)
        
        # Download the best structure
        best_struct, best_meta, best_score = ranked[0]
        pdb_id = best_struct.get("identifier")
        
        print(f"\nSelected best structure:")
        print(f"PDB ID: {pdb_id}")
        print(f"Title: {best_meta.get('title', 'N/A')[:80]}...")
        print(f"Method: {best_meta.get('experimental_method', 'N/A')}")
        print(f"Resolution: {best_meta.get('resolution', 'N/A')} Å")
        print(f"Ligands: {', '.join(best_meta.get('ligands', [])[:3]) or 'None'}")
        print(f"Score: {best_score}")
        
        # Download the structure
        pdb_path = download_pdb_structure(pdb_id, args.output_dir)
        
        # Save results
        result = {
            "pdb_id": pdb_id,
            "pdb_file": pdb_path,
            "source": search_type,
            "search_term": search_term,
            "organism": args.organism,
            "metadata": best_meta,
            "score": best_score,
            "total_found": len(structures),
            "alternatives": [
                {
                    "pdb_id": s.get("identifier"),
                    "score": score,
                    "title": m.get("title", "") if m else ""
                }
                for s, m, score in ranked[1:6]  # Top 5 alternatives
                if m
            ]
        }
        
        # Save JSON output if requested
        if args.json_output:
            with open(args.json_output, "w") as f:
                json.dump(result, f, indent=2)
            print(f"Results saved to: {args.json_output}")
        
        print(f"\nSuccess! Downloaded {pdb_id} to {pdb_path}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
