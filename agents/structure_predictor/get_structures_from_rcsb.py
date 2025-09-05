#!/usr/bin/env python3
"""
Fetch protein structures from RCSB PDB and save them in an AlphaFold 3–like output layout.
Adapted for F.A.D.E integration.
"""

import argparse
import datetime as dt
import json
import os
import sys
import subprocess
from typing import Dict, List, Optional
from pathlib import Path

try:
    import requests
except ImportError:
    sys.stderr.write("This script requires 'requests'. Install with: pip install requests\n")
    raise

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json="
RCSB_SUMMARY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"

# ----------------------------- RCSB search helpers ---------------------------------

def _term_node(attribute: str, value) -> Dict:
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attribute, "operator": "exact_match", "value": str(value)}}

def _exists_node(attribute: str) -> Dict:
    return {"type": "terminal","service": "text","parameters": {"attribute": attribute, "operator": "exists", "value": True}}

def _in_list_node(attribute: str, values: List[str]) -> Dict:
    return {"type": "terminal","service": "text","parameters": {"attribute": attribute, "operator": "in", "value": values}}

def _contains_words(attribute: str, words: str) -> Dict:
    return {"type": "terminal","service":"text","parameters":{"attribute": attribute,"operator":"contains_words","value": words}}

def _group_node(logic: str, nodes: List[Dict]) -> Dict:
    return {"type": "group","logical_operator": logic,"nodes": nodes}

def build_generic_query(
    organism: Optional[str],
    taxid: Optional[int],
    gene: Optional[str],
    uniprot: Optional[str],
    query_text: Optional[str],
    ligands: Optional[List[str]],
    apo: bool,
    max_rows: int = 200
) -> Dict:
    nodes: List[Dict] = []

    # Always protein polymer with coordinates
    nodes.append(_term_node("entity_poly.polymer_type", "Protein"))
    nodes.append(_exists_node("rcsb_entry_info.deposited_polymer_entity_instance_count"))

    # Organism filters
    if taxid is not None:
        nodes.append(_term_node("rcsb_entity_source_organism.taxon_id", str(taxid)))
    if organism:
        nodes.append(_term_node("rcsb_entity_source_organism.taxonomy_lineage.name", organism))

    # Identifiers
    if uniprot:
        nodes.append(_group_node("and", [
            _term_node("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name", "UNIPROT"),
            _term_node("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", uniprot)
        ]))
    if gene:
        nodes.append(_term_node("rcsb_polymer_entity_container_identifiers.gene_name.value", gene))

    # Free-text entity name
    if query_text:
        name_or_desc = _group_node("or", [
            _contains_words("rcsb_polymer_entity_container_identifiers.entity_name", query_text),
            _contains_words("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.gene_name.value", query_text),
            _contains_words("struct.title", query_text)
        ])
        nodes.append(name_or_desc)

    # Ligand constraints
    if apo:
        nodes.append(_exists_node("rcsb_entry_info.nonpolymer_entity_count"))
        nodes.append({"type":"group","logical_operator":"not","nodes":[ _in_list_node("rcsb_nonpolymer_instance_feature_summary.chem_comp_ids", ["GDP","GTP","GNP","GSP","ATP","ADP","AMP","FMN","FAD"]) ]})
    elif ligands:
        nodes.append(_in_list_node("rcsb_nonpolymer_instance_feature_summary.chem_comp_ids", [l.upper() for l in ligands]))

    query = {
        "query": _group_node("and", nodes),
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "scoring_strategy": "combined",
            "sort": [
                {"sort_by": "rcsb_entry_info.experimental_method_priority", "direction": "asc"},
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"},
                {"sort_by": "rcsb_accession_info.initial_release_date", "direction": "desc"}
            ],
            "pager": {"start": 0, "rows": max_rows}
        }
    }
    return query

def rcsb_search_ids(query_json: Dict) -> List[str]:
    resp = requests.get(RCSB_SEARCH_URL + json.dumps(query_json), timeout=30)
    resp.raise_for_status()
    data = resp.json()
    return [item["identifier"] for item in data.get("result_set", [])]

def fetch_entry_summary(pdb_id: str) -> Dict:
    r = requests.get(RCSB_SUMMARY_URL.format(pdb_id), timeout=30)
    r.raise_for_status()
    return r.json()

def try_get(d: Dict, path: List[str], default=None):
    cur = d
    for k in path:
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return default
    return cur

def summarize_entry(entry: Dict) -> Dict:
    pdb_id = try_get(entry, ["rcsb_id"]) or try_get(entry, ["rcsb_id", "entry_id"]) or try_get(entry, ["entry", "id"])
    methods = try_get(entry, ["rcsb_entry_info", "experimental_method"]) or try_get(entry, ["exptl", 0, "method"])
    res = try_get(entry, ["rcsb_entry_info", "resolution_combined"])
    res_val = res[0] if isinstance(res, list) and res else None
    rel_date = try_get(entry, ["rcsb_accession_info", "initial_release_date"])
    ligands = try_get(entry, ["rcsb_entry_info", "nonpolymer_bound_components"])
    return {
        "pdb_id": pdb_id,
        "method": methods,
        "resolution_A": res_val,
        "released": rel_date,
        "ligands": ligands if isinstance(ligands, list) else ([ligands] if ligands else [])
    }

def pick_best(entries: List[Dict]) -> Optional[Dict]:
    def method_rank(m):
        if not m: return 99
        m = m.upper()
        if "X-RAY" in m: return 0
        if "ELECTRON MICROSCOPY" in m or "EM" in m: return 1
        if "NMR" in m: return 2
        return 5

    def parse_date(s):
        try: return dt.datetime.fromisoformat(s.replace("Z",""))
        except Exception: return dt.datetime.min

    if not entries: return None
    return sorted(
        entries,
        key=lambda e: (
            method_rank(e.get("method")),
            e.get("resolution_A") if e.get("resolution_A") is not None else 999.0,
            -parse_date(e.get("released")).timestamp() if e.get("released") else float("inf"),
        )
    )[0]

def download_cif(pdb_id: str, out_path: str):
    url = RCSB_DOWNLOAD_URL.format(pdb_id.upper())
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "wb") as f:
        f.write(r.content)

def write_af3_like_layout(protein_tag: str, output_dir: str, entry: Dict):
    """Write RCSB structure in AF3-like layout for F.A.D.E compatibility"""
    tag_lower = protein_tag.lower()
    base_dir = os.path.join(output_dir, f"{protein_tag}_af3_input")
    prot_dir = os.path.join(base_dir, tag_lower)
    cif_path = os.path.join(prot_dir, f"{tag_lower}_model.cif")
    best_model = os.path.join(output_dir, f"{protein_tag}_best_model.cif")
    summary_json = os.path.join(prot_dir, f"{tag_lower}_summary_confidences.json")
    meta_json = os.path.join(prot_dir, f"{tag_lower}_metadata.json")

    # Download CIF into AF3-style location
    download_cif(entry["pdb_id"], cif_path)

    # Copy to top-level "best model"
    try:
        with open(cif_path, "rb") as src, open(best_model, "wb") as dst:
            dst.write(src.read())
    except Exception as e:
        sys.stderr.write(f"⚠️ Could not write {best_model}: {e}\n")

    # AF3-like summary JSON (F.A.D.E compatible)
    summary_payload = {
        "source": "rcsb_pdb",
        "pdb_id": entry["pdb_id"],
        "experimental_method": entry.get("method"),
        "resolution_A": entry.get("resolution_A"),
        "initial_release_date": entry.get("released"),
        "nonpolymer_bound_components": entry.get("ligands", []),
        "ranking_score": 1.0,  # High score for experimental structures
        "ptm": "N/A",
        "avg_plddt": 95.0,  # High confidence for experimental structures
        "confidence": 1.0
    }
    os.makedirs(os.path.dirname(summary_json), exist_ok=True)
    with open(summary_json, "w") as f:
        json.dump(summary_payload, f, indent=2)

    # Metadata
    with open(meta_json, "w") as f:
        json.dump(entry, f, indent=2)

    return {
        "cif_path": cif_path, 
        "best_model_path": best_model, 
        "summary_json": summary_json, 
        "meta_json": meta_json,
        "entry": entry
    }

def search_and_download(
    uniprot: Optional[str] = None,
    gene: Optional[str] = None,
    organism: Optional[str] = None,
    protein_tag: Optional[str] = None,
    output_dir: str = "./rcsb_out",
    apo: bool = False,
    ligands: Optional[List[str]] = None,
    limit: int = 100
) -> Optional[Dict]:
    """
    Programmatic interface for F.A.D.E integration
    Returns structure info dict if successful, None if no structure found
    """
    
    # Choose protein tag
    protein_tag = (
        protein_tag or
        (uniprot.upper() if uniprot else None) or
        (gene.upper() if gene else None) or
        "PROTEIN"
    )
    
    if not any([uniprot, gene]):
        return None
    
    try:
        # Build and execute search
        query_json = build_generic_query(
            organism=organism,
            taxid=None,
            gene=gene,
            uniprot=uniprot,
            query_text=None,
            ligands=ligands,
            apo=apo,
            max_rows=limit
        )
        
        ids = rcsb_search_ids(query_json)
        if not ids:
            return None
        
        entries = []
        for pid in ids[:limit]:
            try:
                s = fetch_entry_summary(pid)
                entries.append(summarize_entry(s))
            except Exception:
                continue
        
        best_entry = pick_best(entries)
        if not best_entry:
            return None
        
        # Download and write AF3-like layout
        result = write_af3_like_layout(protein_tag, output_dir, best_entry)
        return result
        
    except Exception as e:
        print(f"RCSB search failed: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Fetch protein structures from RCSB and write AlphaFold3-like outputs.")
    parser.add_argument("--uniprot", help="UniProt accession", default=None)
    parser.add_argument("--gene", help="Gene symbol", default=None)
    parser.add_argument("--organism", help="Scientific name", default=None)
    parser.add_argument("--taxid", type=int, help="NCBI Taxonomy ID", default=None)
    parser.add_argument("--query", help="Free-text protein name", default=None)
    parser.add_argument("--pdb-id", help="Specific PDB accession", default=None)
    parser.add_argument("--apo", action="store_true", help="Prefer apo structures")
    parser.add_argument("--ligand", action="append", help="ChemComp ID filter", default=None)
    parser.add_argument("--output-dir", default="./rcsb_af3_out", help="Output directory")
    parser.add_argument("--protein-tag", default=None, help="Tag for filenames")
    parser.add_argument("--limit", type=int, default=100, help="Max entries to consider")

    args = parser.parse_args()

    protein_tag = (
        args.protein_tag or
        (args.uniprot.upper() if args.uniprot else None) or
        (args.gene.upper() if args.gene else None) or
        (args.pdb_id.upper() if args.pdb_id else None) or
        "PROTEIN"
    )

    if args.pdb_id:
        pid = args.pdb_id.upper()
        s = fetch_entry_summary(pid)
        entry = summarize_entry(s)
        paths = write_af3_like_layout(protein_tag, args.output_dir, entry)
    else:
        result = search_and_download(
            uniprot=args.uniprot,
            gene=args.gene,
            organism=args.organism,
            protein_tag=protein_tag,
            output_dir=args.output_dir,
            apo=args.apo,
            ligands=args.ligand,
            limit=args.limit
        )
        
        if not result:
            print("No entries found for the given criteria.")
            return
        
        paths = result

    print("✅ Wrote AF3-like output:")
    for k, v in paths.items():
        if k != "entry":
            print(f"  {k}: {v}")
    print("\nSelected entry metadata:")
    print(json.dumps(paths.get("entry", {}), indent=2))

if __name__ == "__main__":
    main()
