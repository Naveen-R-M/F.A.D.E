#!/usr/bin/env python3
"""
Fetch protein structures from RCSB PDB and save them in an AlphaFold 3–like output layout.
rcsb-api compatible version (uses rcsbapi.search and rcsbapi.data)

Requires:
  - python >= 3.8
  - pip install rcsb-api requests

Example:
  python get_structures_from_rcsb.py \
    --uniprot A0A250XVZ2 \
    --protein-tag KRAS \
    --output-dir structure_test \
    --limit 50
"""

import argparse
import datetime as dt
import json
import os
import sys
from typing import Dict, List, Optional, Iterable

from pathlib import Path

try:
    import requests
except ImportError:
    sys.stderr.write("This script requires 'requests'. Install with: pip install requests\n")
    raise

# --- rcsb-api imports ---
try:
    from rcsbapi.search import AttributeQuery, TextQuery
    from rcsbapi.search import search_attributes as attrs
    from rcsbapi.data import DataQuery as RcsbDataQuery
except ImportError:
    sys.stderr.write("This script now uses 'rcsb-api'. Install with: pip install rcsb-api\n")
    raise


RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"

# ----------------------------- helpers ---------------------------------

def try_get(d: Dict, path: List[str], default=None):
    cur = d
    for k in path:
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return default
    return cur

def _method_rank(m: Optional[str]) -> int:
    if not m:
        return 99
    u = m.upper()
    if "X-RAY" in u:
        return 0
    if "ELECTRON MICROSCOPY" in u or "ELECTRON CRYSTALLOGRAPHY" in u or "EM" in u:
        return 1
    if "NMR" in u:
        return 2
    return 5

def _parse_date(s: Optional[str]) -> dt.datetime:
    if not s:
        return dt.datetime.min
    try:
        return dt.datetime.fromisoformat(s.replace("Z", ""))
    except Exception:
        return dt.datetime.min

# ----------------------------- rcsb-api search & data ---------------------------------

from rcsbapi.search import AttributeQuery, TextQuery, NestedAttributeQuery
from rcsbapi.search import search_attributes as attrs

def build_entity_query(
    organism: Optional[str],
    gene: Optional[str],
    uniprot: Optional[str],
    free_text: Optional[str],
):
    """
    Build a polymer-entity–scoped query using operator syntax and proper nesting.
    Results are polymer entity IDs like '1ABC_1'.
    """
    def AND(a, b):
        return (a & b) if a is not None else b

    q = None  # start empty; add constraints as provided

    # Organism (entity-level)
    if organism:
        q = AND(q, (attrs.rcsb_entity_source_organism.scientific_name == organism))

    # UniProt accession (entity-level) — try two schema styles
    if uniprot:
        u = uniprot.strip().upper()

        # (A) direct list field:
        q_uniprot_a = (attrs.rcsb_polymer_entity_container_identifiers.uniprot_ids == u)

        # (B) nested reference_sequence_identifiers (must use NestedAttributeQuery)
        q_uniprot_b = NestedAttributeQuery(
            AttributeQuery(
                "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                operator="exact_match",
                value="UNIPROT",
            ),
            AttributeQuery(
                "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                operator="exact_match",
                value=u,
            )
        )

        q = AND(q, (q_uniprot_a | q_uniprot_b))

    # Gene symbol (entity-level)
    if gene:
        g = gene.strip().upper()
        q = AND(q, (attrs.rcsb_polymer_entity_container_identifiers.gene_name.value == g))

    # Optional: free-text
    if free_text:
        q = AND(q, TextQuery(value=free_text))

    # If still empty, return a harmless false-like constraint
    if q is None:
        q = (attrs.rcsb_entity_source_organism.scientific_name == "__NO_SUCH_ORG__")

    return q


def execute_entity_query_to_pdb_ids(query, limit: int) -> List[str]:
    """
    Run the entity-scoped query and convert 'entity ids' (e.g., '7P0U_1') to unique PDB IDs ('7P0U').
    """
    # Execute: iterate the query to get result identifiers
    try:
        entity_ids = []
        for i, rid in enumerate(query()):
            if i >= limit:
                break
            entity_ids.append(rid)
    except Exception as e:
        sys.stderr.write(f"RCSB search failed: {e}\n")
        return []

    pdb_ids = sorted({eid.split("_")[0] for eid in entity_ids if isinstance(eid, str) and "_" in eid})
    return pdb_ids


def fetch_entry_metadata_batch(pdb_ids: List[str]) -> List[Dict]:
    """
    Use rcsbapi.data to fetch entry-level metadata for a list of PDB IDs in one shot.
    """
    if not pdb_ids:
        return []

    q = RcsbDataQuery(
        input_type="entries",
        input_ids=pdb_ids,
        return_data_list=[
            "rcsb_id",
            "exptl.method",
            "rcsb_entry_info.resolution_combined",
            "rcsb_accession_info.initial_release_date",
            "rcsb_entry_info.nonpolymer_bound_components",
        ],
    )
    data = q.exec() or {}
    # The JSON shape is: {"data": {"entries": [ {...}, ... ]}}
    entries = try_get(data, ["data", "entries"], default=[]) or []
    return entries


def summarize_entry(entry: Dict) -> Dict:
    pdb_id = try_get(entry, ["rcsb_id"])
    exptl = try_get(entry, ["exptl"], default=[])
    method = None
    if isinstance(exptl, list) and exptl:
        method = try_get(exptl[0], ["method"])
    res = try_get(entry, ["rcsb_entry_info", "resolution_combined"])
    res_val = res[0] if isinstance(res, list) and res else None
    rel_date = try_get(entry, ["rcsb_accession_info", "initial_release_date"])
    ligands = try_get(entry, ["rcsb_entry_info", "nonpolymer_bound_components"], default=[])
    if ligands is None:
        ligands = []
    return {
        "pdb_id": pdb_id,
        "method": method,
        "resolution_A": res_val,
        "released": rel_date,
        "ligands": ligands if isinstance(ligands, list) else ([ligands] if ligands else []),
    }


def filter_entries_by_ligands(
    entries: List[Dict],
    apo: bool = False,
    required_ligands: Optional[List[str]] = None,
) -> List[Dict]:
    """
    Apply ligand / apo logic in Python using the metadata we fetched.
    """
    if not entries:
        return []

    required = {l.upper() for l in (required_ligands or [])}

    # Common nucleotides/analogs to exclude when looking for apo (can adjust)
    exclude_for_apo = {"GDP", "GTP", "GNP", "GSP", "ATP", "ADP", "AMP", "FMN", "FAD"}

    out = []
    for e in entries:
        ligs = set([str(x).upper() for x in (e.get("ligands") or [])])
        if required:
            # keep if any required ligand present
            if ligs & required:
                out.append(e)
        elif apo:
            # keep if no ligand OR (ligands present but none from exclude_for_apo)
            if not ligs or not (ligs & exclude_for_apo):
                out.append(e)
        else:
            out.append(e)
    return out


def pick_best(entries: List[Dict]) -> Optional[Dict]:
    if not entries:
        return None
    return sorted(
        entries,
        key=lambda e: (
            _method_rank(e.get("method")),
            e.get("resolution_A") if e.get("resolution_A") is not None else 999.0,
            -_parse_date(e.get("released")).timestamp() if e.get("released") else float("inf"),
        ),
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
        "confidence": 1.0,
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
        "entry": entry,
    }

# ----------------------------- main orchestration ---------------------------------

def search_and_download(
    uniprot: Optional[str] = None,
    gene: Optional[str] = None,
    organism: Optional[str] = None,
    protein_tag: Optional[str] = None,
    output_dir: str = "./rcsb_out",
    apo: bool = False,
    ligands: Optional[List[str]] = None,
    limit: int = 100,
) -> Optional[Dict]:
    """
    Programmatic interface for F.A.D.E integration using rcsb-api.
    Returns structure info dict if successful, None if no structure found.
    """
    protein_tag = (
        protein_tag
        or (uniprot.upper() if uniprot else None)
        or (gene.upper() if gene else None)
        or "PROTEIN"
    )

    if not any([uniprot, gene]):
        return None

    try:
        # 1) Build & execute entity-scoped query
        q = build_entity_query(
            organism=organism,
            gene=gene,
            uniprot=uniprot,
            free_text=None,
        )
        pdb_ids = execute_entity_query_to_pdb_ids(q, limit=limit)
        if not pdb_ids:
            return None

        # 2) Fetch entry metadata in batch via Data API
        raw_entries = fetch_entry_metadata_batch(pdb_ids)
        summarized = [summarize_entry(e) for e in raw_entries if e]

        # 3) Optional ligand / apo filtering (Python-side)
        filtered = filter_entries_by_ligands(summarized, apo=apo, required_ligands=ligands)

        # 4) Pick best
        best_entry = pick_best(filtered)
        if not best_entry:
            return None

        # 5) Download and write AF3-like layout
        result = write_af3_like_layout(protein_tag, output_dir, best_entry)
        return result

    except Exception as e:
        print(f"RCSB search failed: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Fetch protein structures from RCSB and write AlphaFold3-like outputs (rcsb-api powered).")
    parser.add_argument("--uniprot", help="UniProt accession", default=None)
    parser.add_argument("--gene", help="Gene symbol", default=None)
    parser.add_argument("--organism", help="Scientific name", default=None)
    parser.add_argument("--pdb-id", help="Specific PDB accession", default=None)
    parser.add_argument("--apo", action="store_true", help="Prefer apo structures")
    parser.add_argument("--ligand", action="append", help="ChemComp ID filter (can repeat)", default=None)
    parser.add_argument("--output-dir", default="./rcsb_af3_out", help="Output directory")
    parser.add_argument("--protein-tag", default=None, help="Tag for filenames")
    parser.add_argument("--limit", type=int, default=100, help="Max entries to consider")

    args = parser.parse_args()

    # Protein tag selection priority
    protein_tag = (
        args.protein_tag
        or (args.uniprot.upper() if args.uniprot else None)
        or (args.gene.upper() if args.gene else None)
        or (args.pdb_id.upper() if args.pdb_id else None)
        or "PROTEIN"
    )

    if args.pdb_id:
        # Direct path: fetch metadata via Data API for the provided entry
        pid = args.pdb_id.upper()
        entries = fetch_entry_metadata_batch([pid])
        if not entries:
            print("No metadata found for the given PDB ID.")
            return
        entry = summarize_entry(entries[0])
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
            limit=args.limit,
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
