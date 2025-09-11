import requests
import os
import csv

# --- User settings ---
protein_name = "EGFR"          # Change to "SARS-CoV-2 Spike Delta" if needed
pdb_save_dir = "./pdb_files"    # Directory to save PDB files
csv_file = "ligand_filtered.csv"

os.makedirs(pdb_save_dir, exist_ok=True)

# --- Step 1: Search protein entries ---
search_query = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "struct.title",
            "operator": "contains_words",
            "value": protein_name
        }
    },
    "request_options": {"return_all_hits": True},
    "return_type": "entry"
}

search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
resp = requests.post(search_url, json=search_query)
resp.raise_for_status()
data = resp.json()

all_pdb_ids = [d["identifier"] for d in data.get("result_set", [])]
print(f"Found {len(all_pdb_ids)} entries for {protein_name}")

# --- Step 2: Filter for bound ligands (<100 atoms if info available) ---
filtered_entries = []

for pdb_id in all_pdb_ids:
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    r = requests.get(entry_url)
    if r.status_code != 200:
        continue
    entry = r.json()
    
    ligand_instances = entry.get("rcsb_nonpolymer_instance_feature_summary", [])
    if ligand_instances:
        # Filter ligands by atom count if available
        small_ligands = [
            (l.get("chem_comp_id"), l.get("atom_count")) 
            for l in ligand_instances 
            if l.get("atom_count") is not None and l.get("atom_count") < 100
        ]
    else:
        # Fallback: include ligand IDs from nonpolymer_bound_components, atom count unknown
        ligands_list = entry.get("rcsb_entry_info", {}).get("nonpolymer_bound_components", [])
        small_ligands = [(l, None) for l in ligands_list]

    if small_ligands:
        filtered_entries.append({
            "pdb_id": pdb_id,
            "ligands": small_ligands
        })

print(f"Filtered {len(filtered_entries)} entries with ligands")

# --- Step 3: Download PDB files ---
for entry in filtered_entries:
    pdb_id = entry["pdb_id"]
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_resp = requests.get(pdb_url)
    if pdb_resp.status_code == 200:
        with open(os.path.join(pdb_save_dir, f"{pdb_id}.pdb"), "w") as f:
            f.write(pdb_resp.text)
    else:
        print(f"Failed to download {pdb_id}")

# --- Step 4: Write CSV ---
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["PDB_ID", "Ligand_ID", "Atom_Count"])
    for entry in filtered_entries:
        pdb_id = entry["pdb_id"]
        for ligand_id, atom_count in entry["ligands"]:
            writer.writerow([pdb_id, ligand_id, atom_count])

print(f"CSV written to {csv_file}")
print(f"PDB files saved in {pdb_save_dir}")

