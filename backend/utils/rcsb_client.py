'''
RCSB API Client for F.A.D.E.
'''

import os
import json
import logging
from typing import Any, Dict, List, Optional, Union

import requests

class RCSBClient:
    '''
    Client for interacting with the RCSB API to retrieve queried protein information.
    '''

    def __init__(self) -> None:
        '''Initialize an instance of the RCSB client.'''
        self.base_url = 'https://search.rcsb.org/rcsbsearch/'
        self.base_entry_url = 'https://data.rcsb.org/rest/v1/core/entry/'
        self.base_download_url = 'https://files.rcsb.org/download/'
        self.logger = logging.getLogger('fade.rcsb')

    def search_protein(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        '''
        Search RCSB for relevant proteins from an input query.
        '''
        search_url = f'{self.base_url}v2/query'
        params = {
            'query': query,
            'format': 'json',
            'size': limit
        }

        response = requests.post(search_url, params=params)
        response.raise_for_status()
        results = response.json()

        return results.get('results', [])
    
    # NOTE: Natesan is adding further filtration to this
    # process that would result in only one output entry.
    def filter_bound_complexes(self, results) -> List:
        '''
        Filter returned proteins from query keeping only those
        which are bound to a ligand.
        '''
        result_tags = [d['identifier'] for d in results]
        print(f'Found {len(result_tags)} entries in {results}')

        filtered_entries = []
        for pdb_id in result_tags:
            entry_url = f'{self.base_entry_url}{pdb_id}'
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
        
        print(f'Filtered {len(filtered_entries)} entries with ligands')

        return filtered_entries
    
    def fetch_complexes(self, query: str, limit: int = 10, outdir: str = None) -> None:
        '''
        Pull hit complexes from RCSB and save them in pdb format to
        a specified output directory.
        '''
        if not outdir:
            outdir = os.getcwd()

        results = self.search_protein(query=query, limit=limit)
        filtered_entries = self.filter_bound_complexes(results)

        fetched = []
        for entry in filtered_entries:
            pdb_id = entry['pdb_id']
            down_url = f'{self.base_download_url}{pdb_id}.pdb'
            pdb_resp = requests.get(down_url)
            if pdb_resp.status_code == 200:
                with open(os.path.join(outdir, f'{pdb_id}.pdb'), 'w') as f:
                    f.write(pdb_resp.text)
                fetched.append(entry['pdb_id'])
            else:
                print(f'Warning: Failed to download {pdb_id}')
        
        print(f'Wrote {len(fetched)} complexes to {outdir}')

        return fetched
    

