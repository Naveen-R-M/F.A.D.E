"""Run the production graph with explicit EGFR inputs; no research-model dependency."""
import asyncio
import json
import sys
from pathlib import Path
sys.path.insert(0, '/scratch/rajagopalmohanraj.n/F.A.D.E/backend')
from fade.workflows.drug_discovery import create_drug_discovery_graph, run_drug_discovery_pipeline

ROOT = Path(__file__).parent
branch = int(sys.argv[1])
assert branch in (2, 3)
entity = json.loads((ROOT / '7si1-entity1.json').read_text())
sequence = entity['entity_poly']['pdbx_seq_one_letter_code_can']
target = dict(protein_name='EGFR kinase domain (7SI1 chain A construct)', gene_name='EGFR',
              uniprot_id='P00533', sequence=sequence, sequence_length=len(sequence))
if branch == 2:
    target['existing_structures'] = [dict(pdb_id='7SI1', resolution=1.6,
        title='Crystal structure of apo EGFR kinase domain', has_drug_like_ligand=False,
        ligands=[], drug_like_ligands=[])]

async def run():
    app = await create_drug_discovery_graph()
    try:
        return await run_drug_discovery_pipeline(app,
            f'Generate candidate ligands for EGFR kinase domain, branch {branch}',
            target_info=target, force_prediction=(branch == 3), pocket_detection_mode='fpocket')
    finally:
        await app.checkpointer.conn.close()

state = asyncio.run(run())
(ROOT / f'branch{branch}-state.json').write_text(json.dumps(state, indent=2, default=str))
assert not state.get('error'), state.get('error')
assert state.get('final_report'), 'Missing final report'
assert state['structure_source'].lower() == ('pdb' if branch == 2 else 'boltz2')
if branch == 2:
    assert state['structure']['pdb_id'] == '7SI1'
print('EGFR BRANCH', branch, 'COMPLETE', state.get('run_id'))
print(state['final_report'])
