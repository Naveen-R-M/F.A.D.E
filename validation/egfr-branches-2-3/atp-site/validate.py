"""Map existing fpocket pockets to holo 4G5J; no pocket detection or model reruns."""
import json, re, hashlib
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parent
REPO = Path('/scratch/rajagopalmohanraj.n/F.A.D.E')
parser = PDBParser(QUIET=True)
def load(path): return parser.get_structure(Path(path).stem, str(path))
def protein(chain): return [r for r in chain if r.id[0] == ' ' and 'CA' in r]
def sequence(residues): return ''.join(seq1(r.resname, custom_map={'MSE':'M'}) for r in residues)
def heavy(residue): return [a for a in residue if a.element not in ('H','D')]
def distances(a,b): return np.linalg.norm(np.asarray(a)[:,None,:]-np.asarray(b)[None,:,:],axis=-1)
def rid(r): return f'{r.parent.id}:{r.id[1]}{r.id[2].strip()}:{r.resname}'
def key(r): return (r.parent.id, r.id[1], r.id[2])
def save(structure, path):
 io=PDBIO();io.set_structure(structure);io.save(str(path))

aligner=PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
holo=load(ROOT/'4G5J.pdb')
ligands=[r for r in holo[0].get_residues() if r.resname=='0WN']
assert ligands, 'No 0WN in reference'
# Use the protein chain closest to a reference 0WN. Record the choice explicitly.
choices=[]
for ligand in ligands:
 for chain in holo[0]:
  residues=protein(chain)
  if len(residues)<100: continue
  d=distances([a.coord for r in residues for a in heavy(r)], [a.coord for a in heavy(ligand)]).min()
  choices.append((float(d),chain.id,ligand))
_,chainid,ligand=min(choices,key=lambda x:(x[0],x[1]))
reference=protein(holo[0][chainid]); ligcoords=np.array([a.coord for a in heavy(ligand)])
site=[r for r in reference if distances([a.coord for a in heavy(r)],ligcoords).min()<=5.0]
sitekeys={key(r) for r in site}
assert len(site)>=5
sources=[('branch2_7SI1', REPO/'fpocket_inputs/36906525/7SI1.pdb', REPO/'fpocket_inputs/36906525/7SI1_out'),
 ('branch3_predicted',REPO/'backend/data/structures/boltz2_EGFR_boltz2_1792519941.pdb',REPO/'fpocket_inputs/53ed3c6b/boltz2_EGFR_boltz2_1792519941_out')]
result={'method':{'reference':'4G5J','reference_chain':chainid,'ligand':rid(ligand),
 'site_definition':'Reference protein residues with any heavy atom within 5 A of 0WN',
 'reference_site_residues':[rid(r) for r in site],
 'alignment':'Global sequence alignment, then least-squares rigid CA fit of all aligned residues; no outlier removal',
 'thresholds':'Descriptive overlap only; minimum distance alone does not establish ATP-site identity'},'structures':{}}
manifest={}
def record(path): manifest[str(path)]=hashlib.sha256(Path(path).read_bytes()).hexdigest()
record(ROOT/'4G5J.pdb')
for label,path,out in sources:
 record(path)
 moving=load(path)
 candidates=[]
 for chain in moving[0]:
  residues=protein(chain)
  if not residues: continue
  alignment=aligner.align(sequence(reference),sequence(residues))[0]
  candidates.append((alignment.score,residues,alignment))
 _,residues,alignment=max(candidates,key=lambda x:x[0])
 pairs=[(reference[i],residues[j]) for aa,bb in zip(*alignment.aligned) for i,j in zip(range(*aa),range(*bb))]
 mapping={key(m):r for r,m in pairs}
 sup=Superimposer();sup.set_atoms([r['CA'] for r,m in pairs],[m['CA'] for r,m in pairs])
 # Validate direction and preservation of rigid geometry before using transform.
 rot,tran=sup.rotran
 before=np.array([m['CA'].coord for r,m in pairs]);fixed=np.array([r['CA'].coord for r,m in pairs])
 transformed=before@rot+tran
 assert np.isclose(np.sqrt(np.mean(np.sum((transformed-fixed)**2,axis=1))),sup.rms)
 assert np.allclose(np.linalg.norm(before[1:]-before[:-1],axis=1),np.linalg.norm(transformed[1:]-transformed[:-1],axis=1))
 sup.apply(list(moving.get_atoms()));save(moving,ROOT/f'{label}_aligned.pdb')
 local=[(r,m) for r,m in pairs if key(r) in sitekeys]
 site_rms=float(np.sqrt(np.mean([np.sum((r['CA'].coord-m['CA'].coord)**2) for r,m in local])))
 # Also report local-fit RMSD separately: do NOT use it to score pockets.
 localfit=Superimposer();localfit.set_atoms([r['CA'] for r,m in local],[m['CA'] for r,m in local])
 entry={'input':str(path),'moving_chain':residues[0].parent.id,'aligned_CA_count':len(pairs),
 'sequence_identity':sum(r.resname==m.resname for r,m in pairs)/len(pairs),
 'global_CA_RMSD_A':float(sup.rms),'mapped_site_residue_count':len(local),
 'site_CA_RMSD_after_global_fit_A':site_rms,'site_CA_RMSD_after_local_fit_A':float(localfit.rms),
 'rotation_row_vector':rot.tolist(),'translation_A':tran.tolist(),'pockets':[]}
 info=next(out.glob('*info.txt'));record(info)
 scores={int(n):float(s) for n,s in re.findall(r'Pocket\s+(\d+)\s*:.*?Druggability Score\s*:\s*([\d.]+)',info.read_text(),re.S)}
 for p in sorted((out/'pockets').glob('pocket*_atm.pdb'), key=lambda p:int(re.search(r'pocket(\d+)',p.name)[1])):
  record(p);n=int(re.search(r'pocket(\d+)',p.name)[1]);pocket=load(p)
  prs=[r for r in pocket.get_residues() if r.id[0]==' ']
  mapped={key(mapping[key(r)]) for r in prs if key(r) in mapping}
  overlap=mapped & sitekeys
  pc=np.array([a.coord for r in prs for a in heavy(r)])@rot+tran
  ds=distances(pc,ligcoords)
  vpath=p.with_name(p.name.replace('_atm.pdb','_vert.pqr'))
  spheres=[]
  if vpath.exists():
   record(vpath)
   for line in vpath.read_text().splitlines():
    if line.startswith(('ATOM  ','HETATM')):
     spheres.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
  spherecoverage=None
  if spheres:
   vd=distances(np.asarray(spheres)@rot+tran,ligcoords)
   spherecoverage=float(np.mean(vd.min(axis=0)<=3.0))
  row={'pocket':n,'druggability':scores[n],'pocket_residue_count':len({key(r) for r in prs}),
   'mapped_residue_count':len(mapped),'overlap_count':len(overlap),
   'overlap_reference_residues':[rid(r) for r in site if key(r) in overlap],
   'site_recall':len(overlap)/len(sitekeys),'pocket_precision':len(overlap)/len(mapped) if mapped else 0,
   'residue_Jaccard':len(overlap)/len(mapped|sitekeys),'min_pocket_atom_to_0WN_A':float(ds.min()),
   'ligand_heavy_atom_fraction_within_5A_of_pocket_atoms':float(np.mean(ds.min(axis=0)<=5)),
   'ligand_heavy_atom_fraction_within_3A_of_alpha_sphere_centers':spherecoverage}
  entry['pockets'].append(row)
  # Save selected pockets in the same frame as 4G5J for visual inspection.
  if n in (11,15):
   sup.apply(list(pocket.get_atoms()));save(pocket,ROOT/f'{label}_pocket{n}_aligned.pdb')
 result['structures'][label]=entry
(ROOT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'input_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
lines=['# EGFR ATP-site validation','',f'Reference: **4G5J chain {chainid}, ligand {rid(ligand)}**. Site: {len(site)} protein residues within 5 Å of a ligand heavy atom.','',
 'Existing fpocket outputs only; no pocket detection, prediction, threshold changes, or molecule generation. Global sequence alignment establishes residue correspondence; all paired Cα atoms determine a rigid fit. Distances use this global fit. Pocket numbers are local identifiers, not equivalent sites across structures.','',
 '**Interpretation:** residue overlap and cavity coverage must be considered together. Minimum atom distance alone is not evidence of ATP-site identity. Alpha-sphere coverage below uses a descriptive 3 Å distance to sphere centers, not exact cavity-volume overlap.','']
for label,e in result['structures'].items():
 lines += [f'## {label}','',f"Aligned Cα pairs: {e['aligned_CA_count']}; identity: {e['sequence_identity']:.1%}; global RMSD: {e['global_CA_RMSD_A']:.3f} Å. Site Cα RMSD after global fit: {e['site_CA_RMSD_after_global_fit_A']:.3f} Å; after site-only fit: {e['site_CA_RMSD_after_local_fit_A']:.3f} Å. Site residues mapped: {e['mapped_site_residue_count']}/{len(site)}.",'',
 '| Pocket | Drug score | Site overlap | Site recall | Pocket precision | Min distance (Å) | Ligand coverage: pocket atoms ≤5 Å | Ligand coverage: sphere centers ≤3 Å |',
 '|---|---:|---:|---:|---:|---:|---:|---:|']
 for r in e['pockets']:
  sc=r['ligand_heavy_atom_fraction_within_3A_of_alpha_sphere_centers']
  lines.append(f"| {r['pocket']} | {r['druggability']:.3f} | {r['overlap_count']}/{len(site)} | {r['site_recall']:.1%} | {r['pocket_precision']:.1%} | {r['min_pocket_atom_to_0WN_A']:.2f} | {r['ligand_heavy_atom_fraction_within_5A_of_pocket_atoms']:.1%} | {sc:.1%} |" if sc is not None else str(r))
 lines+=['']
(ROOT/'metrics.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines));print('VALIDATION_PASS: transform assertions and output generation complete')
