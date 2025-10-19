#!/usr/bin/env python3

'''
Run a RosettaLigandEnsemble docking protocol between a queried receptor
and top predicted candidates scored by AI-derived SAR:
'''

from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.ligand_docking import TransformEnsemble, HighResEnsemble

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

import os

def heavy_atom_map(prb, ref):
    """Build (prb_idx, ref_idx) mapping for non-H (heavy) atoms.
    Assumes AddHs preserves the ordering of the original heavy atoms in the produced molecule.
    """
    prb_heavy = [a.GetIdx() for a in prb.GetAtoms() if a.GetAtomicNum() != 1]
    ref_heavy = [a.GetIdx() for a in ref.GetAtoms() if a.GetAtomicNum() != 1]
    if len(prb_heavy) != len(ref_heavy):
        raise ValueError("Mismatch in number of heavy atoms between probe and reference")
    return list(zip(prb_heavy, ref_heavy))

def align_conformers_to_ref(prb, ref, conf_ids):
    """Align each conformer of mol_Hs to the reference molecule (ref) using heavy-atom mapping."""
    mapping = heavy_atom_map(prb, ref)
    ref_conf_id = ref.GetConformer().GetId()
    for cid in conf_ids:
        rdMolAlign.AlignMol(prb, ref, prbCid=cid, refCid=ref_conf_id, atomMap=mapping)

def hits_to_conformers(hits_sdf, max_confs_per_hit=10):
    ''' Output conformers for each hit ligand '''
    
    supplier = Chem.SDMolSupplier(hits_sdf)
    outdir = os.path.dirname(hits_sdf)
    hits_tag = os.path.basename(os.path.splitext(hits_sdf)[0])
    pathlist = []

    for i, mol in enumerate(supplier):
        if mol is None:
            continue

        # Track hit structure data:
        hit_name = str(i+1).zfill(3)
        hit_ref = Chem.Mol(mol)

        # Add hydrogens to the hit structure:
        mol = Chem.AddHs(mol, addCoords=True)

        # Generate conformers for each hit:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=max_confs_per_hit,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True
        )

        # Align each conf to hit coords scaffold:
        align_conformers_to_ref(mol, hit_ref, conf_ids=conf_ids)

        out_file = os.path.join(outdir, f'{hits_tag}_{hit_name}_conformers.sdf')
        with Chem.SDWriter(out_file) as writer:
            for conf_id in conf_ids:
                writer.write(mol, confId=conf_id)
        pathlist.append(out_file)

    return pathlist

