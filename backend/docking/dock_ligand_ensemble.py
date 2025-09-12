#!/usr/bin/env python

'''
Run a RosettaLigandEnsemble Docking protocol between receptor and
top predicted candidates.
'''

import os
import pyrosetta
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.ligand_docking import HighResDocker, InterfaceScoreCalculator
from pyrosetta.rosetta.core.scoring import get_score_function
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from typing import List, Dict

def generate_conformers_and_params(lig_sdf: str, max_confs: int = 10) -> Dict:
    '''
    Generate rotomeric conformers for each ligand input in sdf
    format up to a maximum number of conformer samples and output
    as individual mol files containing each conformer for each ligand.
    Also generates .params files for each ligand.
    '''
    supplier = Chem.SDMolSupplier(lig_sdf)
    lig_data = {}
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        
        lig_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'ligand_{i+1}'
        
        # Generate conformers
        AllChem.EmbedMultipleConfs(mol, numConfs=max_confs, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        
        conformer_files = []
        for conf_id in range(mol.GetNumConformers()):
            conformer_mol_file = f'{lig_name}_conf_{conf_id+1}.mol'
            Chem.MolToMolFile(mol, conformer_mol_file, confId=conf_id)
            conformer_files.append(conformer_mol_file)

        # Generate params file from the first conformer
        params_file = f'{lig_name}.params'
        molfile_to_params_cmd = [
            'python', os.path.join(os.environ['ROSETTA_HOME'], 'main', 'source', 'scripts', 'python', 'public', 'molfile_to_params.py'),
            '-n', lig_name,
            '-p', params_file,
            '--chain=X',
            conformer_files[0]
        ]
        subprocess.run(molfile_to_params_cmd, check=True)

        lig_data[lig_name] = {
            'conformers': conformer_files,
            'params': params_file
        }
        
    return lig_data

def create_complex_pdb(rec_pdb: str, lig_mol_file: str, complex_pdb_file: str):
    '''
    Combines a receptor PDB and a ligand MOL file into a single PDB file.
    '''
    with open(rec_pdb, 'r') as f_rec, open(lig_mol_file, 'r') as f_lig, open(complex_pdb_file, 'w') as f_out:
        # Write receptor atoms
        for line in f_rec:
            if line.startswith('ATOM'):
                f_out.write(line)
        
        # Write TER record to separate receptor and ligand
        f_out.write("TER\n")

        # Write ligand atoms as HETATM
        lig_mol = Chem.MolFromMolFile(lig_mol_file, removeHs=False)
        conformer = lig_mol.GetConformer()
        for atom in lig_mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            atom_name = atom.GetSymbol()
            res_name = lig_mol.GetProp('_Name') if lig_mol.HasProp('_Name') else "LIG"
            res_name = res_name[:3].upper()
            
            # PDB format for HETATM
            # ATOM      1  N   GLY A   1      27.340  34.410  41.380  1.00  0.00           N
            # HETATM record format
            # Record_name Atom_serial_number Atom_name Residue_name Chain_ID Residue_sequence_number X Y Z Occupancy Temperature_factor Element_symbol
            line = f"HETATM{atom.GetIdx()+1:>5} {atom_name:<4}{res_name:<4}X   1    {pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  1.00  0.00          {atom_name:>2}\n"
            f_out.write(line)
        f_out.write("END\n")


def main():

    parser = argparse.ArgumentParser(description='Run a RosettaLigandEnsemble docking protocol between an input receptor and generated ligand candidates.')
    parser.add_argument('--receptor', type=str, required=True,
                        help='The receptor structure in pdb format for docking.')
    parser.add_argument('--ligands', type=str, required=True,
                        help='The ligand ensemble in combined sdf format to dock with respect to the input receptor.')
    parser.add_argument('--n_samples', type=int, default=100,
                        help='The number of scored docked complexes to output for each input ligand in the ensemble.')
    parser.add_argument('--rosetta_home', type=str, default=os.environ.get('ROSETTA_HOME'), help='Path to your Rosetta installation.')

    args = parser.parse_args()

    if not args.rosetta_home:
        raise ValueError("ROSETTA_HOME environment variable is not set and --rosetta_home argument was not provided.")
    os.environ['ROSETTA_HOME'] = args.rosetta_home

    ligand_data = generate_conformers_and_params(args.ligands)

    for lig_name, data in ligand_data.items():
        params_path = data['params']
        
        for conf_mol_file in data['conformers']:
            complex_pdb = f"{os.path.splitext(args.receptor)[0]}_{lig_name}_{os.path.splitext(conf_mol_file)[0]}.pdb"
            create_complex_pdb(args.receptor, conf_mol_file, complex_pdb)

            init_options = [
                '-in:file:extra_res_fa', params_path,
                '-parser:protocol', 'rosetta_scripts/dock.xml', # Placeholder, we will use pyrosetta movers
                '-nstruct', str(args.n_samples),
                '-out:overwrite'
            ]
            init(' '.join(init_options))

            pose = pose_from_pdb(complex_pdb)
            
            # Setup ScoreFunction
            scorefxn = get_score_function(True)

            # Get ligand residue number
            lig_res_num = 0
            for i in range(1, pose.total_residue() + 1):
                if not pose.residue(i).is_protein():
                    lig_res_num = i
                    break
            if lig_res_num == 0:
                print(f"Warning: No ligand found in {complex_pdb}. Skipping.")
                continue

            # Define movers for protocol:
            high_res_docker = HighResDocker()
            intE_calculator = InterfaceScoreCalculator()
            
            for i in range(args.n_samples):
                docked_pose = pose.clone()
                high_res_docker.apply(docked_pose)
                
                # Score the docked pose
                score = scorefxn(docked_pose)
                intE_calculator.apply(docked_pose)
                
                # Save the docked structure
                output_pdb = f"docked_{lig_name}_{os.path.splitext(conf_mol_file)[0]}_sample_{i+1}.pdb"
                docked_pose.dump_pdb(output_pdb)
                
                print(f"Saved {output_pdb} with score: {score}")

if __name__ == '__main__':
    main()

