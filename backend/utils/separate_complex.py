#!/usr/bin/env python

'''
Take an input receptor-ligand complex in PDB format and separate into output
cleaned receptor and ligand files in PDB and SDF format, respectively.
'''

import os
import pymol
from pymol import cmd
import argparse
import sys

# Common solvent and ion residue names to exclude
SOLVENT_RESIDUES = {
    'HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'TIP4', 'SPC', 'SOL',  # Water
    'Na+', 'Cl-', 'K+', 'Mg2+', 'Ca2+', 'Zn2+', 'Fe2+', 'Fe3+',  # Ions
    'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CO', 'NI',  # Ions (alternative names)
    'SO4', 'PO4', 'NO3', 'ACE', 'NMA', 'EDO', 'GOL', 'PEG',  # Small molecules/buffers
    'DMS', 'DMSO', 'MPD', 'BME', 'DTT', 'TCEP', 'HEPES', 'TRIS'  # Common crystallization agents
}

# Common cofactors and prosthetic groups (may want to keep or separate depending on use case)
COFACTOR_RESIDUES = {
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'CTP', 'CDP', 'CMP', 'UTP', 'UDP', 'UMP',  # Nucleotides
    'NAD', 'NADH', 'NAP', 'NADP', 'FAD', 'FMN', 'COA', 'HEM', 'HEME',  # Cofactors
    'MG', 'ZN', 'FE', 'CA', 'MN', 'CO', 'NI', 'CU'  # Metal cofactors
}

def filter_ligands_by_size(ligands_dict, min_heavy_atoms=5, max_heavy_atoms=100):
    """
    Filter ligands based on size criteria to exclude very small molecules 
    (likely artifacts) and very large ones (likely misidentified polymers).
    """
    filtered_ligands = {}
    
    for ligand_key, atoms in ligands_dict.items():
        # Count heavy atoms (non-hydrogen)
        heavy_atoms = sum(1 for atom in atoms if atom.symbol != 'H')
        
        if min_heavy_atoms <= heavy_atoms <= max_heavy_atoms:
            filtered_ligands[ligand_key] = atoms
            print(f"  Keeping ligand {ligand_key[0]} (chain {ligand_key[1]}, residue {ligand_key[2]}): {heavy_atoms} heavy atoms")
        else:
            print(f"  Excluding ligand {ligand_key[0]} (chain {ligand_key[1]}, residue {ligand_key[2]}): {heavy_atoms} heavy atoms (outside size range)")
    
    return filtered_ligands

def extract_ligands(pdb_file, complex_name, include_cofactors=False):
    """
    Extract ligands from a protein complex, filtering out solvents and other artifacts.
    
    Parameters:
    - pdb_file: Path to the PDB file
    - complex_name: Name of the loaded structure in PyMOL
    - include_cofactors: Whether to include cofactors as ligands
    
    Returns:
    - Dictionary of ligands with (resn, chain, resi) as keys and atom lists as values
    """
    
    # Define exclusion set
    exclude_residues = SOLVENT_RESIDUES.copy()
    if not include_cofactors:
        exclude_residues.update(COFACTOR_RESIDUES)
    
    # Create exclusion string for PyMOL selection
    exclude_string = ' and '.join([f'not resn {res}' for res in exclude_residues])
    
    # Select potential ligands: organic molecules that are not proteins, nucleic acids, or excluded residues
    selection = f'{complex_name} and organic and not polymer.protein and not polymer.nucleic and {exclude_string}'
    
    print(f"Ligand selection string: {selection}")
    
    # Get model of potential ligands
    try:
        lig_model = cmd.get_model(selection)
    except:
        print("Warning: Could not get ligand model with organic selection. Trying broader selection...")
        # Fallback to broader selection
        selection = f'{complex_name} and not polymer.protein and not polymer.nucleic and {exclude_string}'
        lig_model = cmd.get_model(selection)
    
    # Group atoms by residue
    ligands = {}
    for atom in lig_model.atom:
        key = (atom.resn, atom.chain, atom.resi)
        ligands.setdefault(key, []).append(atom)
    
    print(f"Found {len(ligands)} potential ligand residues before size filtering")
    
    # Filter by size to remove artifacts and misidentified molecules
    filtered_ligands = filter_ligands_by_size(ligands)
    
    print(f"Retained {len(filtered_ligands)} ligands after filtering")
    
    return filtered_ligands

def save_receptor(complex_name, output_path):
    """Save the protein receptor to a PDB file."""
    receptor_selection = f'{complex_name} and polymer.protein'
    
    # Check if receptor exists
    if cmd.count_atoms(receptor_selection) == 0:
        raise ValueError("No protein atoms found in the complex")
    
    cmd.save(output_path, receptor_selection)
    print(f"Receptor saved to: {output_path}")

def save_ligands(complex_name, ligands_dict, output_prefix):
    """Save each ligand to separate SDF files."""
    ligand_files = []
    
    for i, (ligand_key, atoms) in enumerate(ligands_dict.items(), 1):
        resn, chain, resi = ligand_key
        
        # Create selection for this specific ligand
        if chain.strip():  # If chain is specified
            ligand_selection = f'{complex_name} and resn {resn} and chain {chain} and resi {resi}'
            ligand_name = f"{resn}_{chain}_{resi}"
        else:  # If no chain specified
            ligand_selection = f'{complex_name} and resn {resn} and resi {resi}'
            ligand_name = f"{resn}_{resi}"
        
        # Generate output filename
        if len(ligands_dict) == 1:
            output_file = f"{output_prefix}_ligand.sdf"
        else:
            output_file = f"{output_prefix}_ligand_{i}_{ligand_name}.sdf"
        
        # Save ligand as SDF
        try:
            cmd.save(output_file, ligand_selection)
            ligand_files.append(output_file)
            print(f"Ligand {ligand_name} saved to: {output_file}")
        except Exception as e:
            print(f"Error saving ligand {ligand_name}: {e}")
    
    return ligand_files

def main():
    parser = argparse.ArgumentParser(
        description='Utility for automatically separating receptor-ligand complexes into PDB and SDF files respectively.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python separate_complex.py complex.pdb
  python separate_complex.py complex.pdb --outprefix my_structure
  python separate_complex.py complex.pdb --include-cofactors
        '''
    )
    
    parser.add_argument('pdbcomplex', type=str, help='Input receptor-ligand complex in PDB format.')
    parser.add_argument('--outprefix', type=str, default=None, 
                       help='Prefix for output PDB and SDF files (default: based on input filename).')
    parser.add_argument('--include-cofactors', action='store_true',
                       help='Include cofactors and nucleotides as ligands.')
    parser.add_argument('--min-atoms', type=int, default=5,
                       help='Minimum number of heavy atoms for a molecule to be considered a ligand (default: 5).')
    parser.add_argument('--max-atoms', type=int, default=100,
                       help='Maximum number of heavy atoms for a molecule to be considered a ligand (default: 100).')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.pdbcomplex):
        print(f"Error: Input file '{args.pdbcomplex}' not found.")
        sys.exit(1)
    
    # Initialize PyMOL
    pymol.finish_launching(['pymol', '-c'])
    cmd.reinitialize()
    
    # Generate output prefix
    complex_basename = os.path.basename(args.pdbcomplex)
    complex_name = complex_basename.replace('.pdb', '').replace('.PDB', '')
    
    if args.outprefix:
        output_prefix = args.outprefix
    else:
        output_prefix = complex_name
    
    print(f"Processing complex: {args.pdbcomplex}")
    print(f"Output prefix: {output_prefix}")
    
    try:
        # Load the complex
        cmd.load(args.pdbcomplex, complex_name)
        
        # Extract ligands
        ligands = extract_ligands(args.pdbcomplex, complex_name, args.include_cofactors)
        
        if not ligands:
            print("Warning: No ligands found in the complex.")
            # Still save receptor
            receptor_file = f"{output_prefix}_receptor.pdb"
            save_receptor(complex_name, receptor_file)
            return
        
        # Save receptor
        receptor_file = f"{output_prefix}_receptor.pdb"
        save_receptor(complex_name, receptor_file)
        
        # Save ligands
        ligand_files = save_ligands(complex_name, ligands, output_prefix)
        
        print(f"\nSeparation complete!")
        print(f"Receptor: {receptor_file}")
        print(f"Ligands: {', '.join(ligand_files)}")
        
    except Exception as e:
        print(f"Error processing complex: {e}")
        sys.exit(1)
    finally:
        # Clean up PyMOL
        cmd.reinitialize()

if __name__ == "__main__":
    main()
