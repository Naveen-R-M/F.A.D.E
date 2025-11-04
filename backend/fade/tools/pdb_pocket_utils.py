"""
Utility functions for extracting pocket information from PDB structures.
"""

import re
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

from fade.utils import get_logger

logger = get_logger("tools.pdb_pocket_utils")


def extract_ligand_pocket(
    pdb_path: str,
    ligand_id: str,
    padding: float = 5.0
) -> Dict[str, Any]:
    """
    Extract pocket information from a bound ligand in a PDB structure.
    
    This function:
    1. Finds the ligand in the PDB file
    2. Extracts its center coordinates
    3. Defines a pocket region around the ligand
    
    Args:
        pdb_path: Path to PDB file
        ligand_id: 3-letter code of the ligand (e.g., "ANP", "AEE")
        padding: Distance in Angstroms to expand pocket beyond ligand
        
    Returns:
        Dictionary containing pocket information:
        - center: [x, y, z] coordinates of pocket center
        - residues: List of nearby residues
        - ligand_id: The ligand that defined this pocket
        - pocket_type: "ligand_defined"
    """
    logger.info(f"Extracting pocket from ligand {ligand_id} in {pdb_path}")
    
    if not Path(pdb_path).exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    
    ligand_atoms = []
    protein_atoms = []
    
    # Read PDB file and extract atoms
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Protein atom
                atom_info = parse_pdb_atom_line(line)
                if atom_info:
                    protein_atoms.append(atom_info)
                    
            elif line.startswith('HETATM'):
                # Heteroatom (could be ligand)
                atom_info = parse_pdb_atom_line(line)
                if atom_info and atom_info['residue_name'].strip() == ligand_id:
                    ligand_atoms.append(atom_info)
    
    if not ligand_atoms:
        raise ValueError(f"Ligand {ligand_id} not found in PDB file")
    
    logger.info(f"Found {len(ligand_atoms)} atoms for ligand {ligand_id}")
    
    # Calculate ligand center
    center_x = sum(atom['x'] for atom in ligand_atoms) / len(ligand_atoms)
    center_y = sum(atom['y'] for atom in ligand_atoms) / len(ligand_atoms)
    center_z = sum(atom['z'] for atom in ligand_atoms) / len(ligand_atoms)
    
    logger.info(f"Ligand center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
    
    # Find protein residues near the ligand
    nearby_residues = set()
    cutoff_distance = 8.0  # Angstroms - residues within this distance are part of the pocket
    
    for p_atom in protein_atoms:
        distance = calculate_distance(
            (p_atom['x'], p_atom['y'], p_atom['z']),
            (center_x, center_y, center_z)
        )
        
        if distance <= cutoff_distance:
            residue_id = f"{p_atom['residue_name']}{p_atom['residue_number']}{p_atom['chain_id']}"
            nearby_residues.add(residue_id)
    
    logger.info(f"Found {len(nearby_residues)} residues in pocket")
    
    # Create pocket dictionary
    pocket = {
        "center": [center_x, center_y, center_z],
        "residues": sorted(list(nearby_residues)),
        "ligand_id": ligand_id,
        "pocket_type": "ligand_defined",
        "pocket_score": 1.0,  # Highest confidence since it's from actual ligand
        "pocket_volume": estimate_pocket_volume(ligand_atoms),
        "description": f"Pocket defined by bound ligand {ligand_id}"
    }
    
    return pocket


def parse_pdb_atom_line(line: str) -> Optional[Dict[str, Any]]:
    """
    Parse a single ATOM or HETATM line from a PDB file.
    
    Args:
        line: PDB atom line
        
    Returns:
        Dictionary with atom information or None if parsing fails
    """
    try:
        # Standard PDB format columns
        record_type = line[0:6].strip()
        atom_number = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_number = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        
        return {
            'record_type': record_type,
            'atom_number': atom_number,
            'atom_name': atom_name,
            'residue_name': residue_name,
            'chain_id': chain_id,
            'residue_number': residue_number,
            'x': x,
            'y': y,
            'z': z
        }
    except (ValueError, IndexError):
        return None


def calculate_distance(coord1: Tuple[float, float, float], 
                       coord2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two 3D coordinates."""
    return ((coord1[0] - coord2[0])**2 + 
            (coord1[1] - coord2[1])**2 + 
            (coord1[2] - coord2[2])**2)**0.5


def estimate_pocket_volume(ligand_atoms: List[Dict[str, Any]]) -> float:
    """
    Estimate pocket volume based on ligand atoms.
    
    This is a simple approximation based on the bounding box of ligand atoms.
    
    Args:
        ligand_atoms: List of ligand atom dictionaries
        
    Returns:
        Estimated volume in cubic Angstroms
    """
    if not ligand_atoms:
        return 0.0
    
    x_coords = [atom['x'] for atom in ligand_atoms]
    y_coords = [atom['y'] for atom in ligand_atoms]
    z_coords = [atom['z'] for atom in ligand_atoms]
    
    # Calculate bounding box
    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    z_range = max(z_coords) - min(z_coords)
    
    # Add padding for pocket volume (ligand + binding site)
    padding = 5.0  # Angstroms
    volume = (x_range + padding) * (y_range + padding) * (z_range + padding)
    
    return volume


def get_first_drug_like_ligand(structure_info: Dict[str, Any]) -> Optional[str]:
    """
    Get the ID of the first drug-like ligand from structure info.
    
    Args:
        structure_info: Structure information dictionary from the pipeline
        
    Returns:
        Ligand ID (3-letter code) or None if no drug-like ligands
    """
    drug_like_ligands = structure_info.get("drug_like_ligands", [])
    
    if drug_like_ligands:
        # Return the ID of the first drug-like ligand
        first_ligand = drug_like_ligands[0]
        return first_ligand.get("id", first_ligand.get("name", "")).strip()
    
    return None
