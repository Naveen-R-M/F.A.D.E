"""
PDB Processor for F.A.D.E

This module handles parsing, validation, and preparation of PDB files
for downstream processing, including docking.
"""

import os
from typing import Any, Dict, List, Optional, Tuple, Union
import re
import json
from collections import defaultdict

from Bio.PDB import PDBParser, PDBIO, Select, Structure, Chain, Atom, NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np
import warnings

# Suppress Bio.PDB warnings about missing atoms
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

class PDBProcessor:
    """
    Class for parsing, analyzing, and preparing PDB files.
    """
    
    def __init__(self) -> None:
        """Initialize the PDB processor."""
        self.parser = PDBParser(QUIET=True)
        self.io = PDBIO()
        
    def parse_pdb(self, pdb_file: str) -> Dict[str, Any]:
        """
        Parse a PDB file and extract key information.
        
        Args:
            pdb_file: Path to the PDB file.
            
        Returns:
            Dictionary containing parsed information.
        """
        # Parse the PDB file
        structure = self.parser.get_structure("structure", pdb_file)
        
        # Extract basic information
        chains = {}
        residue_count = 0
        atom_count = 0
        
        # Dictionary to store secondary structure assignments
        secondary_structure = {
            "helix": [],
            "sheet": [],
            "loop": []
        }
        
        # Parse chain information
        for chain in structure.get_chains():
            chain_id = chain.get_id()
            
            # Count residues and atoms
            chain_residues = list(chain.get_residues())
            chain_residue_count = len(chain_residues)
            residue_count += chain_residue_count
            
            # Extract residue information
            residues = []
            
            for residue in chain_residues:
                if residue.get_id()[0] != " ":  # Skip hetero residues
                    continue
                    
                residue_id = residue.get_id()[1]
                residue_name = residue.get_resname()
                
                # Count atoms
                residue_atoms = list(residue.get_atoms())
                atom_count += len(residue_atoms)
                
                # Get residue center
                coords = np.array([atom.get_coord() for atom in residue_atoms])
                residue_center = coords.mean(axis=0) if len(coords) > 0 else np.array([0, 0, 0])
                
                # Assign provisional secondary structure based on phi/psi angles
                # This is a simplistic approach - for production use, use DSSP or similar
                sec_struct = self._estimate_secondary_structure(residue, chain_residues)
                
                # Add to appropriate secondary structure list
                if sec_struct == "helix":
                    secondary_structure["helix"].append(residue_id)
                elif sec_struct == "sheet":
                    secondary_structure["sheet"].append(residue_id)
                else:
                    secondary_structure["loop"].append(residue_id)
                
                residues.append({
                    "id": residue_id,
                    "name": residue_name,
                    "center": residue_center.tolist(),
                    "secondary_structure": sec_struct
                })
            
            # Add chain information
            chains[chain_id] = {
                "residue_count": chain_residue_count,
                "residues": residues
            }
        
        # Extract model confidence scores from B-factors (for AlphaFold models)
        confidence_scores = self._extract_confidence_scores(structure)
        
        # Create result dictionary
        result = {
            "chains": chains,
            "residue_count": residue_count,
            "atom_count": atom_count,
            "secondary_structure": secondary_structure,
            "confidence_scores": confidence_scores
        }
        
        return result
    
    def prepare_for_docking(
        self,
        pdb_file: str,
        output_file: str,
        center: Optional[List[float]] = None,
        radius: float = 20.0
    ) -> str:
        """
        Prepare a protein structure for docking.
        
        Args:
            pdb_file: Path to the input PDB file.
            output_file: Path to save the prepared PDB file.
            center: Optional center coordinates [x, y, z] for the binding site.
                   If not provided, the geometric center of the protein is used.
            radius: Radius around the center to include in the binding site.
            
        Returns:
            Path to the prepared PDB file.
        """
        # Parse the PDB file
        structure = self.parser.get_structure("structure", pdb_file)
        
        # Process structure
        self._remove_hetero_atoms(structure)
        self._add_hydrogens(structure)
        self._assign_charges(structure)
        
        # Write full structure
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        self.io.set_structure(structure)
        self.io.save(output_file)
        
        # If center is provided, create a binding site file
        if center:
            binding_site_file = output_file.replace(".pdb", "_site.pdb")
            
            # Create binding site selector
            binding_site_selector = BindingSiteSelect(structure, center, radius)
            
            # Save binding site
            self.io.set_structure(structure)
            self.io.save(binding_site_file, binding_site_selector)
            
            # Create binding site config
            config_file = output_file.replace(".pdb", "_config.json")
            with open(config_file, "w") as f:
                json.dump({
                    "receptor": output_file,
                    "binding_site": binding_site_file,
                    "center": center,
                    "radius": radius
                }, f, indent=2)
                
            return config_file
        
        return output_file
    
    def extract_binding_sites(
        self,
        pdb_file: str,
        known_sites: Optional[List[Dict[str, Any]]] = None
    ) -> List[Dict[str, Any]]:
        """
        Extract binding sites from a PDB file.
        
        Args:
            pdb_file: Path to the PDB file.
            known_sites: Optional list of known binding sites.
            
        Returns:
            List of dictionaries containing binding site information.
        """
        # Parse the PDB file
        structure = self.parser.get_structure("structure", pdb_file)
        
        binding_sites = []
        
        # Process known binding sites
        if known_sites:
            for site in known_sites:
                site_name = site.get("name", "Unknown")
                site_residues = site.get("residues", [])
                
                # Find residues
                site_coords = []
                residue_ids = []
                
                for chain in structure.get_chains():
                    for residue in chain.get_residues():
                        res_id = residue.get_id()[1]
                        if res_id in site_residues:
                            # Add residue atoms to coordinates
                            for atom in residue.get_atoms():
                                site_coords.append(atom.get_coord())
                            residue_ids.append(res_id)
                
                # Calculate center and radius
                if site_coords:
                    site_coords = np.array(site_coords)
                    center = site_coords.mean(axis=0)
                    distances = np.sqrt(np.sum((site_coords - center)**2, axis=1))
                    radius = np.max(distances)
                    
                    binding_sites.append({
                        "name": site_name,
                        "type": site.get("type", "Known"),
                        "center": {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])},
                        "radius": float(radius),
                        "residue_ids": residue_ids,
                        "score": 1.0,  # Known sites get maximum score
                        "description": f"Known binding site: {site_name}"
                    })
        
        # Auto-detect potential binding sites
        auto_sites = self._detect_pocket_sites(structure)
        
        # Merge with known sites
        for site in auto_sites:
            # Check if this site overlaps with any known site
            overlaps = False
            for known_site in binding_sites:
                if self._sites_overlap(site, known_site):
                    overlaps = True
                    break
                    
            if not overlaps:
                binding_sites.append(site)
        
        return binding_sites
    
    def _remove_hetero_atoms(self, structure: Structure) -> None:
        """
        Remove hetero atoms from a structure.
        
        Args:
            structure: BioPython Structure object.
        """
        # Get all residues
        all_residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    all_residues.append((chain, residue))
        
        # Remove hetero residues
        for chain, residue in all_residues:
            if residue.get_id()[0] != " ":  # Hetero residue indicator
                chain.detach_child(residue.get_id())
    
    def _add_hydrogens(self, structure: Structure) -> None:
        """
        Add hydrogens to a structure.
        
        Args:
            structure: BioPython Structure object.
            
        Note:
            This is a placeholder. In a real implementation, you would use
            a more sophisticated method to add hydrogens, e.g., using PDBFixer,
            OpenBabel, or Schrodinger's PrepWizard.
        """
        # Placeholder for hydrogen addition
        # In a real implementation, you would use external tools
        pass
    
    def _assign_charges(self, structure: Structure) -> None:
        """
        Assign charges to atoms in a structure.
        
        Args:
            structure: BioPython Structure object.
            
        Note:
            This is a placeholder. In a real implementation, you would use
            a more sophisticated method to assign charges, e.g., using
            Schrodinger's PrepWizard or AMBER's Antechamber.
        """
        # Placeholder for charge assignment
        # In a real implementation, you would use external tools
        pass
    
    def _estimate_secondary_structure(
        self, 
        residue, 
        chain_residues
    ) -> str:
        """
        Estimate secondary structure based on crude analysis.
        
        Args:
            residue: BioPython Residue object.
            chain_residues: List of all residues in the chain.
            
        Returns:
            Secondary structure type: "helix", "sheet", or "loop".
            
        Note:
            This is a very simplistic approach. In a real implementation,
            you would use DSSP or a similar algorithm.
        """
        # Default to loop
        return "loop"
    
    def _extract_confidence_scores(self, structure: Structure) -> Dict[str, float]:
        """
        Extract confidence scores from B-factors (for AlphaFold models).
        
        Args:
            structure: BioPython Structure object.
            
        Returns:
            Dictionary containing confidence scores.
        """
        # Get all B-factors
        b_factors = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        for atom in residue:
                            b_factors.append(atom.get_bfactor())
        
        # Calculate statistics
        if b_factors:
            mean_b = np.mean(b_factors)
            max_b = np.max(b_factors)
            min_b = np.min(b_factors)
            
            # For AlphaFold, B-factors represent pLDDT scores (0-100)
            # Higher is better
            confidence = mean_b / 100.0
            
            return {
                "mean_b_factor": float(mean_b),
                "max_b_factor": float(max_b),
                "min_b_factor": float(min_b),
                "confidence": float(confidence)
            }
        
        return {
            "mean_b_factor": 0.0,
            "max_b_factor": 0.0,
            "min_b_factor": 0.0,
            "confidence": 0.0
        }
    
    def _detect_pocket_sites(self, structure: Structure) -> List[Dict[str, Any]]:
        """
        Detect potential binding pockets in the structure.
        
        Args:
            structure: BioPython Structure object.
            
        Returns:
            List of dictionaries containing binding site information.
            
        Note:
            This is a simplistic approach to pocket detection. In a real
            implementation, you would use a more sophisticated algorithm
            like fpocket, SiteMap, or AutoSite.
        """
        # Placeholder for pocket detection
        # In a real implementation, you would use more sophisticated algorithms
        
        # Get all residues
        all_residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        all_residues.append(residue)
        
        # Calculate protein center
        all_coords = []
        for residue in all_residues:
            for atom in residue:
                all_coords.append(atom.get_coord())
                
        if not all_coords:
            return []
            
        all_coords = np.array(all_coords)
        protein_center = all_coords.mean(axis=0)
        
        # Example: return the center of the protein as a potential binding site
        return [{
            "name": "Protein Center",
            "type": "auto",
            "center": {"x": float(protein_center[0]), "y": float(protein_center[1]), "z": float(protein_center[2])},
            "radius": 10.0,  # Default radius
            "residue_ids": [],  # Will be populated by binding site detector
            "score": 0.5,  # Medium confidence
            "description": "Automatically detected binding site at protein center"
        }]
    
    def _sites_overlap(
        self, 
        site1: Dict[str, Any], 
        site2: Dict[str, Any]
    ) -> bool:
        """
        Check if two binding sites overlap.
        
        Args:
            site1: First binding site.
            site2: Second binding site.
            
        Returns:
            True if sites overlap, False otherwise.
        """
        # Extract centers
        center1 = np.array([
            site1["center"]["x"],
            site1["center"]["y"],
            site1["center"]["z"]
        ])
        
        center2 = np.array([
            site2["center"]["x"],
            site2["center"]["y"],
            site2["center"]["z"]
        ])
        
        # Calculate distance
        distance = np.sqrt(np.sum((center1 - center2)**2))
        
        # Check if distance is less than sum of radii
        return distance < (site1["radius"] + site2["radius"])


class BindingSiteSelect(Select):
    """
    Selection class for extracting binding site residues.
    """
    
    def __init__(
        self, 
        structure: Structure,
        center: List[float],
        radius: float
    ) -> None:
        """
        Initialize the binding site selector.
        
        Args:
            structure: BioPython Structure object.
            center: Center coordinates [x, y, z] of the binding site.
            radius: Radius around the center to include.
        """
        super().__init__()
        self.center = np.array(center)
        self.radius = radius
        
        # Build neighbor search
        all_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        for atom in residue:
                            all_atoms.append(atom)
        
        self.neighbor_search = NeighborSearch(all_atoms)
    
    def accept_residue(self, residue: Residue) -> bool:
        """
        Determine whether to include a residue in the output.
        
        Args:
            residue: BioPython Residue object.
            
        Returns:
            True if the residue should be included, False otherwise.
        """
        # Skip hetero residues
        if residue.get_id()[0] != " ":
            return False
            
        # Check if any atom is within radius of center
        for atom in residue:
            coord = atom.get_coord()
            distance = np.sqrt(np.sum((coord - self.center)**2))
            if distance <= self.radius:
                return True
                
        return False
