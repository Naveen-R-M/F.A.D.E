"""
Binding Site Detector for F.A.D.E

This module identifies potential binding sites in protein structures
using geometric analysis, conservation data, and knowledge of known sites.
"""

import os
from typing import Any, Dict, List, Optional, Tuple, Union
import json
import numpy as np

from Bio.PDB import PDBParser, Structure, NeighborSearch
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

# Suppress Bio.PDB warnings about missing atoms
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


class BindingSiteDetector:
    """
    Class for detecting potential binding sites in protein structures.
    """
    
    def __init__(self, llm_client=None) -> None:
        """
        Initialize the binding site detector.
        
        Args:
            llm_client: Optional client for LLM integration.
        """
        self.parser = PDBParser(QUIET=True)
        self.llm_client = llm_client
    
    def detect(
        self,
        pdb_file: str,
        known_binding_sites: Optional[List[Dict[str, Any]]] = None,
        detection_methods: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Detect potential binding sites in a protein structure.
        
        Args:
            pdb_file: Path to the PDB file.
            known_binding_sites: Optional list of known binding sites.
            detection_methods: Optional list of detection methods to use.
                               Defaults to ["geometric", "conservation", "known"].
            
        Returns:
            List of dictionaries containing binding site information.
        """
        # Set default detection methods
        if detection_methods is None:
            detection_methods = ["geometric", "conservation", "known"]
            
        # Parse the PDB file
        structure = self.parser.get_structure("structure", pdb_file)
        
        # Initialize results
        binding_sites = []
        
        # Process known binding sites
        if "known" in detection_methods and known_binding_sites:
            known_sites = self._process_known_sites(structure, known_binding_sites)
            binding_sites.extend(known_sites)
        
        # Detect sites using geometric analysis
        if "geometric" in detection_methods:
            geometric_sites = self._detect_geometric_sites(structure)
            
            # Filter out geometric sites that overlap with known sites
            filtered_geometric_sites = []
            for site in geometric_sites:
                if not self._overlaps_with_existing_sites(site, binding_sites):
                    filtered_geometric_sites.append(site)
                    
            binding_sites.extend(filtered_geometric_sites)
        
        # Detect sites using conservation analysis
        if "conservation" in detection_methods:
            conservation_sites = self._detect_conservation_sites(structure)
            
            # Filter out conservation sites that overlap with existing sites
            filtered_conservation_sites = []
            for site in conservation_sites:
                if not self._overlaps_with_existing_sites(site, binding_sites):
                    filtered_conservation_sites.append(site)
                    
            binding_sites.extend(filtered_conservation_sites)
        
        # Use LLM to analyze results if available
        if self.llm_client and binding_sites:
            binding_sites = self._analyze_with_llm(structure, binding_sites)
        
        # Calculate residue IDs for each site
        binding_sites = self._calculate_residue_ids(structure, binding_sites)
        
        return binding_sites
    
    def _process_known_sites(
        self,
        structure: Structure,
        known_sites: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Process known binding sites.
        
        Args:
            structure: BioPython Structure object.
            known_sites: List of known binding sites.
            
        Returns:
            List of dictionaries containing binding site information.
        """
        processed_sites = []
        
        for site in known_sites:
            # Extract basic info
            site_name = site.get("name", "Unknown")
            site_type = site.get("type", "Known")
            residues = site.get("residues", [])
            
            # Find residues in structure
            site_residues = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[0] == " ":  # Skip hetero residues
                            res_id = residue.get_id()[1]
                            if res_id in residues:
                                site_residues.append(residue)
            
            # Calculate center and radius
            if site_residues:
                # Get all atom coordinates
                all_coords = []
                for residue in site_residues:
                    for atom in residue:
                        all_coords.append(atom.get_coord())
                
                if all_coords:
                    all_coords = np.array(all_coords)
                    center = all_coords.mean(axis=0)
                    distances = np.sqrt(np.sum((all_coords - center)**2, axis=1))
                    radius = np.max(distances)
                    
                    # Create site info
                    processed_sites.append({
                        "name": site_name,
                        "type": site_type,
                        "center": {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])},
                        "radius": float(radius),
                        "score": 1.0,  # Known sites get maximum score
                        "description": f"Known binding site: {site_name}"
                    })
        
        return processed_sites
    
    def _detect_geometric_sites(self, structure: Structure) -> List[Dict[str, Any]]:
        """
        Detect binding sites using geometric analysis.
        
        Args:
            structure: BioPython Structure object.
            
        Returns:
            List of dictionaries containing binding site information.
            
        Note:
            This is a simplistic approach. In a real implementation, you would use
            a more sophisticated algorithm like fpocket, SiteMap, or AutoSite.
        """
        # Get all atoms
        all_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        for atom in residue:
                            all_atoms.append(atom)
        
        if not all_atoms:
            return []
        
        # Calculate protein center
        all_coords = np.array([atom.get_coord() for atom in all_atoms])
        protein_center = all_coords.mean(axis=0)
        
        # Detect cavities
        # This is a placeholder implementation. In a real scenario, you would
        # use more sophisticated algorithms to identify cavities.
        
        # Example: Identify the geometric center as a potential binding site
        sites = [{
            "name": "Geometric Center",
            "type": "geometric",
            "center": {"x": float(protein_center[0]), "y": float(protein_center[1]), "z": float(protein_center[2])},
            "radius": 10.0,  # Default radius
            "score": 0.5,  # Medium confidence
            "description": "Detected binding site at protein geometric center"
        }]
        
        return sites
    
    def _detect_conservation_sites(self, structure: Structure) -> List[Dict[str, Any]]:
        """
        Detect binding sites using conservation analysis.
        
        Args:
            structure: BioPython Structure object.
            
        Returns:
            List of dictionaries containing binding site information.
            
        Note:
            This is a placeholder. In a real implementation, you would use
            conservation data from tools like ConSurf or evolutionary information.
        """
        # Placeholder - in a real implementation, you would use conservation data
        return []
    
    def _overlaps_with_existing_sites(
        self,
        new_site: Dict[str, Any],
        existing_sites: List[Dict[str, Any]]
    ) -> bool:
        """
        Check if a new site overlaps with any existing sites.
        
        Args:
            new_site: New binding site.
            existing_sites: List of existing binding sites.
            
        Returns:
            True if the new site overlaps with any existing site, False otherwise.
        """
        new_center = np.array([
            new_site["center"]["x"],
            new_site["center"]["y"],
            new_site["center"]["z"]
        ])
        new_radius = new_site["radius"]
        
        for site in existing_sites:
            site_center = np.array([
                site["center"]["x"],
                site["center"]["y"],
                site["center"]["z"]
            ])
            site_radius = site["radius"]
            
            # Calculate distance
            distance = np.sqrt(np.sum((new_center - site_center)**2))
            
            # Check if distance is less than sum of radii
            if distance < (new_radius + site_radius):
                return True
                
        return False
    
    def _analyze_with_llm(
        self,
        structure: Structure,
        binding_sites: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Use LLM to analyze and refine binding sites.
        
        Args:
            structure: BioPython Structure object.
            binding_sites: List of detected binding sites.
            
        Returns:
            List of refined binding sites.
        """
        if not self.llm_client:
            return binding_sites
            
        # Extract basic structure info
        chains = []
        residue_count = 0
        atom_count = 0
        
        for model in structure:
            for chain in model:
                chains.append(chain.get_id())
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        residue_count += 1
                        atom_count += len(list(residue.get_atoms()))
        
        # Prepare prompt
        prompt = f"""
        Analyze the following protein structure and detected binding sites:
        
        Protein Structure Info:
        - Chains: {", ".join(chains)}
        - Residue Count: {residue_count}
        - Atom Count: {atom_count}
        
        Detected Binding Sites:
        {json.dumps(binding_sites, indent=2)}
        
        Please analyze these binding sites and:
        1. Rank them by likely biological significance
        2. Provide a confidence score for each (0.0 to 1.0)
        3. Add a detailed description for each site
        4. Note any overlapping sites that should be merged
        
        Format your response as a JSON with the following structure:
        {{
            "ranked_sites": [
                {{
                    "name": "Site name",
                    "type": "Site type",
                    "center": {{"x": X, "y": Y, "z": Z}},
                    "radius": Radius in Angstroms,
                    "score": Confidence score (0.0 to 1.0),
                    "description": "Detailed description"
                }},
                ...
            ],
            "overlapping_sites": [
                [index1, index2, ...]  // Indices of sites that should be merged
            ]
        }}
        """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.2)
            
            # Extract JSON from response
            try:
                # Find JSON in the response
                json_start = response.find('{')
                json_end = response.rfind('}') + 1
                if json_start >= 0 and json_end > json_start:
                    json_str = response[json_start:json_end]
                    analysis = json.loads(json_str)
                    
                    # Update binding sites with ranked sites
                    if "ranked_sites" in analysis:
                        return analysis["ranked_sites"]
                    
                    # If no ranked sites, return original sites
                    return binding_sites
                else:
                    return binding_sites
            except json.JSONDecodeError:
                return binding_sites
        except Exception as e:
            # If anything goes wrong, return original sites
            return binding_sites
    
    def _calculate_residue_ids(
        self,
        structure: Structure,
        binding_sites: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Calculate residue IDs for each binding site.
        
        Args:
            structure: BioPython Structure object.
            binding_sites: List of binding sites.
            
        Returns:
            List of binding sites with residue IDs added.
        """
        # Create neighbor search for the structure
        all_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        for atom in residue:
                            all_atoms.append(atom)
        
        if not all_atoms:
            return binding_sites
            
        ns = NeighborSearch(all_atoms)
        
        # Calculate residue IDs for each site
        for site in binding_sites:
            center = np.array([
                site["center"]["x"],
                site["center"]["y"],
                site["center"]["z"]
            ])
            radius = site["radius"]
            
            # Find atoms within radius
            atoms = ns.search(center, radius)
            
            # Get unique residues
            residues = set()
            for atom in atoms:
                residue = atom.get_parent()
                if residue.get_id()[0] == " ":  # Skip hetero residues
                    residues.add(residue)
            
            # Extract residue IDs
            residue_ids = [residue.get_id()[1] for residue in residues]
            
            # Add residue IDs to site
            site["residue_ids"] = sorted(residue_ids)
        
        return binding_sites
