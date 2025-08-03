"""
Molecule Builder for F.A.D.E Molecule Generator Agent

This component handles the construction of molecules from SMILES strings
and the generation of their 3D structures using RDKit.
"""

import json
import io
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import SDWriter

from utils.logging import get_logger
from utils.gemini_client import GeminiClient


class MoleculeBuilder:
    """
    Builds and optimizes 3D molecular structures from SMILES strings.
    """

    def __init__(self, llm_client: GeminiClient) -> None:
        """
        Initialize the MoleculeBuilder.

        Args:
            llm_client: An initialized GeminiClient instance for LLM interactions.
        """
        self.llm_client = llm_client
        self.logger = get_logger("fade.molecule_generator.molecule_builder")
        self.logger.info("MoleculeBuilder initialized.")

    def generate_smiles_from_binding_site(self, binding_site_info: Dict[str, Any]) -> List[str]:
        """
        Uses LLM to suggest SMILES strings for fragments or scaffolds based on binding site characteristics.

        Args:
            binding_site_info: Dictionary containing information about the binding site.

        Returns:
            A list of suggested SMILES strings.
        """
        prompt = f"""
        Given the following protein binding site characteristics, suggest 3-5 small molecule SMILES strings
        that would be chemically complementary. Focus on fragments or simple scaffolds.

        Binding Site Info: {json.dumps(binding_site_info, indent=2)}

        Provide only a JSON array of SMILES strings, like this: ["CCO", "C1=CC=CC=C1"]
        """
        self.logger.info("Querying LLM for SMILES suggestions based on binding site.")
        try:
            response_text = self.llm_client.generate_text(prompt, temperature=0.7, max_tokens=500)
            smiles_list = json.loads(response_text)
            if not isinstance(smiles_list, list) or not all(isinstance(s, str) for s in smiles_list):
                raise ValueError("LLM response was not a list of strings.")
            self.logger.info(f"LLM suggested {len(smiles_list)} SMILES strings.")
            return smiles_list
        except Exception as e:
            self.logger.error(f"Failed to get SMILES from LLM: {e}")
            return []

    def build_molecule_from_smiles(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Builds an RDKit molecule object from a SMILES string.

        Args:
            smiles: SMILES string of the molecule.

        Returns:
            An RDKit Mol object or None if invalid SMILES.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.logger.warning(f"Invalid SMILES string: {smiles}")
            return mol
        except Exception as e:
            self.logger.error(f"Error building molecule from SMILES {smiles}: {e}")
            return None

    def generate_3d_conformation(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Generates a 3D conformation for an RDKit molecule and optimizes it.

        Args:
            mol: RDKit Mol object.

        Returns:
            RDKit Mol object with 3D conformation or None if generation fails.
        """
        if mol is None:
            return None

        mol = Chem.AddHs(mol)  # Add hydrogens for 3D generation
        try:
            # Generate initial 3D coordinates
            AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Use ETKDG for better conformer generation
            # Optimize the conformation using UFF
            AllChem.UFFOptimizeMolecule(mol)
            self.logger.info("Generated and optimized 3D conformation.")
            return mol
        except Exception as e:
            self.logger.error(f"Error generating 3D conformation: {e}")
            return None

    def mol_to_sdf_string(self, mol: Chem.Mol, mol_name: str = "molecule") -> Optional[str]:
        """
        Converts an RDKit molecule to an SDF string.

        Args:
            mol: RDKit Mol object.
            mol_name: Name to assign to the molecule in the SDF file.

        Returns:
            SDF formatted string or None if conversion fails.
        """
        if mol is None:
            return None
        try:
            # Create an in-memory file-like object
            sdf_stream = io.StringIO()
            writer = SDWriter(sdf_stream)
            writer.SetForceV3000(True) # Force V3000 format
            
            # Write molecule properties to SDF
            for prop_name in mol.GetPropNames():
                mol.SetProp(prop_name, mol.GetProp(prop_name))
            
            mol.SetProp("_Name", mol_name)
            writer.write(mol)
            writer.flush()
            writer.close()
            
            sdf_string = sdf_stream.getvalue()
            return sdf_string
        except Exception as e:
            self.logger.error(f"Error converting molecule to SDF: {e}")
            return None
