"""
Property Calculator for F.A.D.E Molecule Generator Agent

This component calculates various physicochemical properties for molecules
using RDKit.
"""

from typing import Any, Dict, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, Crippen

from utils.logging import get_logger


class PropertyCalculator:
    """
    Calculates a set of common physicochemical properties for RDKit molecules.
    """

    def __init__(self) -> None:
        """
        Initialize the PropertyCalculator.
        """
        self.logger = get_logger("fade.molecule_generator.property_calculator")
        self.logger.info("PropertyCalculator initialized.")

    def calculate_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculates a dictionary of physicochemical properties for a given RDKit molecule.

        Args:
            mol: RDKit Mol object.

        Returns:
            A dictionary where keys are property names and values are their calculated values.
        """
        if mol is None:
            self.logger.warning("Cannot calculate properties for a None molecule.")
            return {}

        properties = {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "h_bond_donors": Lipinski.NumHDonors(mol),
            "h_bond_acceptors": Lipinski.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "num_heteroatoms": Descriptors.NumHeteroatoms(mol),
            "formal_charge": Chem.GetFormalCharge(mol),
        }
        self.logger.debug(f"Calculated properties for molecule: {properties}")
        return properties
