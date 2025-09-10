"""
Molecule Filter for F.A.D.E Molecule Generator Agent

This component applies various filters to generated molecules based on
drug-likeness rules and user-defined criteria.
"""

from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import Lipinski

from utils.logging import get_logger


class MoleculeFilter:
    """
    Filters molecules based on drug-likeness rules (e.g., Lipinski's Rule of Five)
    and other specified criteria.
    """

    def __init__(self) -> None:
        """
        Initialize the MoleculeFilter.
        """
        self.logger = get_logger("fade.molecule_generator.molecule_filter")
        self.logger.info("MoleculeFilter initialized.")

    def apply_lipinski_filter(self, properties: Dict[str, Any]) -> bool:
        """
        Applies Lipinski's Rule of Five to a molecule's properties.

        Rules:
        - Molecular weight < 500 Da
        - LogP < 5
        - Hydrogen bond donors < 5
        - Hydrogen bond acceptors < 10

        Args:
            properties: Dictionary of molecule properties.

        Returns:
            True if the molecule passes Lipinski's Rule of Five, False otherwise.
        """
        mw = properties.get("molecular_weight", float('inf'))
        logp = properties.get("logp", float('inf'))
        hbd = properties.get("h_bond_donors", float('inf'))
        hba = properties.get("h_bond_acceptors", float('inf'))

        passes_lipinski = (
            mw < 500 and
            logp < 5 and
            hbd < 5 and
            hba < 10
        )

        if not passes_lipinski:
            self.logger.debug(f"Molecule failed Lipinski's Rule of Five: MW={mw}, LogP={logp}, HBD={hbd}, HBA={hba}")
        return passes_lipinski

    def apply_custom_filters(self, properties: Dict[str, Any], requirements: Dict[str, Any]) -> bool:
        """
        Applies custom filters based on user-defined requirements.

        Args:
            properties: Dictionary of molecule properties.
            requirements: Dictionary of user-defined molecular property requirements.

        Returns:
            True if the molecule passes all custom filters, False otherwise.
        """
        # Example custom filters (can be expanded based on query parsing)
        # BBB Permeability (simplified, typically more complex models are used)
        if requirements.get("blood_brain_barrier_permeability"):
            tpsa = properties.get("tpsa", float('inf'))
            # Simple heuristic: low TPSA for BBB permeability (e.g., TPSA < 90)
            if not (tpsa < 90):
                self.logger.debug(f"Molecule failed BBB permeability filter: TPSA={tpsa}")
                return False

        # Toxicity (placeholder - in a real scenario, this would involve ML models or databases)
        if requirements.get("toxicity_requirements"):
            # For now, assume no specific toxicity filter is applied here unless a simple rule is defined.
            # This will be handled more robustly by the Evaluator Agent.
            pass

        self.logger.debug("Molecule passed custom filters.")
        return True

    def filter_molecules(self, molecules_with_properties: List[Dict[str, Any]], requirements: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Filters a list of molecules based on Lipinski's Rule of Five and custom requirements.

        Args:
            molecules_with_properties: List of dictionaries, each containing a molecule
                                       (RDKit Mol) and its calculated properties.
            requirements: Dictionary of user-defined molecular property requirements.

        Returns:
            A list of molecules that passed all filters.
        """
        filtered_list = []
        for mol_data in molecules_with_properties:
            mol = mol_data["mol"]
            properties = mol_data["properties"]

            if self.apply_lipinski_filter(properties) and self.apply_custom_filters(properties, requirements):
                filtered_list.append(mol_data)
            else:
                self.logger.info(f"Molecule {mol_data.get('smiles', 'N/A')} failed filtering.")

        self.logger.info(f"Filtered {len(molecules_with_properties)} to {len(filtered_list)} molecules.")
        return filtered_list
