"""
Tests for the PDB Processor
"""

import os
import unittest
import numpy as np
from unittest.mock import MagicMock, patch

from Bio.PDB import PDBParser, Structure

from agents.structure_predictor.pdb_processor import PDBProcessor


class TestPDBProcessor(unittest.TestCase):
    """Test suite for the PDB Processor."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.processor = PDBProcessor()
        
        # Create a mock structure
        self.mock_structure = MagicMock(spec=Structure)
        self.mock_chain = MagicMock()
        self.mock_residue = MagicMock()
        self.mock_atom = MagicMock()
        
        # Configure mocks
        self.mock_atom.get_coord.return_value = np.array([1.0, 2.0, 3.0])
        self.mock_atom.get_bfactor.return_value = 90.0
        self.mock_atom.get_parent.return_value = self.mock_residue
        
        self.mock_residue.get_id.return_value = (" ", 1, " ")
        self.mock_residue.get_resname.return_value = "ALA"
        self.mock_residue.get_atoms.return_value = [self.mock_atom]
        self.mock_residue.__iter__.return_value = iter([self.mock_atom])
        
        self.mock_chain.get_id.return_value = "A"
        self.mock_chain.get_residues.return_value = [self.mock_residue]
        self.mock_chain.__iter__.return_value = iter([self.mock_residue])
        
        self.mock_model = MagicMock()
        self.mock_model.get_chains.return_value = [self.mock_chain]
        self.mock_model.__iter__.return_value = iter([self.mock_chain])
        
        self.mock_structure.get_chains.return_value = [self.mock_chain]
        self.mock_structure.__iter__.return_value = iter([self.mock_model])
        
        # Patch the PDBParser
        self.patcher = patch("Bio.PDB.PDBParser")
        self.mock_parser_class = self.patcher.start()
        self.mock_parser = MagicMock()
        self.mock_parser_class.return_value = self.mock_parser
        self.mock_parser.get_structure.return_value = self.mock_structure
        
    def tearDown(self):
        """Tear down test fixtures."""
        self.patcher.stop()
    
    def test_parse_pdb(self):
        """Test parsing a PDB file."""
        # Create a temporary PDB file
        pdb_file = "tests/data/test.pdb"
        os.makedirs(os.path.dirname(pdb_file), exist_ok=True)
        with open(pdb_file, "w") as f:
            f.write("ATOM      1  N   ALA A   1      11.804  18.255  17.872  1.00 80.00           N\n")
        
        # Run parsing
        result = self.processor.parse_pdb(pdb_file)
        
        # Verify results
        self.assertIsNotNone(result)
        self.assertIn("chains", result)
        self.assertIn("residue_count", result)
        self.assertIn("atom_count", result)
        self.assertIn("secondary_structure", result)
        
        # Verify parser was called
        self.mock_parser.get_structure.assert_called_once_with("structure", pdb_file)
    
    def test_prepare_for_docking(self):
        """Test preparing a structure for docking."""
        # Create a temporary PDB file
        pdb_file = "tests/data/test.pdb"
        output_file = "tests/data/test_prepared.pdb"
        os.makedirs(os.path.dirname(pdb_file), exist_ok=True)
        with open(pdb_file, "w") as f:
            f.write("ATOM      1  N   ALA A   1      11.804  18.255  17.872  1.00 80.00           N\n")
        
        # Run preparation
        center = [0.0, 0.0, 0.0]
        radius = 10.0
        
        # Mock the PDBIO class
        with patch("Bio.PDB.PDBIO") as mock_pdbio_class:
            mock_pdbio = MagicMock()
            mock_pdbio_class.return_value = mock_pdbio
            
            result = self.processor.prepare_for_docking(pdb_file, output_file, center, radius)
        
        # Verify results
        self.assertIsNotNone(result)
        self.assertTrue(isinstance(result, str))
        
        # Verify parser and PDBIO were called
        self.mock_parser.get_structure.assert_called_once_with("structure", pdb_file)
    
    def test_extract_binding_sites(self):
        """Test extracting binding sites."""
        # Create a temporary PDB file
        pdb_file = "tests/data/test.pdb"
        os.makedirs(os.path.dirname(pdb_file), exist_ok=True)
        with open(pdb_file, "w") as f:
            f.write("ATOM      1  N   ALA A   1      11.804  18.255  17.872  1.00 80.00           N\n")
        
        # Define known sites
        known_sites = [
            {
                "name": "Site 1",
                "type": "ligand",
                "residues": [1, 2, 3]
            }
        ]
        
        # Run extraction
        result = self.processor.extract_binding_sites(pdb_file, known_sites)
        
        # Verify results
        self.assertIsNotNone(result)
        self.assertTrue(isinstance(result, list))
        
        # Verify parser was called
        self.mock_parser.get_structure.assert_called_once_with("structure", pdb_file)


if __name__ == "__main__":
    unittest.main()
