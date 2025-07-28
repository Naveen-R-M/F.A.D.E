"""
Tests for the Structure Predictor Agent
"""

import os
import unittest
import json
from unittest.mock import MagicMock, patch

from agents.structure_predictor import StructurePredictor


class TestStructurePredictor(unittest.TestCase):
    """Test suite for the Structure Predictor agent."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a mock config
        self.mock_config = {
            "data_dir": "tests/data"
        }
        
        # Create a mock LLM client
        self.mock_llm_client = MagicMock()
        self.mock_llm_client.generate_text.return_value = json.dumps({
            "quality_assessment": "Good quality structure",
            "recommendations": ["No improvements needed"],
            "suitability": "Suitable for drug discovery",
            "concerns": [],
            "confidence": 0.9
        })
        
        # Create a mock SLURM client
        self.mock_slurm_client = MagicMock()
        self.mock_slurm_client.submit_job.return_value = "12345"
        self.mock_slurm_client.monitor_job.return_value = "COMPLETED"
        
        # Create a mock AlphaFold client
        self.mock_alphafold_client = MagicMock()
        self.mock_alphafold_client.get_best_model_path.return_value = "tests/data/outputs/structures/test_protein/model.pdb"
        
        # Create test directories if they don't exist
        os.makedirs("tests/data/outputs/structures/test_protein", exist_ok=True)
        os.makedirs("tests/data/outputs/docking/test_protein", exist_ok=True)
        
        # Patch the PDB processor, structure validator, and binding site detector
        self.mock_pdb_processor = MagicMock()
        self.mock_pdb_processor.parse_pdb.return_value = {
            "chains": {"A": {"residue_count": 100, "residues": []}},
            "residue_count": 100,
            "atom_count": 1000,
            "secondary_structure": {"helix": [], "sheet": [], "loop": []}
        }
        self.mock_pdb_processor.prepare_for_docking.return_value = "tests/data/outputs/docking/test_protein/test_protein_prepared.pdb"
        
        self.mock_structure_validator = MagicMock()
        self.mock_structure_validator.validate.return_value = {
            "is_valid": True,
            "overall_score": 0.9,
            "validation_scores": {
                "plddt": 0.9,
                "ptm": 0.8,
                "missing_residues_score": 0.95,
                "geometry_score": 0.85
            },
            "missing_residues": [],
            "geometry_issues": []
        }
        
        self.mock_binding_site_detector = MagicMock()
        self.mock_binding_site_detector.detect.return_value = [
            {
                "name": "Site 1",
                "type": "geometric",
                "center": {"x": 0.0, "y": 0.0, "z": 0.0},
                "radius": 10.0,
                "score": 0.8,
                "description": "Detected binding site",
                "residue_ids": [10, 11, 12, 13, 14]
            }
        ]
        
        # Create the agent with mocked components
        with patch("agents.structure_predictor.pdb_processor.PDBProcessor", return_value=self.mock_pdb_processor), \
             patch("agents.structure_predictor.structure_validator.StructureValidator", return_value=self.mock_structure_validator), \
             patch("agents.structure_predictor.binding_site_detector.BindingSiteDetector", return_value=self.mock_binding_site_detector):
            
            self.agent = StructurePredictor(
                config=self.mock_config,
                gemini_api_key="fake_key",
                gemini_model="fake_model"
            )
            
            # Replace clients with mocks
            self.agent.gemini_client = self.mock_llm_client
            self.agent.slurm_client = self.mock_slurm_client
            self.agent.alphafold_client = self.mock_alphafold_client
    
    def test_initialization(self):
        """Test agent initialization."""
        self.assertEqual(self.agent.name, "structure_predictor")
        self.assertEqual(self.agent.config, self.mock_config)
        self.assertIsNotNone(self.agent.gemini_client)
        self.assertIsNotNone(self.agent.slurm_client)
        self.assertIsNotNone(self.agent.alphafold_client)
        self.assertIsNotNone(self.agent.pdb_processor)
        self.assertIsNotNone(self.agent.structure_validator)
        self.assertIsNotNone(self.agent.binding_site_detector)
    
    def test_predict_structure(self):
        """Test structure prediction."""
        # Create test input
        target_name = "test_protein"
        sequence_info = {
            "sequence": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        }
        job_script_path = "tests/data/inputs/job_scripts/test_protein_alphafold.sh"
        
        # Create job script file
        os.makedirs(os.path.dirname(job_script_path), exist_ok=True)
        with open(job_script_path, "w") as f:
            f.write("#!/bin/bash\n# Test job script\n")
        
        # Run prediction
        structure_info = self.agent.predict_structure(target_name, sequence_info, job_script_path)
        
        # Verify results
        self.assertIsNotNone(structure_info)
        self.assertEqual(structure_info["target_name"], target_name)
        self.assertEqual(structure_info["pdb_file"], "tests/data/outputs/structures/test_protein/model.pdb")
        self.assertEqual(structure_info["residue_count"], 100)
        self.assertEqual(structure_info["atom_count"], 1000)
        self.assertAlmostEqual(structure_info["confidence_scores"]["overall"], 0.9)
        
        # Verify method calls
        self.mock_slurm_client.submit_job.assert_called_once_with(job_script_path)
        self.mock_slurm_client.monitor_job.assert_called_once_with("12345")
        self.mock_alphafold_client.get_best_model_path.assert_called_once()
        self.mock_pdb_processor.parse_pdb.assert_called_once()
        self.mock_structure_validator.validate.assert_called_once()
    
    def test_identify_binding_sites(self):
        """Test binding site identification."""
        # Create test input
        target_name = "test_protein"
        structure_info = {
            "target_name": target_name,
            "pdb_file": "tests/data/outputs/structures/test_protein/model.pdb",
            "residue_count": 100,
            "atom_count": 1000
        }
        sequence_info = {
            "sequence": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
            "binding_sites": [
                {
                    "name": "Known Site",
                    "type": "ligand",
                    "residues": [10, 11, 12, 13]
                }
            ]
        }
        
        # Run identification
        binding_sites = self.agent.identify_binding_sites(target_name, structure_info, sequence_info)
        
        # Verify results
        self.assertIsNotNone(binding_sites)
        self.assertEqual(len(binding_sites), 1)
        self.assertEqual(binding_sites[0]["name"], "Site 1")
        self.assertEqual(binding_sites[0]["type"], "geometric")
        self.assertEqual(binding_sites[0]["residue_ids"], [10, 11, 12, 13, 14])
        
        # Verify method calls
        self.mock_binding_site_detector.detect.assert_called_once()
    
    def test_prepare_structure_for_docking(self):
        """Test structure preparation for docking."""
        # Create test input
        target_name = "test_protein"
        structure_info = {
            "target_name": target_name,
            "pdb_file": "tests/data/outputs/structures/test_protein/model.pdb",
            "residue_count": 100,
            "atom_count": 1000
        }
        binding_sites = [
            {
                "name": "Site 1",
                "type": "geometric",
                "center": {"x": 0.0, "y": 0.0, "z": 0.0},
                "radius": 10.0,
                "score": 0.8,
                "description": "Detected binding site",
                "residue_ids": [10, 11, 12, 13, 14]
            }
        ]
        
        # Run preparation
        prepared_structure_path = self.agent.prepare_structure_for_docking(target_name, structure_info, binding_sites)
        
        # Verify results
        self.assertIsNotNone(prepared_structure_path)
        self.assertTrue(os.path.dirname(prepared_structure_path).endswith(f"docking/{target_name}"))
        
        # Verify method calls
        self.mock_pdb_processor.prepare_for_docking.assert_called_once()
    
    def test_process(self):
        """Test the complete process method."""
        # Create test input
        input_data = {
            "sequences": {
                "test_protein": {
                    "sequence": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
                    "binding_sites": [
                        {
                            "name": "Known Site",
                            "type": "ligand",
                            "residues": [10, 11, 12, 13]
                        }
                    ]
                }
            },
            "job_configs": {
                "test_protein_alphafold_job": "tests/data/inputs/job_scripts/test_protein_alphafold.sh"
            }
        }
        
        # Create job script file
        os.makedirs("tests/data/inputs/job_scripts", exist_ok=True)
        with open("tests/data/inputs/job_scripts/test_protein_alphafold.sh", "w") as f:
            f.write("#!/bin/bash\n# Test job script\n")
        
        # Run process
        result = self.agent.process(input_data)
        
        # Verify results
        self.assertIsNotNone(result)
        self.assertIn("structures", result)
        self.assertIn("binding_sites", result)
        self.assertIn("prepared_structures", result)
        self.assertIn("test_protein", result["structures"])
        self.assertIn("test_protein", result["binding_sites"])
        self.assertIn("test_protein", result["prepared_structures"])
        
        # Verify method calls
        self.mock_slurm_client.submit_job.assert_called_once()
        self.mock_slurm_client.monitor_job.assert_called_once()
        self.mock_alphafold_client.get_best_model_path.assert_called_once()
        self.mock_pdb_processor.parse_pdb.assert_called_once()
        self.mock_structure_validator.validate.assert_called_once()
        self.mock_binding_site_detector.detect.assert_called_once()
        self.mock_pdb_processor.prepare_for_docking.assert_called_once()


if __name__ == "__main__":
    unittest.main()
