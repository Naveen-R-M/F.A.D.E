"""
Unit tests for TargetSelector agent
"""

import os
import json
import unittest
from unittest.mock import patch, MagicMock, mock_open

import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agents.target_selector import TargetSelector
from utils.gemini_client import GeminiClient
from utils.uniprot_client import UniProtClient
from utils.config_generator import ConfigGenerator


class TestTargetSelector(unittest.TestCase):
    """Test cases for TargetSelector agent"""
    
    def setUp(self):
        """Set up test fixtures"""
        # Create mock dependencies
        self.gemini_client_mock = MagicMock(spec=GeminiClient)
        self.uniprot_client_mock = MagicMock(spec=UniProtClient)
        self.config_generator_mock = MagicMock(spec=ConfigGenerator)
        
        # Create patcher for external dependencies
        self.gemini_patcher = patch('agents.target_selector.GeminiClient', return_value=self.gemini_client_mock)
        self.uniprot_patcher = patch('agents.target_selector.UniProtClient', return_value=self.uniprot_client_mock)
        self.config_patcher = patch('agents.target_selector.ConfigGenerator', return_value=self.config_generator_mock)
        
        # Start patchers
        self.gemini_mock = self.gemini_patcher.start()
        self.uniprot_mock = self.uniprot_patcher.start()
        self.config_mock = self.config_patcher.start()
        
        # Create test instance
        self.target_selector = TargetSelector(name="test_target_selector")
        
        # Create test data
        self.test_query = "Find molecules targeting KRAS G12D with good BBB permeability"
        self.test_parsed_data = {
            "protein_targets": [
                {
                    "name": "KRAS", 
                    "mutations": [
                        {"original_residue": "G", "position": 12, "mutated_residue": "D"}
                    ]
                }
            ],
            "molecule_properties": {
                "blood_brain_barrier_permeability": True
            }
        }
        self.test_protein = {
            "primaryAccession": "P01116",
            "sequence": {"value": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"},
            "organism": {"scientificName": "Homo sapiens"},
            "proteinDescription": {"recommendedName": {"fullName": {"value": "GTPase KRas"}}}
        }
        self.test_sequence_info = {
            "accession": "P01116",
            "original_sequence": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
            "sequence": "MTEYKLVVVGADGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
            "length": 188,
            "mutations": [{"original_residue": "G", "position": 12, "mutated_residue": "D"}],
            "organism": "Homo sapiens",
            "protein_name": "GTPase KRas",
            "gene_name": "KRAS"
        }
        
    def tearDown(self):
        """Tear down test fixtures"""
        # Stop patchers
        self.gemini_patcher.stop()
        self.uniprot_patcher.stop()
        self.config_patcher.stop()
    
    def test_init(self):
        """Test initialization of TargetSelector agent"""
        self.assertEqual(self.target_selector.name, "test_target_selector")
        self.assertIsNotNone(self.target_selector.gemini_client)
        self.assertIsNotNone(self.target_selector.uniprot_client)
        self.assertIsNotNone(self.target_selector.config_generator)
    
    def test_parse_query(self):
        """Test parsing of natural language query"""
        # Configure mock
        self.gemini_client_mock.extract_protein_info.return_value = self.test_parsed_data
        
        # Call method under test
        result = self.target_selector.parse_query(self.test_query)
        
        # Verify
        self.gemini_client_mock.extract_protein_info.assert_called_once_with(self.test_query)
        self.assertEqual(result, self.test_parsed_data)
    
    def test_fetch_protein_sequence(self):
        """Test fetching of protein sequence"""
        # Configure mocks
        self.uniprot_client_mock.get_protein_by_gene_name.return_value = self.test_protein
        self.uniprot_client_mock.apply_mutation.return_value = self.test_sequence_info["sequence"]
        
        target_info = {
            "name": "KRAS",
            "organism": "human",
            "mutations": [
                {"original_residue": "G", "position": 12, "mutated_residue": "D"}
            ]
        }
        
        # Call method under test
        result = self.target_selector.fetch_protein_sequence(target_info)
        
        # Verify
        self.uniprot_client_mock.get_protein_by_gene_name.assert_called_once_with("KRAS", "human")
        self.uniprot_client_mock.apply_mutation.assert_called_once()
        self.assertEqual(result["accession"], "P01116")
        self.assertEqual(result["gene_name"], "KRAS")
        self.assertIn("G12D", str(result["mutations"]))
    
    def test_fetch_protein_sequence_without_organism(self):
        """Test fetching of protein sequence without organism"""
        # Configure mocks
        self.uniprot_client_mock.get_protein_by_gene_name.return_value = None
        self.uniprot_client_mock.search_protein.return_value = [self.test_protein]
        self.uniprot_client_mock.apply_mutation.return_value = self.test_sequence_info["sequence"]
        
        target_info = {
            "name": "KRAS",
            "mutations": [
                {"original_residue": "G", "position": 12, "mutated_residue": "D"}
            ]
        }
        
        # Call method under test
        result = self.target_selector.fetch_protein_sequence(target_info)
        
        # Verify
        self.uniprot_client_mock.get_protein_by_gene_name.assert_called_once_with("KRAS", None)
        self.uniprot_client_mock.search_protein.assert_called_once_with("KRAS", limit=1)
        self.assertEqual(result["accession"], "P01116")
    
    @patch('os.makedirs')
    @patch('builtins.open', new_callable=mock_open)
    def test_generate_configs(self, mock_file, mock_makedirs):
        """Test generation of configuration files"""
        # Configure mocks
        self.target_selector._get_data_dir = MagicMock(return_value="/tmp/fade/data")
        
        parsed_data = self.test_parsed_data
        sequences = {"KRAS": self.test_sequence_info}
        
        # Call method under test
        result = self.target_selector.generate_configs(parsed_data, sequences)
        
        # Verify
        self.assertIn("KRAS_fasta", result)
        self.assertIn("KRAS_alphafold_config", result)
        self.assertIn("KRAS_alphafold_job", result)
        self.assertIn("parsed_query", result)
        mock_makedirs.assert_called()
        self.uniprot_client_mock.save_fasta.assert_called_once()
        self.config_generator_mock.create_alphafold_job.assert_called_once()
    
    @patch('os.path.dirname')
    def test_get_data_dir(self, mock_dirname):
        """Test getting the data directory"""
        # Configure mocks
        mock_dirname.side_effect = [
            "/tmp/fade/agents",
            "/tmp/fade"
        ]
        
        # Call method under test
        result = self.target_selector._get_data_dir()
        
        # Verify
        self.assertEqual(result, "/tmp/fade/data")
    
    @patch('os.makedirs')
    @patch('builtins.open', new_callable=mock_open)
    def test_process(self, mock_file, mock_makedirs):
        """Test the entire process flow"""
        # Configure mocks
        self.target_selector.parse_query = MagicMock(return_value=self.test_parsed_data)
        self.target_selector.fetch_protein_sequence = MagicMock(return_value=self.test_sequence_info)
        self.target_selector.generate_configs = MagicMock(return_value={
            "KRAS_fasta": "/tmp/fade/data/inputs/sequences/KRAS.fasta",
            "KRAS_alphafold_config": "/tmp/fade/data/inputs/configs/KRAS_alphafold.json",
            "KRAS_alphafold_job": "/tmp/fade/data/inputs/job_scripts/KRAS_alphafold.sh",
            "parsed_query": "/tmp/fade/data/inputs/configs/parsed_query.json"
        })
        
        # Call method under test
        result = self.target_selector.process(self.test_query)
        
        # Verify
        self.target_selector.parse_query.assert_called_once_with(self.test_query)
        self.target_selector.fetch_protein_sequence.assert_called_once()
        self.target_selector.generate_configs.assert_called_once()
        self.assertIn("parsed_data", result)
        self.assertIn("sequences", result)
        self.assertIn("config_files", result)
        self.assertEqual(result["parsed_data"], self.test_parsed_data)
        self.assertEqual(result["sequences"]["KRAS"], self.test_sequence_info)
        

if __name__ == "__main__":
    unittest.main()
