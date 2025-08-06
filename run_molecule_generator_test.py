import os
import sys
import json
from unittest.mock import MagicMock, patch

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from agents.molecule_generator.molecule_generator import MoleculeGenerator
from utils.gemini_client import GeminiClient
from utils.logging import setup_logging, get_logger

# Setup basic logging for the script
setup_logging(log_level="INFO", logs_dir="/tmp/molgen_test_logs")
logger = get_logger("test_molecule_generator_script")

def run_molecule_generator_test():
    logger.info("Starting Molecule Generator test script...")

    # Mock GeminiClient for tests
    mock_gemini_client = MagicMock(spec=GeminiClient)
    # Configure mock to return a valid JSON list of SMILES strings
    mock_gemini_client.generate_text.return_value = json.dumps(["CCO", "CC(=O)O", "C1=CC=C(C=C1)C(=O)O"]) # Ethanol, Acetic Acid, Benzoic Acid

    # Define mock inputs
    target_name = "KRAS_test_molgen"
    base_output_dir = "/scratch/rajagopalmohanraj.n/F.A.D.E/data/outputs/docking"
    
    mock_prepared_structures = {
        target_name: os.path.join(base_output_dir, target_name, "KRAS_prepared.json")
    }

    mock_binding_sites = {
        target_name: [
            {
                "name": "GTP_pocket",
                "center": {"x": 10, "y": 20, "z": 30},
                "radius": 10.0,
                "residue_ids": [10, 12, 13, 16, 117, 146], # Example KRAS residues
                "description": "Known GTP binding pocket"
            }
        ]
    }

    mock_parsed_query_data = {
        "molecule_properties": {
            "blood_brain_barrier_permeability": True,
            "low_toxicity": True,
            "lipinski_rule_of_five": True
        }
    }

    input_data = {
        "prepared_structures": mock_prepared_structures,
        "binding_sites": mock_binding_sites,
        "parsed_query_data": mock_parsed_query_data
    }

    # Initialize MoleculeGenerator and inject the mocked GeminiClient
    with patch('agents.molecule_generator.molecule_generator.GeminiClient', return_value=mock_gemini_client):
        generator = MoleculeGenerator(
            gemini_api_key="dummy_key", 
            gemini_model="dummy_model", 
            config={"data_dir": "/scratch/rajagopalmohanraj.n/F.A.D.E/data"}
        )

        logger.info("Running MoleculeGenerator.process...")
        results = generator.process(input_data)

        logger.info("Molecule Generator test completed. Results:")
        print(json.dumps(results, indent=2))

    # Optional: Clean up mock files if desired after verification
    # os.remove(mock_prepared_structures[target_name])
    # os.remove(os.path.join(base_output_dir, target_name, "KRAS_GTP_pocket_prepared.pdb"))
    # os.rmdir(os.path.join(base_output_dir, target_name))

if __name__ == "__main__":
    run_molecule_generator_test()