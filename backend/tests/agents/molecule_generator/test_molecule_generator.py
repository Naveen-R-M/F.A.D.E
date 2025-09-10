import pytest
import os
import json
from unittest.mock import MagicMock, patch

from rdkit import Chem
from rdkit.Chem import AllChem

from agents.molecule_generator.molecule_builder import MoleculeBuilder
from agents.molecule_generator.property_calculator import PropertyCalculator
from agents.molecule_generator.molecule_filter import MoleculeFilter
from agents.molecule_generator.molecule_generator import MoleculeGenerator
from utils.gemini_client import GeminiClient

# Mock GeminiClient for tests
@pytest.fixture
def mock_gemini_client():
    mock_client = MagicMock(spec=GeminiClient)
    # Default mock for generate_text, can be overridden in specific tests
    mock_client.generate_text.return_value = json.dumps(["CCO", "CC(=O)O"])
    return mock_client

# --- Tests for MoleculeBuilder ---
def test_molecule_builder_init(mock_gemini_client):
    builder = MoleculeBuilder(llm_client=mock_gemini_client)
    assert builder.llm_client == mock_gemini_client

def test_generate_smiles_from_binding_site(mock_gemini_client):
    mock_gemini_client.generate_text.return_value = json.dumps(["C1CCCCC1", "C=O"])
    builder = MoleculeBuilder(llm_client=mock_gemini_client)
    binding_site_info = {"description": "hydrophobic pocket"}
    smiles_list = builder.generate_smiles_from_binding_site(binding_site_info)
    assert isinstance(smiles_list, list)
    assert "C1CCCCC1" in smiles_list
    assert "C=O" in smiles_list
    mock_gemini_client.generate_text.assert_called_once()

def test_generate_smiles_from_binding_site_invalid_llm_response(mock_gemini_client):
    mock_gemini_client.generate_text.return_value = "not a json list"
    builder = MoleculeBuilder(llm_client=mock_gemini_client)
    smiles_list = builder.generate_smiles_from_binding_site({})
    assert smiles_list == []

def test_build_molecule_from_smiles_valid():
    builder = MoleculeBuilder(llm_client=MagicMock()) # LLM not used here
    mol = builder.build_molecule_from_smiles("CCO")
    assert mol is not None
    assert Chem.MolToSmiles(mol) == "CCO"

def test_build_molecule_from_smiles_invalid():
    builder = MoleculeBuilder(llm_client=MagicMock())
    mol = builder.build_molecule_from_smiles("invalid_smiles")
    assert mol is None

def test_generate_3d_conformation():
    builder = MoleculeBuilder(llm_client=MagicMock())
    mol = Chem.MolFromSmiles("CCO")
    mol_3d = builder.generate_3d_conformation(mol)
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() > 0

def test_generate_3d_conformation_none_input():
    builder = MoleculeBuilder(llm_client=MagicMock())
    mol_3d = builder.generate_3d_conformation(None)
    assert mol_3d is None

def test_mol_to_sdf_string():
    builder = MoleculeBuilder(llm_client=MagicMock())
    mol = Chem.MolFromSmiles("CCO")
    mol_3d = builder.generate_3d_conformation(mol)
    sdf_string = builder.mol_to_sdf_string(mol_3d, "ethanol")
    assert sdf_string is not None
    assert "ethanol" in sdf_string
    assert "V3000" in sdf_string # Indicates 3D SDF

def test_mol_to_sdf_string_none_input():
    builder = MoleculeBuilder(llm_client=MagicMock())
    sdf_string = builder.mol_to_sdf_string(None)
    assert sdf_string is None

# --- Tests for PropertyCalculator ---
def test_property_calculator_init():
    calculator = PropertyCalculator()
    assert isinstance(calculator, PropertyCalculator)

def test_calculate_properties():
    calculator = PropertyCalculator()
    mol = Chem.MolFromSmiles("CCO") # Ethanol
    properties = calculator.calculate_properties(mol)
    assert "molecular_weight" in properties
    assert "logp" in properties
    assert properties["molecular_weight"] == pytest.approx(46.07, 0.01)
    assert properties["logp"] == pytest.approx(-0.05, abs=0.05) # Increased tolerance for LogP

def test_calculate_properties_none_input():
    calculator = PropertyCalculator()
    properties = calculator.calculate_properties(None)
    assert properties == {}

# --- Tests for MoleculeFilter ---
def test_molecule_filter_init():
    mol_filter = MoleculeFilter()
    assert isinstance(mol_filter, MoleculeFilter)

def test_apply_lipinski_filter_pass():
    mol_filter = MoleculeFilter()
    properties = {"molecular_weight": 200, "logp": 2, "h_bond_donors": 1, "h_bond_acceptors": 3}
    assert mol_filter.apply_lipinski_filter(properties) is True

def test_apply_lipinski_filter_fail_mw():
    mol_filter = MoleculeFilter()
    properties = {"molecular_weight": 600, "logp": 2, "h_bond_donors": 1, "h_bond_acceptors": 3}
    assert mol_filter.apply_lipinski_filter(properties) is False

def test_apply_lipinski_filter_fail_logp():
    mol_filter = MoleculeFilter()
    properties = {"molecular_weight": 200, "logp": 6, "h_bond_donors": 1, "h_bond_acceptors": 3}
    assert mol_filter.apply_lipinski_filter(properties) is False

def test_apply_custom_filters_pass():
    mol_filter = MoleculeFilter()
    properties = {"tpsa": 80, "logp": 2.5}
    requirements = {"blood_brain_barrier_permeability": True}
    assert mol_filter.apply_custom_filters(properties, requirements) is True

def test_apply_custom_filters_fail_bbb():
    mol_filter = MoleculeFilter()
    properties = {"tpsa": 100, "logp": 2.5} # High TPSA
    requirements = {"blood_brain_barrier_permeability": True}
    assert mol_filter.apply_custom_filters(properties, requirements) is False

def test_filter_molecules():
    mol_filter = MoleculeFilter()
    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("c1ccccc1CCCCCCCCCC") # High MW, LogP
    
    mol1_props = {"molecular_weight": 46, "logp": 0.05, "h_bond_donors": 1, "h_bond_acceptors": 1, "tpsa": 20}
    mol2_props = {"molecular_weight": 200, "logp": 6.0, "h_bond_donors": 0, "h_bond_acceptors": 0, "tpsa": 0}

    molecules_with_properties = [
        {"mol": mol1, "properties": mol1_props, "smiles": "CCO"},
        {"mol": mol2, "properties": mol2_props, "smiles": "c1ccccc1CCCCCCCCCC"}
    ]
    requirements = {"blood_brain_barrier_permeability": True}

    filtered = mol_filter.filter_molecules(molecules_with_properties, requirements)
    assert len(filtered) == 1
    assert filtered[0]["smiles"] == "CCO"

# --- Tests for MoleculeGenerator ---
@patch('agents.molecule_generator.molecule_generator.MoleculeBuilder')
@patch('agents.molecule_generator.molecule_generator.PropertyCalculator')
@patch('agents.molecule_generator.molecule_generator.MoleculeFilter')
@patch('agents.molecule_generator.molecule_generator.os.makedirs')
@patch('agents.molecule_generator.molecule_generator.os.path.exists', return_value=True)
@patch('agents.molecule_generator.molecule_generator.json.load', return_value={}) # Mock memory load
def test_molecule_generator_process(mock_json_load, mock_exists, mock_makedirs, MockMoleculeFilter, MockPropertyCalculator, MockMoleculeBuilder, mock_gemini_client):
    # Mock instances of the components
    mock_builder_instance = MockMoleculeBuilder.return_value
    mock_calculator_instance = MockPropertyCalculator.return_value
    mock_filter_instance = MockMoleculeFilter.return_value

    # Configure mocks for the process flow
    mock_builder_instance.generate_smiles_from_binding_site.return_value = ["CCO"]
    mock_builder_instance.build_molecule_from_smiles.return_value = Chem.MolFromSmiles("CCO")
    mock_builder_instance.generate_3d_conformation.return_value = Chem.MolFromSmiles("CCO")
    mock_builder_instance.mol_to_sdf_string.return_value = "MOL_SDF_STRING"

    mock_calculator_instance.calculate_properties.return_value = {"molecular_weight": 46, "logp": 0.05}
    
    # Mock filter to pass all molecules for simplicity in this test
    mock_filter_instance.filter_molecules.side_effect = lambda mols, reqs: mols 

    generator = MoleculeGenerator(gemini_api_key="dummy_key", gemini_model="dummy_model", config={"data_dir": "/tmp/data"})
    generator.gemini_client = mock_gemini_client # Assign the fixture mock

    input_data = {
        "prepared_structures": {"KRAS": "/path/to/kras_prepared.json"},
        "binding_sites": {"KRAS": [{"name": "GTP_pocket", "center": {"x":0,"y":0,"z":0}, "radius":10}]},
        "parsed_query_data": {"molecule_properties": {"blood_brain_barrier_permeability": True}}
    }

    results = generator.process(input_data)

    assert "generated_molecules" in results
    assert "KRAS" in results["generated_molecules"]
    assert results["generated_molecules"]["KRAS"]["status"] == "molecules_generated"
    assert len(results["generated_molecules"]["KRAS"]["candidates"]) > 0
    assert results["generated_molecules"]["KRAS"]["candidates"][0]["smiles"] == "CCO"
    assert results["generated_molecules"]["KRAS"]["candidates"][0]["sdf_string"] == "MOL_SDF_STRING"

    mock_builder_instance.generate_smiles_from_binding_site.assert_called_once()
    mock_builder_instance.build_molecule_from_smiles.assert_called_once()
    mock_builder_instance.generate_3d_conformation.assert_called_once()
    mock_calculator_instance.calculate_properties.assert_called_once()
    mock_filter_instance.filter_molecules.assert_called_once()
    mock_builder_instance.mol_to_sdf_string.assert_called_once()

@patch('agents.molecule_generator.molecule_generator.os.makedirs')
@patch('agents.molecule_generator.molecule_generator.os.path.exists', return_value=True)
@patch('agents.molecule_generator.molecule_generator.json.load', return_value={}) # Mock memory load
def test_molecule_generator_process_no_binding_sites(mock_json_load, mock_exists, mock_makedirs, mock_gemini_client):
    generator = MoleculeGenerator(gemini_api_key="dummy_key", gemini_model="dummy_model", config={"data_dir": "/tmp/data"})
    generator.gemini_client = mock_gemini_client

    input_data = {
        "prepared_structures": {"KRAS": "/path/to/kras_prepared.json"},
        "binding_sites": {"KRAS": []}, # No binding sites
        "parsed_query_data": {"molecule_properties": {}}
    }

    results = generator.process(input_data)
    assert "generated_molecules" in results
    assert "KRAS" not in results["generated_molecules"] # KRAS should be skipped
    assert results["generated_molecules"] == {} # No molecules generated
