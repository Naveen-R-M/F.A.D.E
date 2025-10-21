"""
Central State Definition for Drug Discovery Pipeline

This module defines the TypedDict that represents the complete state
of the drug discovery pipeline as it flows through different agents.
"""

from typing import TypedDict, Optional, List, Dict, Any, Literal
from datetime import datetime


class ProteinTarget(TypedDict):
    """Information about the protein target"""
    uniprot_id: Optional[str]
    protein_name: Optional[str]
    gene_name: Optional[str]
    organism: Optional[str]
    sequence: Optional[str]
    sequence_length: Optional[int]
    function_description: Optional[str]
    disease_associations: Optional[List[str]]
    mutations: Optional[List[str]]  # e.g., ["G12C", "G12D"]


class KnownCompound(TypedDict):
    """Information about known compounds/drugs"""
    compound_id: str  # ChEMBL ID or PDB ligand ID
    name: Optional[str]
    smiles: str
    binding_affinity: Optional[float]  # in nM or µM
    affinity_unit: Optional[str]
    source: str  # "ChEMBL", "PDB", "DrugBank", etc.
    clinical_phase: Optional[str]  # "Preclinical", "Phase I", "Phase II", "Phase III", "Approved"
    mechanism: Optional[str]


class ProteinStructure(TypedDict):
    """Information about the protein structure"""
    structure_path: str  # Path to PDB file
    pdb_id: Optional[str]
    resolution: Optional[float]  # in Angstroms
    source: Literal["PDB", "AlphaFold", "Boltz2"]
    confidence_score: Optional[float]  # pLDDT for AlphaFold/Boltz2
    chain_id: Optional[str]
    has_ligand: bool
    ligand_names: Optional[List[str]]
    prepared_structure_path: Optional[str]  # Path to prepared/cleaned structure


class Pocket(TypedDict):
    """Information about a binding pocket"""
    pocket_id: str
    center: List[float]  # [x, y, z] coordinates
    residues: List[str]  # List of residue IDs
    volume: float  # in Å³
    surface_area: float  # in Ų
    druggability_score: float  # 0-1 scale
    depth: float  # in Å
    hydrophobicity: float
    electrostatic_potential: Optional[float]
    is_known_site: bool  # Whether this overlaps with known binding sites


class GeneratedMolecule(TypedDict):
    """Information about a generated molecule"""
    molecule_id: str
    smiles: str
    generation_method: str  # "DiffSBDD", "REINVENT", "etc"
    mol_weight: float
    logp: float
    tpsa: float  # Topological Polar Surface Area
    hbd: int  # Hydrogen Bond Donors
    hba: int  # Hydrogen Bond Acceptors
    rotatable_bonds: int
    qed_score: float  # Quantitative Estimate of Drug-likeness
    sa_score: float  # Synthetic Accessibility
    lipinski_violations: int
    pains_alerts: List[str]
    passed_filters: bool
    sdf_path: Optional[str]  # Path to 3D structure file


class ScreeningResult(TypedDict):
    """Results from Boltz2 screening"""
    molecule_id: str
    smiles: str
    binding_affinity: float  # Predicted binding affinity (kcal/mol)
    affinity_confidence: float  # Confidence score from Boltz2
    pose_path: Optional[str]  # Path to predicted pose PDB file
    interactions: Optional[Dict[str, List[str]]]  # e.g., {"H-bonds": ["ASP189", "SER195"]}
    rmsd_from_pocket: Optional[float]  # RMSD from pocket center
    strain_energy: Optional[float]


class AnalysisResult(TypedDict):
    """Final analysis results"""
    top_molecule: Optional[ScreeningResult]
    top_10_molecules: List[ScreeningResult]
    comparison_with_known: Dict[str, Any]
    novel_scaffolds: List[str]
    key_insights: List[str]
    synthesis_feasibility: Dict[str, Any]
    report_path: Optional[str]  # Path to generated report
    visualization_paths: Optional[List[str]]  # Paths to generated images


class DrugDiscoveryState(TypedDict):
    """
    Central state for the drug discovery pipeline.
    This state flows through all agents in the LangGraph workflow.
    """
    # Metadata
    run_id: str
    timestamp: datetime
    user_id: Optional[str]
    
    # Input
    query: str  # Original natural language query
    
    # Target Information (from Research Module)
    target_info: Optional[ProteinTarget]
    known_compounds: Optional[List[KnownCompound]]
    
    # Structure Information (from Structure Module)
    structure: Optional[ProteinStructure]
    structure_validation: Optional[Dict[str, Any]]
    
    # Pocket Information (from Targeting Module)
    pockets: Optional[List[Pocket]]
    selected_pocket: Optional[Pocket]
    pocket_selection_rationale: Optional[str]
    
    # Generated Molecules (from Invention Module)
    generated_molecules: Optional[List[GeneratedMolecule]]
    filtered_molecules: Optional[List[GeneratedMolecule]]
    generation_parameters: Optional[Dict[str, Any]]
    
    # Screening Results (from Screening Module)
    screening_results: Optional[List[ScreeningResult]]
    screening_metadata: Optional[Dict[str, Any]]
    
    # Final Analysis (from Analysis Module)
    analysis: Optional[AnalysisResult]
    final_report: Optional[str]  # Natural language report
    
    # Pipeline Control
    current_step: Optional[str]
    error_messages: Optional[List[str]]
    warnings: Optional[List[str]]
    execution_log: Optional[List[Dict[str, Any]]]
    should_continue: bool
    max_molecules_to_generate: int
    max_molecules_to_screen: int
    
    # Performance Metrics
    execution_times: Optional[Dict[str, float]]  # Time taken by each module
    api_calls_count: Optional[Dict[str, int]]  # API calls made
    total_cost: Optional[float]  # If using paid APIs


# Additional type hints for complex nested structures
class ExecutionLogEntry(TypedDict):
    """Entry in the execution log"""
    timestamp: datetime
    agent: str
    action: str
    status: Literal["started", "completed", "failed", "skipped"]
    message: Optional[str]
    data: Optional[Dict[str, Any]]


class StateValidation(TypedDict):
    """Validation result for state transitions"""
    is_valid: bool
    errors: List[str]
    warnings: List[str]
