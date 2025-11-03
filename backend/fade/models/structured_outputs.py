"""
Pydantic models for structured LLM outputs in F.A.D.E pipeline.

This module defines all the structured output schemas used by LLM agents
throughout the drug discovery pipeline.
"""

from typing import Optional, List, Dict, Any, Literal
from pydantic import BaseModel, Field, model_validator


# ============================================================================
# Research Module Models
# ============================================================================

class TargetExtractionOutput(BaseModel):
    """
    Structured output for extracting protein target information from queries.
    Used in the research module for parsing natural language drug discovery requests.
    """
    gene_name: Optional[str] = Field(
        None, 
        description="Gene symbol (e.g., KRAS, EGFR, BTK, JAK2)"
    )
    protein_name: Optional[str] = Field(
        None, 
        description="Common protein name (focus on KINASES, PROTEASES, RECEPTORS)"
    )
    pdb_id: Optional[str] = Field(
        None,
        description="A 4-character PDB ID if one was explicitly mentioned (e.g., 4OBE, 8AFB)"
    )
    uniprot_id: Optional[str] = Field(
        None,
        description="A 6 or 10-character UniProt ID if one was explicitly mentioned (e.g., P01116, A0A1L1T3F0)"
    )
    mutations: List[str] = Field(
        default_factory=list,
        description="List of specific mutations (e.g., ['G12C', 'T790M', 'L858R'])"
    )
    drug_keywords: List[str] = Field(
        default_factory=list,
        description="Small molecule drug names mentioned (e.g., ['erlotinib', 'gefitinib'])"
    )
    target_type: Literal["kinase", "protease", "receptor", "enzyme", "unknown"] = Field(
        "unknown",
        description="Type of drug target"
    )
    disease_context: Optional[str] = Field(
        None,
        description="Disease or condition mentioned in the query"
    )
    additional_requirements: Optional[str] = Field(
        None,
        description="Other requirements mentioned (BBB penetrant, oral bioavailability, etc.)"
    )
    
    @model_validator(mode='after')
    def check_at_least_one_identifier(self) -> 'TargetExtractionOutput':
        """Ensure we have at least a gene name, protein name, PDB ID, or UniProt ID."""
        if not self.gene_name and not self.protein_name and not self.pdb_id and not self.uniprot_id:
            raise ValueError("Must have either gene_name, protein_name, pdb_id, or uniprot_id")
        return self

class QueryRefinementSuggestion(BaseModel):
    """
    Structured output for query refinement when no structures are found.
    Provides user-friendly guidance and alternative query suggestions.
    """
    explanation: str = Field(
        ...,
        description="Clear explanation of why no results were found (2-3 sentences)"
    )
    suggestions: List[str] = Field(
        ...,
        description="2-3 alternative query formulations that might work better",
        min_items=1,
        max_items=5
    )
    known_complexes: Optional[List[str]] = Field(
        None,
        description="Known complexes for this target if applicable"
    )
    target_characteristics: Optional[str] = Field(
        None,
        description="Brief note about the target type (e.g., 'peptidyl-prolyl isomerase with limited inhibitors')"
    )
    confidence: float = Field(
        ...,
        description="Confidence in the suggestions (0-1)",
        ge=0.0,
        le=1.0
    )


# ============================================================================
# Structure Module Models
# ============================================================================

class StructureSelectionOutput(BaseModel):
    """
    Structured output for selecting the best structure from multiple options.
    Used in the structure module for intelligent structure prioritization.
    """
    selected_pdb_id: Optional[str] = Field(
        None,
        description="PDB ID of the selected structure"
    )
    selection_reason: str = Field(
        ...,
        description="Reason for selecting this structure"
    )
    structure_quality: Literal["high", "medium", "low"] = Field(
        ...,
        description="Quality assessment of the selected structure"
    )
    has_relevant_ligand: bool = Field(
        False,
        description="Whether the structure has a relevant drug-like ligand"
    )
    ligand_similarity_to_query: Optional[float] = Field(
        None,
        description="Similarity score (0-1) of bound ligand to query requirements",
        ge=0.0,
        le=1.0
    )
    recommended_action: Literal["use_as_is", "remove_ligand", "predict_new", "use_alphafold"] = Field(
        ...,
        description="Recommended action for this structure"
    )


# ============================================================================
# Targeting Module Models
# ============================================================================

class PocketRankingOutput(BaseModel):
    """
    Structured output for ranking and selecting binding pockets.
    Used in the targeting module for intelligent pocket selection.
    """
    selected_pocket_id: str = Field(
        ...,
        description="ID of the selected pocket"
    )
    ranking_rationale: str = Field(
        ...,
        description="Explanation for why this pocket was selected"
    )
    druggability_assessment: Literal["excellent", "good", "moderate", "poor"] = Field(
        ...,
        description="Overall druggability assessment"
    )
    key_residues: List[str] = Field(
        default_factory=list,
        description="Key residues for drug binding in this pocket"
    )
    pocket_characteristics: Dict[str, Any] = Field(
        default_factory=dict,
        description="Important characteristics (hydrophobic, polar, size, etc.)"
    )
    confidence_score: float = Field(
        ...,
        description="Confidence in this pocket selection (0-1)",
        ge=0.0,
        le=1.0
    )


# ============================================================================
# Invention Module Models
# ============================================================================

class MoleculeDesignStrategy(BaseModel):
    """
    Structured output for molecule design strategy.
    Used to guide molecule generation based on pocket and target characteristics.
    """
    design_approach: Literal["fragment_based", "scaffold_hopping", "de_novo", "analog_design"] = Field(
        ...,
        description="Overall design strategy"
    )
    key_pharmacophores: List[str] = Field(
        default_factory=list,
        description="Key pharmacophore features to include"
    )
    avoid_features: List[str] = Field(
        default_factory=list,
        description="Chemical features to avoid"
    )
    molecular_weight_range: Dict[str, float] = Field(
        default={"min": 250, "max": 500},
        description="Target molecular weight range"
    )
    logp_range: Dict[str, float] = Field(
        default={"min": 1.0, "max": 3.5},
        description="Target LogP range for good drug-likeness"
    )
    specific_interactions: List[str] = Field(
        default_factory=list,
        description="Specific interactions to target (H-bonds with Ser139, π-stack with Tyr32, etc.)"
    )


class MoleculeFilteringOutput(BaseModel):
    """
    Structured output for molecule filtering decisions.
    Used to assess and filter generated molecules.
    """
    molecule_id: str = Field(
        ...,
        description="Identifier for the molecule"
    )
    passed_filters: bool = Field(
        ...,
        description="Whether the molecule passed all filters"
    )
    lipinski_violations: int = Field(
        0,
        description="Number of Lipinski's Rule of Five violations",
        ge=0,
        le=5
    )
    qed_score: float = Field(
        ...,
        description="Quantitative Estimate of Drug-likeness (0-1)",
        ge=0.0,
        le=1.0
    )
    synthetic_accessibility: float = Field(
        ...,
        description="Synthetic accessibility score (1=easy, 10=hard)",
        ge=1.0,
        le=10.0
    )
    warnings: List[str] = Field(
        default_factory=list,
        description="List of warnings or concerns about the molecule"
    )
    strengths: List[str] = Field(
        default_factory=list,
        description="Positive attributes of the molecule"
    )


# ============================================================================
# Screening Module Models
# ============================================================================

class BindingPredictionOutput(BaseModel):
    """
    Structured output for binding affinity predictions.
    Used in the screening module to assess molecule-target interactions.
    """
    molecule_id: str = Field(
        ...,
        description="Identifier for the molecule"
    )
    predicted_affinity: float = Field(
        ...,
        description="Predicted binding affinity (kcal/mol)"
    )
    binding_mode: Literal["competitive", "allosteric", "covalent", "unknown"] = Field(
        ...,
        description="Predicted binding mode"
    )
    key_interactions: List[str] = Field(
        default_factory=list,
        description="Key interactions identified"
    )
    pose_confidence: float = Field(
        ...,
        description="Confidence in the predicted pose (0-1)",
        ge=0.0,
        le=1.0
    )
    selectivity_prediction: Optional[str] = Field(
        None,
        description="Predicted selectivity profile"
    )


# ============================================================================
# Analysis Module Models
# ============================================================================

class DrugCandidateAnalysis(BaseModel):
    """
    Structured output for final drug candidate analysis.
    Used to generate comprehensive assessment of top candidates.
    """
    candidate_rank: int = Field(
        ...,
        description="Overall rank among all candidates",
        ge=1
    )
    molecule_id: str = Field(
        ...,
        description="Identifier for the molecule"
    )
    overall_score: float = Field(
        ...,
        description="Composite score combining all factors (0-100)",
        ge=0.0,
        le=100.0
    )
    strengths: List[str] = Field(
        default_factory=list,
        description="Key strengths of this candidate"
    )
    weaknesses: List[str] = Field(
        default_factory=list,
        description="Potential concerns or weaknesses"
    )
    optimization_suggestions: List[str] = Field(
        default_factory=list,
        description="Suggestions for further optimization"
    )
    development_potential: Literal["high", "medium", "low"] = Field(
        ...,
        description="Overall development potential"
    )
    recommended_next_steps: List[str] = Field(
        default_factory=list,
        description="Recommended next experimental or computational steps"
    )


class FinalReportSummary(BaseModel):
    """
    Structured output for the final natural language report.
    Used to structure the final report generation.
    """
    executive_summary: str = Field(
        ...,
        description="Brief executive summary of the drug discovery results"
    )
    target_summary: str = Field(
        ...,
        description="Summary of the target and its importance"
    )
    methodology_highlights: List[str] = Field(
        default_factory=list,
        description="Key methodological approaches used"
    )
    top_candidates: List[Dict[str, Any]] = Field(
        default_factory=list,
        description="Summary of top 3-5 candidates"
    )
    key_findings: List[str] = Field(
        default_factory=list,
        description="Most important findings"
    )
    limitations: List[str] = Field(
        default_factory=list,
        description="Study limitations to acknowledge"
    )
    future_directions: List[str] = Field(
        default_factory=list,
        description="Suggested future work"
    )


# ============================================================================
# Chat/Query Processing Models
# ============================================================================

class QueryClassification(BaseModel):
    """
    Structured output for classifying user queries.
    Used to route queries to appropriate handlers.
    """
    query_type: Literal["drug_discovery", "information", "clarification", "other"] = Field(
        ...,
        description="Type of query"
    )
    requires_computation: bool = Field(
        ...,
        description="Whether this query requires computational pipeline"
    )
    confidence: float = Field(
        ...,
        description="Confidence in classification (0-1)",
        ge=0.0,
        le=1.0
    )
    suggested_action: str = Field(
        ...,
        description="Suggested action to take"
    )
    extracted_parameters: Optional[Dict[str, Any]] = Field(
        None,
        description="Any parameters extracted from the query"
    )
