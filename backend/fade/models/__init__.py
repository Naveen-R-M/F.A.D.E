"""
Models package for F.A.D.E structured outputs.
"""

from .structured_outputs import (
    # Research models
    TargetExtractionOutput,
    
    # Structure models
    StructureSelectionOutput,
    
    # Targeting models
    PocketRankingOutput,
    
    # Invention models
    MoleculeDesignStrategy,
    MoleculeFilteringOutput,
    
    # Screening models
    BindingPredictionOutput,
    
    # Analysis models
    DrugCandidateAnalysis,
    FinalReportSummary,
    
    # Query processing models
    QueryClassification
)

__all__ = [
    "TargetExtractionOutput",
    "StructureSelectionOutput",
    "PocketRankingOutput",
    "MoleculeDesignStrategy",
    "MoleculeFilteringOutput",
    "BindingPredictionOutput",
    "DrugCandidateAnalysis",
    "FinalReportSummary",
    "QueryClassification"
]
