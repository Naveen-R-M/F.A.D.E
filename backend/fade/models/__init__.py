"""
Models package for F.A.D.E structured outputs.
"""

from .structured_outputs import (
    # Research models
    TargetExtractionOutput,
    QueryRefinementSuggestion,
    
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

from .guidance import (
    # Enhanced guidance models
    QueryRefinementSuggestion as EnhancedQueryRefinement,
    TargetClassification,
    SearchMetadata,
    UserGuidance
)

__all__ = [
    "TargetExtractionOutput",
    "QueryRefinementSuggestion",
    "StructureSelectionOutput",
    "PocketRankingOutput",
    "MoleculeDesignStrategy",
    "MoleculeFilteringOutput",
    "BindingPredictionOutput",
    "DrugCandidateAnalysis",
    "FinalReportSummary",
    "QueryClassification",
    "EnhancedQueryRefinement",
    "TargetClassification",
    "SearchMetadata",
    "UserGuidance"
]
