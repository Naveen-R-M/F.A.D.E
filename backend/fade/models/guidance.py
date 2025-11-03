"""
Pydantic models for structured guidance when queries fail or need refinement.
"""

from typing import List, Optional, Literal
from pydantic import BaseModel, Field


class QueryRefinementSuggestion(BaseModel):
    """Structured output for query refinement when searches fail."""
    
    problem_analysis: str = Field(
        ..., 
        description="Clear explanation of why the search failed"
    )
    
    target_context: str = Field(
        ..., 
        description="Information about this target type and why it might be challenging"
    )
    
    suggestions: List[str] = Field(
        ..., 
        max_items=5,
        min_items=1,
        description="Specific query reformulations that might work better"
    )
    
    known_alternatives: Optional[List[str]] = Field(
        None,
        description="Known complexes or alternative targets if applicable"
    )
    
    difficulty_level: Literal["easy", "moderate", "challenging"] = Field(
        ...,
        description="How challenging this target is for drug discovery"
    )
    
    next_steps: str = Field(
        ...,
        description="What the user should do next"
    )
    
    confidence: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Confidence in the suggestions (0-1)"
    )


class TargetClassification(BaseModel):
    """Classification of the target type for specialized guidance."""
    
    target_type: Literal["kinase", "enzyme", "receptor", "transcription_factor", "other"] = Field(
        ...,
        description="The broad class of the protein target"
    )
    
    druggability: Literal["high", "moderate", "low", "undruggable"] = Field(
        ...,
        description="Assessment of how druggable this target is"
    )
    
    common_issues: List[str] = Field(
        ...,
        description="Common issues when searching for this target type"
    )
    
    recommended_databases: List[str] = Field(
        default=["PDB", "ChEMBL", "UniProt"],
        description="Best databases for this target type"
    )


class SearchMetadata(BaseModel):
    """Metadata about failed search attempts."""
    
    original_query: str = Field(
        ...,
        description="The original user query"
    )
    
    search_attempts: List[dict] = Field(
        default_factory=list,
        description="Details of each search attempt made"
    )
    
    pdb_ids_checked: List[str] = Field(
        default_factory=list,
        description="PDB IDs that were checked but didn't have small molecules"
    )
    
    failure_reason: Literal[
        "no_results",
        "no_small_molecules", 
        "ambiguous_target",
        "missing_target_info",
        "invalid_query"
    ] = Field(
        ...,
        description="The specific reason for failure"
    )


class UserGuidance(BaseModel):
    """Complete guidance package for the user."""
    
    refinement: QueryRefinementSuggestion = Field(
        ...,
        description="Structured refinement suggestions"
    )
    
    target_classification: Optional[TargetClassification] = Field(
        None,
        description="Classification of the target if identified"
    )
    
    search_metadata: SearchMetadata = Field(
        ...,
        description="Metadata about what was searched"
    )
    
    user_message: str = Field(
        ...,
        description="Formatted message for display to user"
    )
    
    retry_queries: List[str] = Field(
        ...,
        max_items=3,
        description="Top queries user can immediately retry"
    )
