"""
Guidance generator for creating helpful user guidance when searches fail.
Uses LLM to analyze failures and suggest improvements.
"""

import logging
from typing import Dict, Any, Optional, List
from langchain_core.language_models import BaseLanguageModel

from fade.models.guidance import (
    QueryRefinementSuggestion,
    TargetClassification,
    SearchMetadata,
    UserGuidance
)
from fade.utils import get_logger

logger = get_logger("agents.research.guidance_generator")


def generate_query_guidance(
    llm: BaseLanguageModel,
    original_query: str,
    target_info: Dict[str, Any],
    search_metadata: Dict[str, Any]
) -> UserGuidance:
    """
    Generate intelligent guidance when searches fail.
    
    Args:
        llm: Language model to use
        original_query: The original user query
        target_info: Extracted target information
        search_metadata: Metadata about the failed search
        
    Returns:
        UserGuidance object with structured guidance
    """
    logger.info(f"Generating guidance for failed query: {original_query}")
    
    # First, classify the target
    target_classification = _classify_target(llm, target_info)
    
    # Build search metadata
    search_meta = SearchMetadata(
        original_query=original_query,
        search_attempts=search_metadata.get("search_attempts", []),
        pdb_ids_checked=search_metadata.get("pdb_ids_checked", []),
        failure_reason=_determine_failure_reason(search_metadata)
    )
    
    # Generate refinement suggestions based on target type
    refinement = _generate_refinement_suggestions(
        llm, 
        original_query, 
        target_info, 
        target_classification,
        search_meta
    )
    
    # Format user-friendly message
    user_message = _format_user_message(refinement, target_classification)
    
    # Extract top retry queries
    retry_queries = refinement.suggestions[:3]
    
    return UserGuidance(
        refinement=refinement,
        target_classification=target_classification,
        search_metadata=search_meta,
        user_message=user_message,
        retry_queries=retry_queries
    )


def _classify_target(llm: BaseLanguageModel, target_info: Dict[str, Any]) -> Optional[TargetClassification]:
    """
    Classify the target type for specialized guidance.
    
    Args:
        llm: Language model to use
        target_info: Extracted target information
        
    Returns:
        TargetClassification object or None
    """
    if not target_info or not target_info.get("protein_name"):
        return None
    
    prompt = f"""
    Classify this protein target for drug discovery:
    
    Protein: {target_info.get('protein_name', 'Unknown')}
    Gene: {target_info.get('gene_name', 'Unknown')}
    Type hint: {target_info.get('target_type', 'Unknown')}
    
    Classify into one of: kinase, enzyme, receptor, transcription_factor, other
    Assess druggability: high, moderate, low, undruggable
    List 2-3 common issues when searching for this target type.
    """
    
    try:
        classification = llm.with_structured_output(TargetClassification).invoke(prompt)
        logger.debug(f"Target classified as: {classification.target_type} with {classification.druggability} druggability")
        return classification
    except Exception as e:
        logger.warning(f"Could not classify target: {e}")
        return None


def _determine_failure_reason(metadata: Dict[str, Any]) -> str:
    """
    Determine the specific reason for failure.
    
    Args:
        metadata: Search metadata
        
    Returns:
        Failure reason string
    """
    if metadata.get("no_results"):
        if metadata.get("no_small_molecules"):
            return "no_small_molecules"
        else:
            return "no_results"
    elif metadata.get("partial_results"):
        return "missing_target_info"
    else:
        return "invalid_query"


def _generate_refinement_suggestions(
    llm: BaseLanguageModel,
    original_query: str,
    target_info: Dict[str, Any],
    target_classification: Optional[TargetClassification],
    search_metadata: SearchMetadata
) -> QueryRefinementSuggestion:
    """
    Generate specific refinement suggestions based on target type.
    
    Args:
        llm: Language model to use
        original_query: The original user query
        target_info: Extracted target information
        target_classification: Target classification if available
        search_metadata: Search metadata
        
    Returns:
        QueryRefinementSuggestion with structured suggestions
    """
    # Build context-aware prompt based on target type
    if target_classification:
        prompt = _build_specialized_prompt(
            original_query,
            target_info,
            target_classification,
            search_metadata
        )
    else:
        prompt = _build_general_prompt(
            original_query,
            target_info,
            search_metadata
        )
    
    # Use structured output to get suggestions
    try:
        suggestions = llm.with_structured_output(QueryRefinementSuggestion).invoke(prompt)
        logger.info(f"Generated {len(suggestions.suggestions)} refinement suggestions")
        return suggestions
    except Exception as e:
        logger.error(f"Error generating refinement suggestions: {e}")
        # Return fallback suggestions
        return QueryRefinementSuggestion(
            problem_analysis="The search for protein structures with small molecule inhibitors did not return results.",
            target_context="This target may have limited structural data available.",
            suggestions=[
                f"{target_info.get('gene_name', 'target')} with known inhibitor",
                f"{target_info.get('protein_name', 'protein')} crystal structure",
                f"{target_info.get('gene_name', 'target')} small molecule complex"
            ],
            difficulty_level="moderate",
            next_steps="Try one of the suggested queries or provide more specific inhibitor names.",
            confidence=0.5
        )


def _build_specialized_prompt(
    original_query: str,
    target_info: Dict[str, Any],
    target_classification: TargetClassification,
    search_metadata: SearchMetadata
) -> str:
    """
    Build a specialized prompt based on target type.
    
    Args:
        original_query: Original user query
        target_info: Target information
        target_classification: Target classification
        search_metadata: Search metadata
        
    Returns:
        Formatted prompt string
    """
    target_type = target_classification.target_type
    
    # Specialized prompts for different target types
    if target_type == "enzyme" and "prolyl" in target_info.get("protein_name", "").lower():
        # Special case for prolyl isomerases like Cyclophilin
        return f"""
        The user searched for: "{original_query}"
        
        Target identified: {target_info.get('protein_name')} ({target_info.get('gene_name')})
        This is a peptidyl-prolyl isomerase enzyme.
        
        No small molecule inhibitors were found in PDB structures.
        
        IMPORTANT: Cyclophilins (like PPIA) have very few small molecule inhibitors.
        The main inhibitors are cyclosporin A (CsA) and sanglifehrin A.
        
        Generate helpful guidance that:
        1. Explains that this enzyme class has limited small molecule inhibitors
        2. Suggests searching for known complexes with cyclosporin or sanglifehrin
        3. Provides 2-3 specific query formulations
        
        For PPIA/Cyclophilin A specifically, suggest:
        - "Cyclophilin A cyclosporin complex"
        - "PPIA CsA structure"
        - "Human cyclophilin sanglifehrin"
        
        Return a structured response with clear explanation and actionable suggestions.
        """
    
    elif target_type == "kinase":
        return f"""
        The user searched for: "{original_query}"
        
        Target identified: {target_info.get('protein_name')} ({target_info.get('gene_name')})
        This is a kinase target.
        
        No structures with small molecule inhibitors were found.
        
        Kinases usually have many inhibitor structures available.
        The issue might be:
        1. Using an uncommon name or synonym
        2. Missing the specific mutation if applicable
        3. The specific kinase might be less studied
        
        Generate helpful guidance with 3-4 specific query suggestions.
        Include mutation if mentioned: {target_info.get('mutations', [])}
        
        Return a structured response focusing on kinase-specific search strategies.
        """
    
    else:
        # Generic specialized prompt
        return f"""
        The user searched for: "{original_query}"
        
        Target identified: {target_info.get('protein_name')} ({target_info.get('gene_name')})
        Target type: {target_type}
        Druggability: {target_classification.druggability}
        
        No structures with small molecule inhibitors were found.
        Common issues for this target type: {target_classification.common_issues}
        
        Generate helpful guidance that:
        1. Explains why this specific target type might lack structures
        2. Suggests alternative query formulations
        3. Mentions any known inhibitors if applicable
        
        Return a structured response with actionable suggestions.
        """


def _build_general_prompt(
    original_query: str,
    target_info: Dict[str, Any],
    search_metadata: SearchMetadata
) -> str:
    """
    Build a general prompt when target classification is not available.
    
    Args:
        original_query: Original user query
        target_info: Target information
        search_metadata: Search metadata
        
    Returns:
        Formatted prompt string
    """
    pdb_checked = len(search_metadata.pdb_ids_checked)
    
    return f"""
    The user searched for: "{original_query}"
    
    Target information extracted:
    - Protein: {target_info.get('protein_name', 'Unknown')}
    - Gene: {target_info.get('gene_name', 'Unknown')}
    - Type: {target_info.get('target_type', 'Unknown')}
    
    Search result: {search_metadata.failure_reason}
    {f"Checked {pdb_checked} PDB structures but none had drug-like small molecules" if pdb_checked > 0 else "No PDB structures found"}
    
    Generate helpful guidance that:
    1. Explains why the search might have failed
    2. Suggests 2-3 specific alternative queries
    3. Provides actionable next steps
    
    Be concise but helpful. Focus on what the user can do differently.
    
    Return a structured response with clear suggestions.
    """


def _format_user_message(
    refinement: QueryRefinementSuggestion,
    target_classification: Optional[TargetClassification]
) -> str:
    """
    Format a user-friendly message from the refinement suggestions.
    
    Args:
        refinement: Refinement suggestions
        target_classification: Target classification if available
        
    Returns:
        Formatted user message string
    """
    lines = []
    
    # Add problem analysis
    lines.append(f"💡 {refinement.problem_analysis}")
    lines.append("")
    
    # Add target context if available
    if refinement.target_context:
        lines.append(f"ℹ️ {refinement.target_context}")
        lines.append("")
    
    # Add suggestions
    lines.append("**Suggested queries:**")
    for suggestion in refinement.suggestions[:3]:
        lines.append(f"  • {suggestion}")
    
    # Add known complexes if available
    if refinement.known_alternatives:
        lines.append("")
        lines.append("**Known complexes:**")
        for complex_name in refinement.known_alternatives[:3]:
            lines.append(f"  • {complex_name}")
    
    # Add next steps
    lines.append("")
    lines.append(f"**Next steps:** {refinement.next_steps}")
    
    # Add difficulty indicator if classification available
    if target_classification:
        lines.append("")
        difficulty_emoji = {
            "challenging": "🔴",
            "moderate": "🟡",
            "easy": "🟢"
        }.get(refinement.difficulty_level, "⚪")
        lines.append(f"{difficulty_emoji} Target difficulty: {refinement.difficulty_level}")
    
    return "\n".join(lines)
