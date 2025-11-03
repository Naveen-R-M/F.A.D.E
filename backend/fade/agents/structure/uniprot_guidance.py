"""
UniProt failure guidance generator using LLM.

This module generates helpful guidance when UniProt lookups fail,
explaining the issue and suggesting alternative queries.
"""

from typing import Dict, Any, List, Optional
from langchain_core.output_parsers import PydanticOutputParser
from pydantic import BaseModel, Field
from fade.tools.llm_client import get_llm_client
from fade.utils import get_logger

logger = get_logger("guidance.uniprot")


class UniProtGuidance(BaseModel):
    """Structured guidance for UniProt lookup failures."""
    
    explanation: str = Field(
        description="Clear explanation of why the UniProt lookup failed"
    )
    
    suggestions: List[str] = Field(
        description="3-5 specific alternative queries the user can try",
        min_items=3,
        max_items=5
    )
    
    tips: List[str] = Field(
        description="General tips for successful UniProt queries",
        min_items=2,
        max_items=4
    )
    
    confidence: float = Field(
        description="Confidence score (0-1) in the suggestions",
        ge=0.0,
        le=1.0
    )


def generate_uniprot_guidance(
    query: str,
    uniprot_id: str,
    error_type: str = "not_found"
) -> Dict[str, Any]:
    """
    Generate intelligent guidance for UniProt lookup failures.
    
    Args:
        query: Original user query
        uniprot_id: The UniProt ID that failed
        error_type: Type of error ("not_found", "invalid_format", "no_sequence")
        
    Returns:
        Dict with guidance information
    """
    logger.info(f"Generating UniProt guidance for failed ID: {uniprot_id}")
    
    try:
        llm = get_llm_client()
        
        # Create structured LLM for guidance
        structured_llm = llm.with_structured_output(UniProtGuidance)
        
        # Build context-aware prompt
        if error_type == "not_found":
            error_context = f"""
The UniProt ID '{uniprot_id}' was not found in the UniProt database.
This could mean:
1. The ID is outdated or has been deprecated
2. It's a typo or incorrect format
3. It's from a different database (not UniProt)
4. It's a protein from a non-standard organism
"""
        elif error_type == "no_sequence":
            error_context = f"""
The UniProt entry '{uniprot_id}' exists but has no sequence data.
This is unusual and might indicate:
1. A deprecated or merged entry
2. A fragment or partial sequence entry
3. A predicted protein that hasn't been validated
"""
        else:
            error_context = f"The UniProt ID '{uniprot_id}' has an issue: {error_type}"
        
        prompt = f"""
A user tried to query: "{query}"

{error_context}

Generate helpful guidance for the user. Consider:
- The ID might be from a different database (PDB, RefSeq, Ensembl)
- Suggest using gene names or protein names instead
- Recommend well-studied proteins if they want examples
- For A0A prefixes, these are TrEMBL (unreviewed) entries that might not exist

Provide:
1. A clear, friendly explanation of what went wrong
2. 3-5 specific alternative queries they can try
3. 2-4 general tips for successful queries
4. Your confidence in the suggestions (0-1)

Make suggestions specific and actionable. Include example UniProt IDs or gene names they can try.
"""
        
        # Get structured guidance
        result = structured_llm.invoke(prompt)
        
        # Format the guidance
        guidance = {
            "explanation": result.explanation,
            "suggestions": result.suggestions,
            "tips": result.tips,
            "confidence": result.confidence,
            "original_query": query,
            "failed_id": uniprot_id
        }
        
        logger.info(f"Generated {len(result.suggestions)} suggestions with {result.confidence:.0%} confidence")
        return guidance
        
    except Exception as e:
        logger.error(f"Failed to generate UniProt guidance: {e}")
        
        # Provide fallback guidance
        return {
            "explanation": f"The UniProt ID '{uniprot_id}' could not be found in the database.",
            "suggestions": [
                "Try searching with a gene name like 'EGFR' or 'KRAS'",
                "Use a validated UniProt ID like 'P00533' (EGFR) or 'P01116' (KRAS)",
                "Search for a protein name like 'epidermal growth factor receptor'",
                f"Verify the ID '{uniprot_id}' is correct and from UniProt (not PDB or other databases)"
            ],
            "tips": [
                "UniProt IDs starting with 'A0A' are unreviewed and may not always exist",
                "Gene names (like EGFR, BRCA1) are more reliable than protein IDs",
                "You can check valid IDs at uniprot.org"
            ],
            "confidence": 0.5,
            "original_query": query,
            "failed_id": uniprot_id
        }


def format_uniprot_guidance_message(guidance: Dict[str, Any]) -> str:
    """
    Format UniProt guidance into a user-friendly message.
    
    Args:
        guidance: Guidance dictionary from generate_uniprot_guidance
        
    Returns:
        Formatted message string
    """
    lines = []
    
    # Add explanation
    lines.append("💡 UNIPROT LOOKUP GUIDANCE")
    lines.append("="*60)
    lines.append(f"\n{guidance['explanation']}")
    
    # Add failed ID info
    lines.append(f"\nFailed ID: {guidance['failed_id']}")
    lines.append(f"Original query: {guidance['original_query']}")
    
    # Add suggestions
    lines.append("\n📝 Try one of these alternative queries:")
    for i, suggestion in enumerate(guidance['suggestions'], 1):
        # Make suggestions directly executable
        if any(keyword in suggestion.lower() for keyword in ['try', 'use', 'search']):
            # Extract the actual query from the suggestion
            if "'" in suggestion or '"' in suggestion:
                import re
                match = re.search(r"['\"]([^'\"]+)['\"]", suggestion)
                if match:
                    query_part = match.group(1)
                    lines.append(f"   {i}. python main.py \"{query_part}\"")
                else:
                    lines.append(f"   {i}. {suggestion}")
            else:
                lines.append(f"   {i}. {suggestion}")
        else:
            lines.append(f"   {i}. python main.py \"{suggestion}\"")
    
    # Add tips
    lines.append("\n💭 General Tips:")
    for tip in guidance['tips']:
        lines.append(f"   • {tip}")
    
    # Add confidence if high
    if guidance['confidence'] >= 0.7:
        lines.append(f"\n🎯 Confidence: {guidance['confidence']:.0%}")
    
    return "\n".join(lines)
