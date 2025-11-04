"""
Drug Discovery State for LangGraph workflow.

This module defines the state that flows through the LangGraph nodes.
"""

# --- Import NotRequired ---
from typing import TypedDict, Optional, List, Dict, Any, Annotated, NotRequired
from datetime import datetime
from langgraph.graph.message import add_messages


class DrugDiscoveryState(TypedDict):
    """
    Central state for the drug discovery pipeline using LangGraph.
    This state flows through all nodes in the graph.
    """
    # --- Required fields (must be in initial_state) ---
    run_id: str
    timestamp: datetime
    query: str
    messages: Annotated[list, add_messages]
    should_continue: bool
    current_step: Optional[str] # Included as required since it's in initial_state
    
    # --- NotRequired fields (key can be omitted) ---
    job_id: NotRequired[Optional[str]]
    user_id: NotRequired[Optional[str]]
    
    target_info: NotRequired[Optional[Dict[str, Any]]]
    known_compounds: NotRequired[Optional[List[Dict[str, Any]]]]
    
    structure: NotRequired[Optional[Dict[str, Any]]]
    structure_source: NotRequired[Optional[str]]
    
    pockets: NotRequired[Optional[List[Dict[str, Any]]]]
    selected_pocket: NotRequired[Optional[Dict[str, Any]]]
    
    generated_molecules: NotRequired[Optional[List[Dict[str, Any]]]]
    filtered_molecules: NotRequired[Optional[List[Dict[str, Any]]]]
    
    screening_results: NotRequired[Optional[List[Dict[str, Any]]]]
    
    analysis: NotRequired[Optional[Dict[str, Any]]]
    final_report: NotRequired[Optional[str]]
    
    error: NotRequired[Optional[str]]
    is_continuation: NotRequired[Optional[bool]]  # Flag for continuation queries
    routing_decision: NotRequired[Optional[str]]  # For entry router to specify next node
    
    # For intelligent guidance on failures
    guidance: NotRequired[Optional[str]]
    suggested_queries: NotRequired[Optional[List[str]]]