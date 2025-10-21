"""
Drug Discovery State for LangGraph workflow.

This module defines the state that flows through the LangGraph nodes.
"""

from typing import TypedDict, Optional, List, Dict, Any, Annotated
from datetime import datetime
from langgraph.graph.message import add_messages


class DrugDiscoveryState(TypedDict):
    """
    Central state for the drug discovery pipeline using LangGraph.
    This state flows through all nodes in the graph.
    """
    # Metadata
    run_id: str
    job_id: Optional[str]  # Shared job ID for HPC operations
    timestamp: datetime
    user_id: Optional[str]
    
    # Input
    query: str
    
    # Messages for conversation history (LangGraph pattern)
    messages: Annotated[list, add_messages]
    
    # Target Information (from Research Module)
    target_info: Optional[Dict[str, Any]]
    known_compounds: Optional[List[Dict[str, Any]]]
    
    # Structure Information (from Structure Module)
    structure: Optional[Dict[str, Any]]
    structure_source: Optional[str]  # "PDB", "AlphaFold", "Boltz2"
    
    # Pocket Information (from Targeting Module)
    pockets: Optional[List[Dict[str, Any]]]
    selected_pocket: Optional[Dict[str, Any]]
    
    # Generated Molecules (from Invention Module)
    generated_molecules: Optional[List[Dict[str, Any]]]
    filtered_molecules: Optional[List[Dict[str, Any]]]
    
    # Screening Results (from Screening Module)
    screening_results: Optional[List[Dict[str, Any]]]
    
    # Final Analysis (from Analysis Module)
    analysis: Optional[Dict[str, Any]]
    final_report: Optional[str]
    
    # Pipeline Control
    current_step: Optional[str]
    error: Optional[str]
    should_continue: bool
