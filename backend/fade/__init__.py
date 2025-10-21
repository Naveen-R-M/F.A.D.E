"""
F.A.D.E (Fully Agentic Drug Engine) - LangGraph Implementation

A comprehensive AI-powered drug discovery platform using LangGraph for orchestration.
This implementation uses LangGraph nodes instead of custom agent classes.
"""

__version__ = "2.0.0"
__author__ = "F.A.D.E Team"

# Core imports
from fade.state.langgraph_state import DrugDiscoveryState
from fade.workflows.drug_discovery import (
    create_drug_discovery_graph,
    run_drug_discovery_pipeline
)

__all__ = [
    "DrugDiscoveryState",
    "create_drug_discovery_graph",
    "run_drug_discovery_pipeline",
]
