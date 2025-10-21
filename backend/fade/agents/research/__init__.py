"""
Research nodes for LangGraph pipeline.

These are LangGraph node functions, not agent classes.
"""

from fade.agents.research.langgraph_nodes import (
    target_research_node,
    validate_target_node
)

__all__ = [
    "target_research_node",
    "validate_target_node",
]
