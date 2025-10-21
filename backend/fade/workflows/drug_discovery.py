"""
Main LangGraph workflow for drug discovery pipeline.

This module defines the StateGraph that orchestrates all agents.
"""

import uuid
from datetime import datetime
from typing import Dict, Any

from langgraph.graph import StateGraph, END
from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.agents.research.langgraph_nodes import (
    target_research_node, 
    validate_target_node
)
from fade.agents.structure.langgraph_nodes import (
    structure_resolver_node,
    structure_wait_node
)
from fade.agents.targeting.langgraph_nodes import (
    pocket_mapper_node,
    pocket_ranker_node
)
from fade.agents.invention.langgraph_nodes import (
    molecule_generator_node,
    medicinal_filter_node
)
from fade.config import config
from fade.utils import get_logger

logger = get_logger("workflow.drug_discovery")


def create_drug_discovery_graph():
    """
    Create the LangGraph workflow for drug discovery.
    
    Returns:
        Compiled LangGraph application
    """
    # Initialize the graph with our state
    workflow = StateGraph(DrugDiscoveryState)
    
    # Add nodes for each step
    workflow.add_node("target_research", target_research_node)
    workflow.add_node("validate_target", validate_target_node)
    workflow.add_node("structure_resolver", structure_resolver_node)
    workflow.add_node("structure_wait", structure_wait_node)
    workflow.add_node("pocket_mapper", pocket_mapper_node)
    workflow.add_node("pocket_ranker", pocket_ranker_node)
    workflow.add_node("molecule_generator", molecule_generator_node)
    workflow.add_node("medicinal_filter", medicinal_filter_node)
    
    # TODO: Add these nodes when implemented
    # workflow.add_node("boltz2_screener", boltz2_screener_node)
    # workflow.add_node("results_analyzer", results_analyzer_node)
    
    # Define the edges (workflow sequence)
    workflow.set_entry_point("target_research")
    
    # Research -> Validation
    workflow.add_conditional_edges(
        "target_research",
        should_continue_research,
        {
            "continue": "validate_target",
            "end": END
        }
    )
    
    # Validation -> Structure
    workflow.add_conditional_edges(
        "validate_target",
        should_continue_validation,
        {
            "continue": "structure_resolver",
            "end": END
        }
    )
    
    # Structure Resolver -> Wait (if Boltz-2) or Pocket Mapping
    workflow.add_conditional_edges(
        "structure_resolver",
        should_wait_for_structure,
        {
            "wait": "structure_wait",
            "continue": "pocket_mapper",
            "end": END
        }
    )
    
    # Structure Wait -> Pocket Mapping
    workflow.add_conditional_edges(
        "structure_wait",
        should_continue_structure,
        {
            "continue": "pocket_mapper",
            "end": END
        }
    )
    
    # Pocket Mapper -> Pocket Ranker
    workflow.add_conditional_edges(
        "pocket_mapper",
        should_continue_pockets,
        {
            "continue": "pocket_ranker",
            "end": END
        }
    )
    
    # Pocket Ranker -> Molecule Generator
    workflow.add_conditional_edges(
        "pocket_ranker",
        should_continue_targeting,
        {
            "continue": "molecule_generator",
            "end": END
        }
    )
    
    # Molecule Generator -> Medicinal Filter
    workflow.add_conditional_edges(
        "molecule_generator",
        should_continue_generation,
        {
            "continue": "medicinal_filter",
            "end": END
        }
    )
    
    # Medicinal Filter -> Next step (screening)
    workflow.add_conditional_edges(
        "medicinal_filter",
        should_continue_filtering,
        {
            "continue": END,  # Will be boltz2_screener when implemented
            "end": END
        }
    )
    
    # TODO: Add more edges as nodes are implemented
    # workflow.add_edge("medicinal_filter", "boltz2_screener")
    # workflow.add_edge("boltz2_screener", "results_analyzer")
    # workflow.add_edge("results_analyzer", END)
    
    # Compile the graph
    app = workflow.compile()
    
    return app


def should_continue_research(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after research."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Research failed: {state.get('error')}")
        return "end"
    
    if not state.get("target_info"):
        logger.warning("No target information found")
        return "end"
    
    return "continue"


def should_continue_validation(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after validation."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Validation failed: {state.get('error')}")
        return "end"
    
    return "continue"


def should_wait_for_structure(state: DrugDiscoveryState) -> str:
    """Determine next step after structure resolver."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Structure resolution failed: {state.get('error')}")
        return "end"
    
    structure = state.get("structure", {})
    
    # If Boltz-2 job is pending, wait for it
    if structure.get("source") == "Boltz2" and structure.get("status") == "pending":
        return "wait"
    
    # Otherwise continue to pocket mapping
    return "continue"


def should_continue_structure(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after structure wait."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Structure prediction failed: {state.get('error')}")
        return "end"
    
    return "continue"


def should_continue_pockets(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after pocket mapping."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Pocket detection failed: {state.get('error')}")
        return "end"
    
    if not state.get("pockets"):
        logger.warning("No pockets detected")
        return "end"
    
    return "continue"


def should_continue_targeting(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after pocket ranking."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Pocket selection failed: {state.get('error')}")
        return "end"
    
    if not state.get("selected_pocket"):
        logger.warning("No pocket selected")
        return "end"
    
    return "continue"


def should_continue_generation(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after molecule generation."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Molecule generation failed: {state.get('error')}")
        return "end"
    
    if not state.get("generated_molecules"):
        logger.warning("No molecules generated")
        return "end"
    
    return "continue"


def should_continue_filtering(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after medicinal filtering."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Molecule filtering failed: {state.get('error')}")
        return "end"
    
    if not state.get("filtered_molecules"):
        logger.warning("No molecules passed filters")
        return "end"
    
    return "continue"


def run_drug_discovery_pipeline(query: str, user_id: str = None) -> Dict[str, Any]:
    """
    Run the complete drug discovery pipeline.
    
    Args:
        query: Natural language drug discovery query
        user_id: Optional user identifier
        
    Returns:
        Final state with results
    """
    logger.info(f"Starting pipeline for query: {query}")
    
    # Create initial state
    run_id = str(uuid.uuid4())[:8]
    initial_state: DrugDiscoveryState = {
        "run_id": run_id,
        "job_id": None,  # Will be set by first HPC operation
        "timestamp": datetime.now(),
        "user_id": user_id,
        "query": query,
        "messages": [HumanMessage(content=f"Starting drug discovery for: {query}")],
        "target_info": None,
        "known_compounds": None,
        "structure": None,
        "structure_source": None,
        "pockets": None,
        "selected_pocket": None,
        "generated_molecules": None,
        "filtered_molecules": None,
        "screening_results": None,
        "analysis": None,
        "final_report": None,
        "current_step": "initialization",
        "error": None,
        "should_continue": True
    }
    
    # Create and run the graph
    app = create_drug_discovery_graph()
    
    try:
        # Run the workflow
        final_state = app.invoke(initial_state)
        
        # Generate summary
        if final_state.get("target_info"):
            target = final_state["target_info"]
            summary = f"Pipeline completed for {target.get('protein_name')} ({target.get('uniprot_id')})"
            
            if final_state.get("known_compounds"):
                summary += f"\n- Found {len(final_state['known_compounds'])} known compounds"
            
            if final_state.get("structure"):
                structure = final_state["structure"]
                summary += f"\n- Structure: {structure.get('source')} source"
                if structure.get("pdb_id"):
                    summary += f" (PDB: {structure['pdb_id']})"
            
            if final_state.get("selected_pocket"):
                pocket = final_state["selected_pocket"]
                summary += f"\n- Selected pocket: {pocket['pocket_id']} (druggability: {pocket.get('druggability_score', 0):.2f})"
            
            if final_state.get("filtered_molecules"):
                summary += f"\n- Generated {len(final_state.get('generated_molecules', []))} molecules"
                summary += f"\n- Filtered to {len(final_state['filtered_molecules'])} drug-like molecules"
            
            logger.info(summary)
        
        return final_state
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        initial_state["error"] = str(e)
        initial_state["should_continue"] = False
        return initial_state


if __name__ == "__main__":
    # Test the workflow
    test_query = "Find inhibitors for KRAS G12C mutation"
    result = run_drug_discovery_pipeline(test_query)
    
    if result.get("target_info"):
        print(f"✓ Target: {result['target_info'].get('protein_name')}")
        print(f"✓ UniProt: {result['target_info'].get('uniprot_id')}")
        
        if result.get("known_compounds"):
            print(f"✓ Known compounds: {len(result['known_compounds'])}")
        
        if result.get("structure"):
            print(f"✓ Structure: {result['structure'].get('source')}")
            if result['structure'].get('structure_path'):
                print(f"  Path: {result['structure']['structure_path']}")
        
        if result.get("pockets"):
            print(f"✓ Pockets detected: {len(result['pockets'])}")
        
        if result.get("selected_pocket"):
            pocket = result['selected_pocket']
            print(f"✓ Selected pocket: {pocket['pocket_id']}")
            print(f"  Druggability: {pocket.get('druggability_score', 0):.2f}")
            print(f"  Volume: {pocket.get('volume', 0):.1f} Ų")
        
        if result.get("generated_molecules"):
            print(f"✓ Generated molecules: {len(result['generated_molecules'])}")
        
        if result.get("filtered_molecules"):
            print(f"✓ Filtered molecules: {len(result['filtered_molecules'])}")
