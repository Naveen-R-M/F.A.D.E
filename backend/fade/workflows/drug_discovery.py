"""
Main LangGraph workflow for drug discovery pipeline with simplified research.

This module defines the StateGraph that orchestrates all agents.
"""

import uuid
import asyncio
import aiosqlite
from datetime import datetime
from typing import Dict, Any, cast

from langgraph.graph import StateGraph, END
from langchain_core.messages import HumanMessage
from langchain_core.runnables import RunnableConfig

from fade.state.langgraph_state import DrugDiscoveryState
from langgraph.checkpoint.sqlite.aio import AsyncSqliteSaver

# Use the simplified research nodes
from fade.agents.research.langgraph_nodes_simplified import (
    target_research_node
)

from fade.agents.structure.langgraph_nodes import (
    structure_resolver_node,
    structure_wait_node
)
from fade.agents.targeting.langgraph_nodes import (
    pocket_mapper_node,
    pocket_ranker_node,
    pocket_from_ligand_node  # New node for holo structures
)
from fade.agents.invention.langgraph_nodes import (
    molecule_generator_node,
    medicinal_filter_node
)
from fade.config import config
from fade.utils import get_logger

logger = get_logger("workflow.drug_discovery")

def entry_router_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Entry router that decides where to start in the graph based on the query and state.
    For continuations, it routes to the appropriate node instead of always starting at research.
    """
    query = state.get("query", "")
    is_continuation = state.get("is_continuation", False)
    
    logger.info(f"[Entry Router] Query: {query[:100]}...")
    logger.info(f"[Entry Router] Is continuation: {is_continuation}")
    
    if is_continuation:
        # Check what we have in state
        has_target = state.get("target_info") is not None
        has_structure = state.get("structure") is not None
        has_pockets = state.get("pockets") is not None or state.get("selected_pocket") is not None
        
        logger.info(f"[Entry Router] State check - Target: {has_target}, Structure: {has_structure}, Pockets: {has_pockets}")
        
        query_lower = query.lower()
        
        # Route based on query intent and available state
        if "properties" in query_lower or "physicochemical" in query_lower:
            logger.info("[Entry Router] Detected property analysis request")
            # Property analysis request
            if not has_target:
                logger.info("[Entry Router] No target info - cannot do property analysis")
                return {
                    "messages": [HumanMessage(content="I need target information first. Please specify a protein target.")],
                    "should_continue": False
                }
            # We have target, provide helpful guidance
            target_info = state['target_info']
            logger.info(f"[Entry Router] Target info keys: {list(target_info.keys()) if isinstance(target_info, dict) else 'not a dict'}")
            
            # Try different fields to get the best display name
            display_name = (target_info.get('protein_name') or 
                          target_info.get('gene_name') or 
                          target_info.get('name') or 
                          'the target')
            uniprot_id = target_info.get('uniprot_id', 'Unknown')
            seq_length = target_info.get('sequence_length', 'Unknown')
            
            message = f"""I understand you want to analyze the physicochemical properties of {display_name} (UniProt: {uniprot_id}).

Here's what we know about {display_name}:
- Protein: {display_name}
- UniProt ID: {uniprot_id}
- Sequence length: {seq_length} amino acids"""
            
            if target_info.get('existing_structures'):
                message += f"\n- Available PDB structures: {len(target_info['existing_structures'])} structures found"
            
            message += f"""

Note: Direct physicochemical property analysis of the protein target is currently under development.

However, I can help you with:
1. **Structure Analysis**: Load and analyze 3D structures from PDB
2. **Pocket Detection**: Identify and characterize binding pockets  
3. **Molecule Generation**: Generate small molecules with specific physicochemical properties
4. **Lead Optimization**: Design molecules with desired ADMET properties

What would you like to explore? For example:
- "Find binding pockets for {display_name}"
- "Load structure for {display_name}"
- "Generate molecules for {display_name}"
"""
            
            logger.info(f"[Entry Router] Providing property analysis guidance for {display_name}")
            
            return {
                "messages": [HumanMessage(content=message)],
                "should_continue": False,
                "current_step": "guidance_provided"
            }
        
        elif "pocket" in query_lower or "binding site" in query_lower:
            # Pocket analysis request
            if not has_structure:
                logger.info("[Entry Router] Pocket request but no structure - routing to structure resolver")
                return {"routing_decision": "structure_resolver"}
            else:
                logger.info("[Entry Router] Pocket request with structure - routing to pocket mapper")
                return {"routing_decision": "pocket_mapper"}
        
        elif "structure" in query_lower or "pdb" in query_lower:
            # Structure request
            logger.info("[Entry Router] Structure request - routing to structure resolver")
            return {"routing_decision": "structure_resolver"}
        
        elif "molecule" in query_lower or "generate" in query_lower:
            # Molecule generation request
            if not has_pockets:
                if not has_structure:
                    logger.info("[Entry Router] Molecule gen but no structure - routing to structure resolver")
                    return {"routing_decision": "structure_resolver"}
                else:
                    logger.info("[Entry Router] Molecule gen but no pockets - routing to pocket mapper")
                    return {"routing_decision": "pocket_mapper"}
            else:
                logger.info("[Entry Router] Molecule gen with pockets - would route to molecule generator")
                return {
                    "messages": [HumanMessage(content="Molecule generation is coming soon!")],
                    "should_continue": False
                }
        
        else:
            # Generic continuation - provide context-aware guidance
            if has_target:
                target_name = state['target_info'].get('protein_name') or state['target_info'].get('gene_name', 'your target')
                return {
                    "messages": [HumanMessage(content=(
                        f"I can help you continue working with {target_name}:\n"
                        "- Find and analyze 3D structures\n"
                        "- Detect binding pockets\n" 
                        "- Generate drug-like molecules\n\n"
                        "Please specify what you'd like to do."
                    ))],
                    "should_continue": False
                }
            else:
                return {
                    "messages": [HumanMessage(content=(
                        "I can help you with:\n"
                        "- Finding structures for your target\n"
                        "- Detecting binding pockets\n"
                        "- Generating drug-like molecules\n\n"
                        "Please specify what you'd like to do."
                    ))],
                    "should_continue": False
                }
    
    else:
        # New query - start from target research
        logger.info("[Entry Router] New query - routing to target research")
        return {"routing_decision": "target_research"}

def route_from_entry(state: DrugDiscoveryState) -> str:
    """
    Routing function for the entry router.
    """
    # Check if we should stop (for guidance messages)
    if not state.get("should_continue", True):
        logger.info("[Entry Router] should_continue=False, ending")
        return "guidance"
    
    # Check for explicit routing decision
    routing = state.get("routing_decision")
    if routing:
        logger.info(f"[Entry Router] Routing to: {routing}")
        return routing
    
    # Default to research for new queries
    logger.info("[Entry Router] No explicit routing, defaulting to target_research")
    return "target_research"

async def create_drug_discovery_graph():
    """
    Create the LangGraph workflow for drug discovery.
    
    Returns:
        Compiled LangGraph application
    """
    conn = await aiosqlite.connect("fade_memory.sqlite")
    memory = AsyncSqliteSaver(conn=conn)
    # Initialize the graph with our state
    workflow = StateGraph(DrugDiscoveryState)
    
    # Add nodes for each step
    workflow.add_node("target_research", target_research_node)
    # workflow.add_node("validate_target", validate_target_node)  # Not implemented yet
    workflow.add_node("structure_resolver", structure_resolver_node)
    workflow.add_node("structure_wait", structure_wait_node)
    workflow.add_node("pocket_mapper", pocket_mapper_node)
    workflow.add_node("pocket_ranker", pocket_ranker_node)
    workflow.add_node("pocket_from_ligand", pocket_from_ligand_node)  # New node
    workflow.add_node("molecule_generator", molecule_generator_node)
    workflow.add_node("medicinal_filter", medicinal_filter_node)
    
    # TODO: Add these nodes when implemented
    # workflow.add_node("boltz2_screener", boltz2_screener_node)
    # workflow.add_node("results_analyzer", results_analyzer_node)

    # Add a router node that decides where to start based on the query
    workflow.add_node("entry_router", entry_router_node)
    
    # Define the edges (workflow sequence)
    workflow.set_entry_point("entry_router")
    
    # Entry router decides where to go
    workflow.add_conditional_edges(
        "entry_router",
        route_from_entry,
        {
            "target_research": "target_research",
            "structure_resolver": "structure_resolver",
            "pocket_mapper": "pocket_mapper",
            "guidance": END,
            "end": END
        }
    )
    
    # Research -> Structure (skipping validation for now)
    workflow.add_conditional_edges(
        "target_research",
        should_continue_research,
        {
            "continue": "structure_resolver",
            "end": END
        }
    )
    
    # # Validation -> Structure (commented out until validate_target is implemented)
    # workflow.add_conditional_edges(
    #     "validate_target",
    #     should_continue_validation,
    #     {
    #         "continue": "structure_resolver",
    #         "end": END
    #     }
    # )
    
    # Structure Resolver -> Pocket from Ligand (Holo) or Pocket Mapping (Apo)
    workflow.add_conditional_edges(
        "structure_resolver",
        # This function decides Holo vs Apo vs Wait
        route_after_structure, 
        {
            "pocket_from_ligand": "pocket_from_ligand", # Holo path
            "pocket_mapper": "pocket_mapper",         # Apo path
            "wait": "structure_wait",                 # Boltz2 path
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
            # "continue": "molecule_generator",
            "continue": END,
            "end": END
        }
    )
    
    # Pocket from Ligand -> Molecule Generator (Holo path)
    workflow.add_conditional_edges(
        "pocket_from_ligand",
        should_continue_targeting,  # Reuse same logic
        {
            # "continue": "molecule_generator",
            "continue": END,
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

    # Set up the memory. This will create a 'fade_memory.sqlite' file.
    app = workflow.compile(checkpointer=memory) # type: ignore
    
    return app


def should_continue_research(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after research."""
    if state.get("error") or not state.get("should_continue", True):
        if state.get("error"):
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

def should_continue_structure(state: DrugDiscoveryState) -> str:
    """Determine if we should continue after structure wait."""
    if state.get("error") or not state.get("should_continue", True):
        logger.warning(f"Structure prediction failed: {state.get('error')}")
        return "end"
    
    return "continue"

def route_after_structure(state: DrugDiscoveryState) -> str:
    """
    Decides whether to extract pocket from ligand (Holo), 
    find pockets (Apo), or wait (Boltz2).
    """
    if state.get("error"):
        logger.warning(f"Structure resolution failed: {state.get('error')}")
        return "end"
        
    structure = state.get("structure", {})
    
    # If Boltz-2 job is pending, wait for it
    if structure.get("source") == "Boltz2" and structure.get("status") == "pending":
        logger.info("Structure source is Boltz2 and status is pending, moving to 'wait' state.")
        return "wait"
    
    # If we have a Holo structure (with a drug-like ligand)
    if structure.get("has_drug_like_ligand"):
        logger.info("Holo structure found. Proceeding to extract pocket from bound ligand.")
        return "pocket_from_ligand" 
    
    # If we have an Apo structure (or AFDB/Boltz2 structure that's ready)
    else:
        logger.info("Apo or predicted structure found. Proceeding to pocket mapping.")
        return "pocket_mapper"

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


async def run_drug_discovery_pipeline(app: Any, query: str, thread_config: RunnableConfig | None = None) -> Dict[str, Any]:
    """
    Run the complete drug discovery pipeline.
    
    Args:
        query: Natural language drug discovery query
        user_id: Optional user identifier
        
    Returns:
        Final state with results
    """
    logger.info(f"Starting pipeline for query: {query}")
    logger.info("Using simplified workflow: Query → RCSB → Structure/Pockets → Molecules")
    
    # Create initial state
    run_id = str(uuid.uuid4())[:8]
    initial_state: DrugDiscoveryState = {
        "run_id": run_id,
        "job_id": None,  # Will be set by first HPC operation
        "timestamp": datetime.now(),
        "user_id": None,
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
    
    # # Create and run the graph
    # app = create_drug_discovery_graph()
    
    try:
        # Run the workflow
        final_state = await app.ainvoke(initial_state, config=thread_config)
        
        # Generate summary
        if final_state.get("target_info"):
            target = final_state["target_info"]
            summary = f"Pipeline completed for {target.get('protein_name')} ({target.get('uniprot_id') or 'PDB-derived'})"
            
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
        
        return cast(Dict[str, Any], final_state)
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        initial_state["error"] = str(e)
        initial_state["should_continue"] = False
        return cast(Dict[str, Any], initial_state)


async def main():
    # Test the workflow
    test_query = "Find inhibitors for KRAS G12C mutation"
    
    # Await the creation of the app
    app = await create_drug_discovery_graph()
    
    # Create a test thread config
    test_thread_id = f"test_run_{uuid.uuid4()}"
    thread_config: RunnableConfig = {"configurable": {"thread_id": test_thread_id}}

    # Await the async function and pass the app
    result = await run_drug_discovery_pipeline(app, test_query, thread_config)
    
    if result.get("target_info"):
        print(f"✓ Target: {result['target_info'].get('protein_name')}")
        print(f"✓ UniProt: {result['target_info'].get('uniprot_id', 'N/A')}")
        
        if result['target_info'].get('existing_structures'):
            structures = result['target_info']['existing_structures']
            print(f"✓ PDB structures: {len(structures)}")
            for s in structures[:3]:
                print(f"  - {s['pdb_id']}: {s.get('title', 'N/A')[:50]}...")
        
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

if __name__ == "__main__":
    asyncio.run(main())