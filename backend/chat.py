import sys
import json
from pathlib import Path

import uuid
from pydantic import BaseModel, Field
from typing import Literal, Any, List

# Import your existing pipeline and config logic
from fade.workflows.drug_discovery import run_drug_discovery_pipeline
from fade.config import config
from fade.utils import setup_logging, get_logger
from fade.tools.llm_client import get_llm_client

# --- Import BaseMessage for correct list typing ---
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage, BaseMessage

# --- 1. Pydantic model for Intent Classification ---
class UserIntent(BaseModel):
    """Classifies the user's intent."""
    intent: Literal["chitchat", "drug_discovery", "about_system"] = Field(
        ...,
        description=(
            "Classify the user's intent. "
            "'drug_discovery' for any query about proteins, genes, PDB IDs, drugs, or finding molecules (e.g., 'Find KRAS G12C'). "
            "'chitchat' for greetings, goodbyes, thanks, or how are you (e.g., 'Hello', 'Thanks'). "
            "'about_system' for when the user asks what you are, what you can do, or for help (e.g., 'What is F.A.D.E?', 'help', 'What can you do?')."
        )
    )

# --- 2. Function to load RAG context ---
def load_context_files(logger):
    """
    Loads project context and sample queries from text files.
    """
    data_dir = Path("data")
    project_context = data_dir / "project_context"
    sample_queries = data_dir / "sample_queries.txt"
    try:
        with open(project_context, "r") as f:
            project_context = f.read()
        logger.info("Loaded project_context.txt for RAG.")
    except FileNotFoundError:
        logger.warning("project_context.txt not found. Using default context.")
        project_context = "F.A.D.E is a Fully Agentic Drug Engine for computational drug discovery."

    try:
        with open(sample_queries, "r") as f:
            sample_queries = f.read()
        logger.info("Loaded sample_queries.txt for RAG.")
    except FileNotFoundError:
        logger.warning("sample_queries.txt not found. Using default samples.")
        sample_queries = "Example: 'Find inhibitors for KRAS G12C' or 'Generate structure for P01116'."

    return project_context, sample_queries

# --- 3. Function to get user intent (FIXED) ---
def get_user_intent(llm: Any, query: str, history: list) -> str:
    """
    Uses a fast LLM to classify the user's query intent.
    """
    llm_with_tool = llm.with_structured_output(UserIntent)
    
    # Convert simple list history to LangChain messages
    history_messages: List[BaseMessage] = [] # Explicitly type history list
    for role, content in history[-4:]: # Use last 4 messages for context
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))

    # --- Explicitly type the prompt list as list[BaseMessage] ---
    prompt: List[BaseMessage] = [
        SystemMessage(content="""Classify the user's intent based on their latest query and the chat history.

    - 'drug_discovery' is for specific scientific tasks like "Find KRAS G12C".
    - 'chitchat' is for simple greetings or conversation like "Hello".
    - 'about_system' is when the user asks what you are, for help, or "What did you just say?".
    """),
    ]
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=f"Classify this query: \"{query}\"")) # This is now valid

    try:
        result = llm_with_tool.invoke(prompt)
        return result.intent
    except Exception as e:
        logger.warning(f"Intent classification failed: {e}. Defaulting to 'drug_discovery'.")
        return "drug_discovery"

# --- 4. Function to handle RAG queries (FIXED) ---
def get_rag_response(llm: Any, query: str, context: str, samples: str, history: list) -> str:
    """
    Generates a helpful, context-aware response for 'about_system' queries.
    """
    history_messages: List[BaseMessage] = []
    for role, content in history[-4:]:
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))

    # --- FIX: Explicitly type the prompt list as list[BaseMessage] ---
    prompt: List[BaseMessage] = [
        SystemMessage(content=f"""You are F.A.D.E, a helpful AI assistant for drug discovery.
    Use the following context to answer the user's question about yourself.

    --- Project Context (What you are) ---
    {context}
    ---
    
    --- Sample Queries (What you can do) ---
    {samples}
    ---
    """),
    ]
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=query)) # This is now valid

    try:
        response = llm.invoke(prompt)
        return f"F.A.D.E: {response.content}"
    except Exception as e:
        logger.error(f"RAG response generation failed: {e}")
        return "F.A.D.E: I am F.A.D.E, an agentic workflow for drug discovery."

# --- 5. Function to handle chitchat (FIXED) ---
def get_chitchat_response(llm: Any, query: str, history: list) -> str:
    """
    Generates a simple, friendly response for non-pipeline queries.
    """
    history_messages: List[BaseMessage] = []
    for role, content in history[-4:]:
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))
    
    # --- Explicitly type the prompt list as list[BaseMessage] ---
    prompt: List[BaseMessage] = [
        SystemMessage(content="""You are F.A.D.E, a helpful AI assistant for drug discovery.
    A user is making simple chitchat.
    Give a very short (1 sentence), friendly, professional response.
    """),
    ]
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=query)) # This is now valid
    
    try:
        response = llm.invoke(prompt)
        return f"F.A.D.E: {response.content}"
    except Exception as e:
        logger.error(f"Chitchat response generation failed: {e}")
        return "F.A.D.E: Hello! How can I help?"

# We can reuse your exact text formatting logic
def print_results(final_state: dict):
    """
    Prints the final state in the user-friendly text format.
    (This code is copied directly from your main.py)
    """
    print("\n" + "="*80)
    print("F.A.D.E DRUG DISCOVERY RESULTS")
    print("="*80)
    
    if final_state.get("error"):
        # Check if this is guidance rather than a hard error
        if final_state.get("error_type") in ["no_structures_found", "uniprot_not_found"]:
            # Both RCSB and UniProt guidance use similar format
            print("\n" + "="*80)
            if final_state.get("error_type") == "uniprot_not_found":
                # UniProt guidance is already well-formatted
                print(final_state['error'])
            else:
                # RCSB guidance
                print("💡 QUERY REFINEMENT NEEDED")
                print("="*80)
                print(f"\n{final_state['error']}")
                
                # Show suggested queries if available
                if final_state.get("suggested_queries"):
                    print("\n📝 Try one of these queries:")
                    for i, suggestion in enumerate(final_state["suggested_queries"], 1):
                        print(f"   {i}. {suggestion}") # Changed to be friendlier for chat
                
                # Show detailed guidance if available
                if final_state.get("guidance"):
                    guidance = final_state["guidance"]
                    if guidance.get("refinement") and guidance["refinement"].get("confidence"):
                        confidence = guidance["refinement"]["confidence"]
                        print(f"\n🎯 Confidence in suggestions: {confidence:.0%}")
        else:
            print(f"\n❌ Pipeline Failed: {final_state['error']}")
            print(f"   Step: {final_state.get('current_step', 'unknown')}")
    
    elif final_state.get("target_info"):
        target = final_state["target_info"]
        print(f"\n✓ Target Identified:")
        print(f"  • Protein: {target.get('protein_name', 'Unknown')}")
        print(f"  • UniProt ID: {target.get('unipit_id', 'N/A')}")
        print(f"  • Gene: {target.get('gene_name', 'N/A')}")
        
        if target.get("mutations"):
            print(f"  • Mutations: {', '.join(target['mutations'])}")
        
        if target.get("sequence_length"):
            print(f"  • Sequence Length: {target['sequence_length']} amino acids")
        

        if target.get("existing_structures"):
            structures = target["existing_structures"]
            print(f"\n✓ Existing Structures: {len(structures)} PDB entries")
            for struct in structures[:3]:
                print(f"  • {struct['pdb_id']}: {struct.get('title', 'N/A')[:60]}...")
                if struct.get("resolution"):
                    print(f"    Resolution: {struct['resolution']} Å")
        
        # Show pocket information
        if final_state.get("pockets"):
            pockets = final_state["pockets"]
            print(f"\n✓ Binding Pockets Detected: {len(pockets)}")
            for pocket in pockets[:3]:
                print(f"  • {pocket['pocket_id']}:")
                print(f"    Druggability: {pocket.get('druggability_score', 0):.2f}")
                print(f"    Volume: {pocket.get('volume', 0):.1f} Å³") # Corrected unit
                if pocket.get("description"):
                    print(f"    Description: {pocket['description']}")
        
        if final_state.get("selected_pocket"):
            pocket = final_state["selected_pocket"]
            print(f"\n✓ Selected Target Pocket: {pocket['pocket_id']}")
            if final_state.get("pocket_selection_rationale"):
                print(f"  Rationale: {final_state['pocket_selection_rationale']}")
    
    else:
        print("\n⚠ No target information found")
        print(f"Current step: {final_state.get('current_step', 'unknown')}")
    
    print("="*80 + "\n")

HELP_MESSAGE = """
F.A.D.E: I'm F.A.D.E, an agentic workflow for drug discovery.
You can ask me to find or generate structures and molecules.

Here are some example queries:
  - "Find inhibitors for KRAS G12C"
  - "Structure of EGFR with L858R mutation"
  - "Find PDB structure 4OBE"
  - "Generate molecules for human olfactory receptor 52N4"
  - "Generate structure for P0C0S5"

What drug discovery target can I help you with?
"""

def chat_main():
    """Main entry point for the interactive chatbot."""
    
    # Set up logging (use 'INFO' and 'no_rich' for a cleaner chat interface)
    setup_logging(level="INFO", use_rich=False)
    global logger
    logger = get_logger("chat")
    
    # Validate configuration
    try:
        config.validate()
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        logger.info("Please check your .env file")
        sys.exit(1)
    
    # --- Load RAG context ---
    project_context, sample_queries = load_context_files(logger)

    # --- Initialize the LLM for routing/chitchat ---
    try:
        router_llm = get_llm_client()
    except Exception as e:
        logger.error(f"Could not initialize LLM for routing: {e}")
        sys.exit(1)    
    
    print("="*80)
    print("F.A.D.E - Agentic Drug Engine (Chat Mode)")
    print("="*80)
    print("Type your query to begin. Type 'quit' or 'exit' to stop.")
    
    # --- 1. Create a unique ID for this entire chat session ---
    session_id = f"fade_chat_{uuid.uuid4()}"
    print(f"(Session ID: {session_id})")

    # This config will be passed to LangGraph on every pipeline run
    thread_config = {"configurable": {"thread_id": session_id}}
    
    # --- 2. Create chat_history for conversational memory ---
    chat_history = []
    
    while True:
        try:
            # 1. Get query from user
            query = input("You: ")
            
            # 2. Check for exit/help commands
            if not query:
                continue
            if query.lower() in ["quit", "exit"]:
                print("F.A.D.E: Goodbye!")
                break
            if query.lower() == "help":
                print(HELP_MESSAGE)
                continue

            # --- 3. Route the query ---
            print("F.A.D.E: Classifying intent...")
            intent = get_user_intent(router_llm, query, chat_history)
            
            if intent == "drug_discovery":
                # If it's a real query, run the full pipeline
                print("F.A.D.E: Starting drug discovery pipeline...")
                
                # --- Removed redundant user_id argument ---
                final_state = run_drug_discovery_pipeline(query, thread_config) # type: ignore
                print_results(final_state)
                
            elif intent == "chitchat":
                response_text = get_chitchat_response(router_llm, query, chat_history)
                print(response_text)
                chat_history.append(("user", query))
                chat_history.append(("ai", response_text.replace("F.A.D.E: ", "")))
                
            elif intent == "about_system":
                response_text = get_rag_response(router_llm, query, project_context, sample_queries, chat_history)
                print(response_text)
                chat_history.append(("user", query))
                chat_history.append(("ai", response_text.replace("F.A.D.E: ", "")))
        
        except (KeyboardInterrupt, EOFError):
            print("\nF.A.D.E: Goodbye!")
            break
        except Exception as e:
            # Catch any unexpected errors from the pipeline
            logger.error(f"An unexpected error occurred: {e}", exc_info=True)
            print(f"\n❌ An unexpected error occurred: {e}")


if __name__ == "__main__":
    chat_main()