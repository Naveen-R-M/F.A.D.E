import sys
import uuid
import json
import asyncio
from pathlib import Path
from pydantic import BaseModel, Field
from typing import Any, Dict, Optional, List, Literal, AsyncGenerator, Union
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from langchain_core.runnables import RunnableConfig
from langchain_core.messages import HumanMessage, AIMessage, SystemMessage, BaseMessage
from datetime import datetime
from contextlib import asynccontextmanager

# --- 1. Import all your F.A.D.E. code ---
from fade.config import config
from fade.utils import setup_logging, get_logger
from fade.tools.llm_client import get_llm_client

# --- : Import the functions, NOT the compiled app ---
from fade.workflows.drug_discovery import (
    create_drug_discovery_graph, # Import the async function that creates the app
    run_drug_discovery_pipeline, 
    DrugDiscoveryState
)

# --- Setup Logging and Config ---
setup_logging(level="INFO", use_rich=False)
logger = get_logger("api")

try:
    config.validate()
    logger.info("Configuration validated.")
except ValueError as e:
    logger.error(f"Configuration error: {e}")
    sys.exit(1)

# This dictionary will hold our globally initialized objects
app_state = {}
# Use a more robust session cache that logs operations
class SessionCache:
    def __init__(self):
        self._cache: Dict[str, List[tuple[str, str]]] = {}
    
    def get(self, session_id: str, default=None):
        result = self._cache.get(session_id, default if default is not None else [])
        logger.info(f"[SessionCache] Retrieved session {session_id}: {len(result)} entries")
        return result
    
    def set(self, session_id: str, history: List[tuple[str, str]]):
        self._cache[session_id] = history
        logger.info(f"[SessionCache] Saved session {session_id}: {len(history)} entries")
        # Log last entry for debugging
        if history:
            role, msg = history[-1]
            logger.info(f"[SessionCache] Last entry: {role}: {msg[:100]}...")
    
    def __setitem__(self, key: str, value: List[tuple[str, str]]):
        self.set(key, value)
    
    def __getitem__(self, key: str):
        return self.get(key, [])

session_cache = SessionCache()

@asynccontextmanager
async def lifespan(app: FastAPI):
    """
    On API startup, load models and create the compiled graph.
    """
    logger.info("Application startup...")
    
    # Load RAG context
    app_state["project_context"], app_state["sample_queries"] = load_context_files(logger)
    
    # Load Router LLM
    try:
        app_state["router_llm"] = get_llm_client()
        logger.info("Router LLM initialized.")
    except Exception as e:
        logger.error(f"Could not initialize router LLM: {e}")
        sys.exit(1)
        
    # Await the creation of the LangGraph app
    try:
        app_state["langgraph_app"] = await create_drug_discovery_graph()
        logger.info("Compiled LangGraph app initialized with AsyncSqliteSaver.")
    except Exception as e:
        logger.error(f"Could not compile LangGraph app: {e}")
        sys.exit(1)
    
    # --- This is where the application runs ---
    yield
    # ----------------------------------------
    
    # Code here would run on shutdown
    logger.info("Application shutdown...")

# --- 3. Load Models & Context (at startup, not in the endpoint) ---
app_fastapi = FastAPI(
    title="F.A.D.E - Agentic Drug Engine",
    description="API for the F.A.D.E drug discovery pipeline.",
    lifespan=lifespan
)

app_fastapi.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins (e.g., localhost:3000)
    allow_credentials=True,
    allow_methods=["*"],  # Allows all methods (GET, POST, OPTIONS, etc.)
    allow_headers=["*"],  # Allows all headers
)

# --- (Pydantic model) ---
class UserIntent(BaseModel):
    """Classifies the user's intent."""
    intent: Literal["chitchat", "drug_discovery", "about_system", "continuation"] = Field(...)
    context: Optional[str] = Field(default=None, description="Additional context for the intent")

# --- (Helper functions - Now defined as ASYNC) ---
def load_context_files(logger):
    """
    Loads project context and sample queries from text files.
    """
    data_dir = Path("data")
    project_context = data_dir / "project_context.txt" # Added .txt
    sample_queries = data_dir / "sample_queries.txt"
    try:
        with open(project_context, "r") as f:
            project_context_text = f.read()
        logger.info("Loaded project_context.txt for RAG.")
    except FileNotFoundError:
        logger.warning("project_context.txt not found. Using default context.")
        project_context_text = "F.A.D.E is a Fully Agentic Drug Engine for computational drug discovery."

    try:
        with open(sample_queries, "r") as f:
            sample_queries_text = f.read()
        logger.info("Loaded sample_queries.txt for RAG.")
    except FileNotFoundError:
        logger.warning("sample_queries.txt not found. Using default samples.")
        sample_queries_text = "Example: 'Find inhibitors for KRAS G12C' or 'Generate structure for P01116'."

    return project_context_text, sample_queries_text

async def get_user_intent(llm: Any, query: str, history: list) -> tuple[str, Optional[str]]:
    """
    Uses a fast LLM to classify the user's query intent.
    Returns (intent, context) tuple.
    """
    llm_with_tool = llm.with_structured_output(UserIntent)
    
    # Check if this might be a continuation of a previous conversation
    is_potential_continuation = False
    last_ai_message = None
    
    logger.info(f"[Intent] Analyzing query: '{query}'")
    logger.info(f"[Intent] History has {len(history)} entries")
    
    # Check if the query itself suggests continuation (uses pronouns without context)
    query_lower = query.lower()
    query_words = query_lower.split()
    
    # Strong indicators that this is a continuation
    pronoun_indicators = ["its", "it's", "it", "this", "that", "these", "those", "their", "them"]
    for pronoun in pronoun_indicators:
        if pronoun in query_words:  # Check for whole word match
            logger.info(f"[Intent] Query contains pronoun '{pronoun}' - likely continuation")
            is_potential_continuation = True
            break
    
    # Special case: "analyze its" is ALWAYS a continuation if there's history
    if "analyze its" in query_lower and history:
        logger.info(f"[Intent] 'analyze its' pattern detected - DEFINITELY continuation")
        is_potential_continuation = True
    
    if history:
        # Look for the last AI message
        for role, content in reversed(history):
            if role == "ai":
                last_ai_message = content
                logger.info(f"[Intent] Last AI message (truncated): {content[:200]}...")
                # Check if the last AI message ended with a question or suggestion
                if any(indicator in content.lower() for indicator in [
                    "what would you like", "perhaps", "would you like to", 
                    "next step", "how would you like to proceed", "please let me know",
                    "you can now", "shall we", "would you prefer", "choose",
                    "options", "alternatively", "or would you", "what aspect",
                    "explore next", "next?", "continue", "proceed", "direction"
                ]):
                    is_potential_continuation = True
                    logger.info(f"[Intent] CONTINUATION DETECTED: AI message contains question/suggestion")
                else:
                    logger.info(f"[Intent] No continuation indicators found in AI message")
                break
    else:
        logger.info(f"[Intent] No history available")
    
    history_messages: List[BaseMessage] = []
    for role, content in history[-4:]:
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))

    prompt: List[BaseMessage] = [
        SystemMessage(content="""You are a router. Classify the user's intent based on their latest query and conversation history.

- **'continuation'**: This is when the user is responding to a previous conversation or using pronouns that refer to something already discussed.
  Key indicators:
  - Uses pronouns like "it", "its", "this", "that" without clear antecedents
  - Responds to questions/suggestions from the AI
  - References something from the previous message
  - Short queries that don't make sense without context
  
  Examples after discussing EGFR:
  "analyze its properties", "what about its structure", "generate molecules for it", "tell me more", "continue"
  
  CRITICAL: "analyze its physicochemical properties" is ALWAYS a continuation if there was a prior target discussed!

- **'drug_discovery'**: This is for NEW specific scientific tasks that start a fresh analysis.
  Examples: "Find inhibitors for KRAS G12C", "Structure of P00533", "Generate molecules for 4OBE".
  These are COMPLETE, STANDALONE requests with explicit target names.

- **'about_system'**: Questions ABOUT the F.A.D.E system itself.
  Examples: "What is F.A.D.E?", "What can you do?", "help".

- **'chitchat'**: Simple greetings or casual conversation.
  Examples: "Hello", "Thanks", "How are you?".

DECISION RULES:
1. If the query uses "it/its" and there's a previous conversation about a target, it's ALWAYS 'continuation'.
2. "analyze its physicochemical properties" → ALWAYS 'continuation' if there's any prior conversation
3. Any query starting with "analyze its" or "analyze it" → ALWAYS 'continuation' if there's history
4. Short queries with pronouns → 'continuation'
5. Queries that don't make sense without context → 'continuation'
""" + (f"\n\nLast AI message: '{last_ai_message[:200]}...'" if last_ai_message else "\n\nNo previous conversation.") + 
            ("\n\nQUERY ANALYSIS: This query uses pronouns ('its') without context, strongly suggesting continuation." if is_potential_continuation else "")),
    ]
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=f"Classify this query: \"{query}\""))

    try:
        result = await llm_with_tool.ainvoke(prompt)
        logger.info(f"Intent classified as: {result.intent} with context: {result.context}")
        return result.intent, result.context
    except Exception as e:
        logger.warning(f"Intent classification failed: {e}. Using fallback logic.")
        
        # Fallback logic with aggressive continuation detection
        query_lower = query.lower()
        query_words = query_lower.split()
        
        # CRITICAL: "analyze its" pattern is ALWAYS a continuation with history
        if "analyze its" in query_lower and history:
            logger.info("[Intent] FALLBACK: 'analyze its' → continuation")
            return "continuation", None
        
        # Strong indicators of continuation
        pronoun_found = any(word in query_words for word in ["its", "it", "this", "that", "these", "those"])
        if pronoun_found and history:
            logger.info(f"[Intent] FALLBACK: Pronoun + history → continuation")
            return "continuation", None
        
        # If we already detected it might be a continuation, trust that
        if is_potential_continuation:
            logger.info("[Intent] FALLBACK: Pre-detected as continuation")
            return "continuation", None
        
        # Check for about_system
        if "help" in query_lower or "what is fade" in query_lower or "who are" in query_lower:
            logger.info("[Intent] FALLBACK: Detected as about_system")
            return "about_system", None
        
        # Short messages are usually chitchat
        if len(query.split()) < 3:
            logger.info("[Intent] FALLBACK: Detected as chitchat")
            return "chitchat", None
        
        # Default to drug_discovery for complex queries
        logger.info("[Intent] FALLBACK: Defaulting to drug_discovery")
        return "drug_discovery", None

async def get_rag_response(llm: Any, query: str, context: str, samples: str, history: list) -> str:
    """
    Generates a helpful, context-aware response for 'about_system' queries.
    """
    history_messages: List[BaseMessage] = []
    for role, content in history[-4:]:
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))

    prompt: List[BaseMessage] = [
        SystemMessage(content=f"""You are F.A.D.E, a helpful AI assistant...
    --- Project Context ---
    {context}
    --- Sample Queries ---
    {samples}
    ---
    """),
    ]
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=query))

    try:
        response = await llm.ainvoke(prompt) # Use .ainvoke()
        return f"F.A.D.E: {response.content}"
    except Exception as e:
        logger.error(f"RAG response generation failed: {e}")
        return "F.A.D.E: I am F.A.D.E, an agentic workflow for drug discovery."

async def get_chitchat_response(llm: Any, query: str, history: list) -> str:
    """
    Generates a simple, friendly response for non-pipeline queries.
    Handles follow-up requests like "keep it concise" or "say that in French".
    """
    history_messages: List[BaseMessage] = []
    for role, content in history[-6:]:  # Include more history for better context
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))
    
    prompt: List[BaseMessage] = [
        SystemMessage(content="""You are F.A.D.E, a helpful AI assistant for drug discovery.

IMPORTANT: Pay close attention to the conversation history. 
- If the user asks you to modify a previous response (e.g., "keep it concise", "say that in French"), 
  you should take your previous response and modify it according to their request.
- If the user is asking for a translation, translate your PREVIOUS response, not their request.
- If the user asks for a shorter/longer version, modify your PREVIOUS response accordingly.
- Always maintain context from the conversation history.

For general chitchat, be friendly and helpful."""),
    ]
    # FIX: Add the current query to the prompt for the LLM
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=query)) 
    
    try:
        response = await llm.ainvoke(prompt) # Use .ainvoke()
        return f"F.A.D.E: {response.content}"
    except Exception as e:
        logger.error(f"Chitchat response generation failed: {e}")
        return "F.A.D.E: Hello! How can I help?"

def _create_summary_context(final_state: dict) -> str:
    """
    Creates a simplified, clean JSON string of the pipeline's key results.
    """
    generated_molecules = final_state.get("generated_molecules")
    filtered_molecules = final_state.get("filtered_molecules")
    summary_data = {
        "query": final_state.get("query"),
        "error": final_state.get("error"),
        "current_step": final_state.get("current_step"),
        "target_info": final_state.get("target_info", {}),
        "structure": final_state.get("structure", {}),
        "selected_pocket": final_state.get("selected_pocket"),
        "generated_molecules_count": len(generated_molecules) if generated_molecules else 0,
        "filtered_molecules_count": len(filtered_molecules) if filtered_molecules else 0,
        "guidance": final_state.get("guidance"),
        "suggested_queries": final_state.get("suggested_queries")
    }
    if "sequence" in summary_data.get("target_info", {}):
        del summary_data["target_info"]["sequence"]
    if "existing_structures" in summary_data.get("target_info", {}):
        del summary_data["target_info"]["existing_structures"]
    return json.dumps(summary_data, indent=2, default=str)

async def summarize_with_llm(llm: Any, final_state: dict) -> str:
    """
    Uses an LLM to generate a natural language summary of the pipeline state.
    """
    summary_context = _create_summary_context(final_state)
    query = final_state.get("query", "the user's request")
    prompt = f"""
    You are F.A.D.E, a conversational AI assistant for drug discovery...
    [Your full summarization prompt]
    {summary_context}
    Your response:
    """
    try:
        response = await llm.ainvoke(prompt)
        content = response.content if hasattr(response, 'content') else str(response)
        return content
    except Exception as e:
        logger.error(f"LLM summarization failed: {e}")
        return final_state.get("error", "An unknown error occurred.")

# --- Define Request/Response Models ---
class ChatRequest(BaseModel):
    query: str
    session_id: Optional[str] = None 

class ChatResponse(BaseModel):
    session_id: str
    intent: str
    response_text: Optional[str] = None
    pipeline_state: Optional[Dict[str, Any]] = None

# --- Async Generator for Streaming Pipeline Events ---
async def debug_state(langgraph_app, thread_config):
    """Debug function to inspect what's in the state."""
    try:
        # Method 1: Get state directly
        state = await langgraph_app.aget_state(thread_config)
        if state:
            logger.info(f"[DEBUG] State from aget_state:")
            if hasattr(state, 'values'):
                values = state.values
                logger.info(f"[DEBUG]   Keys: {list(values.keys()) if values else 'None'}")
                if values and 'target_info' in values:
                    target = values['target_info']
                    if target:
                        logger.info(f"[DEBUG]   Target: {target.get('protein_name')} ({target.get('uniprot_id')})")
                    else:
                        logger.info(f"[DEBUG]   Target is None")
        
        # Method 2: Get checkpoint
        checkpoint_tuple = await langgraph_app.checkpointer.aget_tuple(thread_config)
        if checkpoint_tuple:
            logger.info(f"[DEBUG] Checkpoint tuple found")
            if checkpoint_tuple.checkpoint:
                cp = checkpoint_tuple.checkpoint
                logger.info(f"[DEBUG]   Checkpoint keys: {list(cp.keys())}")
                if 'channel_values' in cp:
                    cv = cp['channel_values']
                    logger.info(f"[DEBUG]   Channel values keys: {list(cv.keys()) if cv else 'None'}")
                    if cv and 'target_info' in cv:
                        logger.info(f"[DEBUG]   Has target_info: {cv['target_info'] is not None}")
    except Exception as e:
        logger.error(f"[DEBUG] Error inspecting state: {e}")

async def stream_pipeline_events(app_state: Dict[str, Any], query: str, thread_config: RunnableConfig, is_continuation: bool = False) -> AsyncGenerator[str, None]:
    """
    Calls the LangGraph app's .astream_events() method and yields 
    SSE-formatted data for messages in real-time.
    """
    
    # Get objects from the shared app_state
    langgraph_app = app_state["langgraph_app"]
    router_llm = app_state["router_llm"]
    
    # Check if we should load existing state for continuation
    existing_state = None
    if is_continuation:
        logger.info(f"[State] Loading state for continuation...")
        
        # Debug what's in the state before loading
        await debug_state(langgraph_app, thread_config)
        
        try:
            # Try to get the latest state from the checkpointer
            logger.info(f"[State] Attempting to load state for continuation, thread: {thread_config}")
            
            # Get the checkpoint
            checkpoint_tuple = await langgraph_app.checkpointer.aget_tuple(thread_config)
            
            if checkpoint_tuple and checkpoint_tuple.checkpoint:
                checkpoint = checkpoint_tuple.checkpoint
                logger.info(f"[State] Found checkpoint with keys: {list(checkpoint.keys())}")
                
                # Try different ways to get the state
                if "channel_values" in checkpoint:
                    existing_state = checkpoint["channel_values"]
                    logger.info(f"[State] Loaded state from channel_values")
                elif "values" in checkpoint:
                    existing_state = checkpoint["values"][0] if checkpoint["values"] else {}
                    logger.info(f"[State] Loaded state from values")
                else:
                    existing_state = checkpoint
                    logger.info(f"[State] Using checkpoint directly as state")
                
                # Log what we found
                if existing_state:
                    logger.info(f"[State] Loaded existing state with keys: {list(existing_state.keys()) if isinstance(existing_state, dict) else 'not a dict'}")
                    if isinstance(existing_state, dict):
                        target_info = existing_state.get('target_info')
                        if target_info:
                            logger.info(f"[State] Found target_info: {target_info.get('protein_name')} ({target_info.get('uniprot_id')})")
                        else:
                            logger.warning(f"[State] No target_info in existing state!")
                            # Try to get the last known good state
                            logger.info(f"[State] Attempting to get state from graph...")
                            graph_state = await langgraph_app.aget_state(thread_config)
                            if graph_state and hasattr(graph_state, 'values'):
                                existing_state = graph_state.values
                                logger.info(f"[State] Got state from graph with keys: {list(existing_state.keys())}")
                                target_info = existing_state.get('target_info')
                                if target_info:
                                    logger.info(f"[State] Found target_info from graph: {target_info.get('protein_name')}")
            else:
                logger.warning(f"[State] No checkpoint found for thread")
                # Try alternative method
                graph_state = await langgraph_app.aget_state(thread_config)
                if graph_state and hasattr(graph_state, 'values'):
                    existing_state = graph_state.values
                    logger.info(f"[State] Got state from aget_state with keys: {list(existing_state.keys())}")
        except Exception as e:
            logger.error(f"[State] Failed to load existing state: {e}", exc_info=True)
    
    # Create initial state
    run_id = str(uuid.uuid4())[:8]
    
    if existing_state and is_continuation:
        # For continuation, preserve existing data
        initial_state = existing_state.copy()
        initial_state["query"] = query
        initial_state["messages"] = initial_state.get("messages", []) + [
            HumanMessage(content=f"Continuation: {query}")
        ]
        initial_state["is_continuation"] = True
        initial_state["should_continue"] = True
        initial_state["error"] = None
        logger.info(f"Using existing state with target: {initial_state.get('target_info', {}).get('protein_name')}")
    else:
        # New query - start fresh
        initial_state: DrugDiscoveryState = {
            "run_id": run_id,
            "job_id": None,
            "timestamp": datetime.now(),
            "query": query,
            "messages": [HumanMessage(content=f"Starting drug discovery for: {query}")],
            "current_step": "initialization",
            "should_continue": True,
            "user_id": None,
            "is_continuation": False
        }

    # Track messages we've already sent
    sent_messages = set()
    # For continuations, start with the existing state
    accumulated_state = initial_state.copy() if is_continuation else {}
    all_messages = []
    
    # Stream the initial message immediately
    if is_continuation:
        initial_msg = f"Continuing analysis: {query}"
    else:
        initial_msg = f"Starting drug discovery for: {query}"
    
    data_payload = json.dumps({
        "type": "message",
        "content": initial_msg,
        "step": "initialization"
    })
    yield f"event: message\ndata: {data_payload}\n\n"
    logger.info(f"[STREAMED] Initial: {initial_msg}")
    all_messages.append(initial_msg)
    sent_messages.add(f"initial:{initial_msg[:50]}")
    
    try:
        # Use astream_events for real-time streaming
        async for event in langgraph_app.astream_events(initial_state, config=thread_config, version="v1"):
            event_kind = event.get("event", "")
            
            if event_kind == "on_chain_stream":
                chunk = event.get("data", {}).get("chunk", {})
                if not chunk:
                    continue
                
                logger.debug(f"[STREAM DEBUG] Chunk keys: {list(chunk.keys())}")
                
                for node_name, node_output in chunk.items():
                    if not isinstance(node_output, dict):
                        continue
                    
                    # Log what we're accumulating - but only for non-continuation runs  
                    if 'target_info' in node_output and not is_continuation:
                        logger.info(f"[STREAM DEBUG] Node {node_name} provided target_info")
                        target = node_output['target_info']
                        if target:
                            logger.info(f"[STREAM DEBUG]   Target: {target.get('protein_name')} ({target.get('uniprot_id')})")
                    
                    # For continuations, check if state is preserved
                    if is_continuation and node_name == "entry_router" and accumulated_state.get('target_info'):
                        logger.info(f"[STREAM DEBUG] Continuation has preserved target_info: {accumulated_state['target_info'].get('protein_name')}")
                    
                    # Store node output in accumulated state
                    accumulated_state.update(node_output)
                    
                    # Check for messages in the node output
                    if "messages" in node_output:
                        messages = node_output["messages"]
                        if not isinstance(messages, list):
                            messages = [messages]
                        
                        for msg in messages:
                            content = None
                            if hasattr(msg, 'content'):
                                content = msg.content
                            elif isinstance(msg, str):
                                import re
                                match = re.search(r"content='([^']+)'", msg)
                                if match:
                                    content = match.group(1)
                                else:
                                    content = msg
                            elif isinstance(msg, dict) and 'content' in msg:
                                content = msg['content']
                            
                            if content:
                                msg_key = f"{node_name}:{content[:50]}"
                                if msg_key not in sent_messages:
                                    sent_messages.add(msg_key)
                                    all_messages.append(content)
                                    
                                    data_payload = json.dumps({
                                        "type": "message",
                                        "content": content,
                                        "step": node_output.get("current_step", node_name)
                                    })
                                    yield f"event: message\ndata: {data_payload}\n\n"
                                    logger.info(f"[STREAMED] {node_name}: {content[:100]}...")
                    
                    # Also check for current_step to send status updates
                    if "current_step" in node_output:
                        status_payload = json.dumps({
                            "type": "status",
                            "step": node_output["current_step"],
                            "node": node_name,
                            "error": node_output.get("error")
                        })
                        yield f"event: status\ndata: {status_payload}\n\n"
                        logger.info(f"[STATUS] {node_name}: {node_output['current_step']}")
                    
                    # Check for errors
                    if "error" in node_output and node_output["error"]:
                        error_payload = json.dumps({
                            "type": "error",
                            "node": node_name,
                            "error": node_output["error"],
                            "step": node_output.get("current_step", "unknown")
                        })
                        yield f"event: error\ndata: {error_payload}\n\n"
                        logger.info(f"[ERROR] {node_name}: {node_output['error']}")
            
            elif event_kind == "on_node_start":
                node_name = event.get("name", "unknown")
                logger.info(f"[NODE START] {node_name}")
            
            elif event_kind == "on_node_end":
                node_name = event.get("name", "unknown")
                logger.info(f"[NODE END] {node_name}")
            
            await asyncio.sleep(0.01)

    except Exception as e:
        logger.error(f"Streaming error: {e}")
        error_payload = json.dumps({"type": "error", "content": str(e)})
        yield f"event: error\ndata: {error_payload}\n\n"
        accumulated_state["error"] = str(e) # Add error to final state
    
    # Prepare final state
    final_state = accumulated_state
    final_state["messages"] = all_messages
    
    # Debug what's being saved after completion
    if not is_continuation:  # Only for initial queries
        logger.info(f"[DEBUG] After initial query completion:")
        logger.info(f"[DEBUG] Final state keys: {list(final_state.keys())}")
        if 'target_info' in final_state:
            target = final_state['target_info']
            if target:
                logger.info(f"[DEBUG] Target saved: {target.get('protein_name')} ({target.get('uniprot_id')})")
            else:
                logger.info(f"[DEBUG] Target is None in final state")
        else:
            logger.info(f"[DEBUG] No target_info in final state!")
        
        # Check what actually got saved to checkpointer
        await debug_state(langgraph_app, thread_config)
    
    if final_state:
        thread_id = thread_config.get("configurable", {}).get("thread_id", "unknown_thread")
        logger.info(f"Pipeline stream finished for session {thread_id}")

        try:
            summary = await summarize_with_llm(router_llm, final_state)
        except Exception as e:
            logger.error(f"Summary generation failed: {e}")
            summary = "Pipeline processing completed."
        
        # Create serializable version of messages
        serializable_messages = []
        for msg in all_messages:
            if isinstance(msg, str):
                serializable_messages.append({"type": "message", "content": msg})
            elif isinstance(msg, dict):
                serializable_messages.append(msg)
            else:
                serializable_messages.append({"type": "message", "content": str(msg)})
        
        # Create clean state for serialization
        clean_state = {}
        for key, value in final_state.items():
            if key == "messages":
                clean_state[key] = serializable_messages
            else:
                clean_state[key] = value
        
        # CRITICAL: Manually save the state to the checkpointer if not a continuation
        if not is_continuation:
            try:
                # Determine what needs to be saved
                state_to_save = {}
                
                if final_state.get('target_info'):
                    state_to_save['target_info'] = final_state['target_info']
                    logger.info(f"[SAVE] Will save target_info: {final_state['target_info'].get('protein_name')} ({final_state['target_info'].get('uniprot_id')})")
                
                if final_state.get('structure'):
                    state_to_save['structure'] = final_state['structure']
                    logger.info(f"[SAVE] Will save structure info")
                
                if final_state.get('selected_pocket'):
                    state_to_save['selected_pocket'] = final_state['selected_pocket']
                    logger.info(f"[SAVE] Will save selected pocket")
                
                if state_to_save:
                    logger.info(f"[SAVE] Manually saving state to checkpointer")
                    
                    # Update the graph state
                    await langgraph_app.aupdate_state(thread_config, state_to_save)
                    logger.info(f"[SAVE] State saved successfully")
                    
                    # Verify it was saved
                    await debug_state(langgraph_app, thread_config)
                else:
                    logger.info(f"[SAVE] No critical state to save")
            except Exception as e:
                logger.error(f"[SAVE] Failed to manually save state: {e}", exc_info=True)
        
        summary_payload = json.dumps({
            "type": "summary",
            "content": summary,
            "state": clean_state
        }, default=str)
        yield f"event: summary\ndata: {summary_payload}\n\n"
    
    yield "event: done\ndata: {}\n\n"

# --- The /chat-stream Endpoint ---
@app_fastapi.post("/chat-stream")
async def handle_chat_stream(request: ChatRequest):
    """
    Main endpoint for the F.A.D.E chatbot (streaming).
    It classifies intent and returns a StreamingResponse.
    """
    
    # CRITICAL: Ensure session_id is maintained across requests
    if not request.session_id:
        session_id = f"api_session_{uuid.uuid4()}"
        logger.warning(f"[Stream] No session_id provided, generated new: {session_id}")
    else:
        session_id = request.session_id
        logger.info(f"[Stream] Using provided session_id: {session_id}")
    
    thread_config: RunnableConfig = {"configurable": {"thread_id": session_id}}
    query = request.query
    
    router_llm = app_state["router_llm"]
    project_context = app_state["project_context"]
    sample_queries = app_state["sample_queries"]

    chat_history = session_cache.get(session_id, [])
    
    logger.info(f"[CRITICAL] Session: {session_id}")
    logger.info(f"[CRITICAL] Query: {query}")
    logger.info(f"[CRITICAL] History length: {len(chat_history)}")
    
    # CRITICAL DEBUG: Show what's in the history
    if chat_history:
        logger.info(f"[CRITICAL] Session {session_id} has existing history:")
        for i, (role, msg) in enumerate(chat_history[-2:]):
            logger.info(f"  Entry {i}: {role}: {msg[:150]}...")
    else:
        logger.warning(f"[CRITICAL] Session {session_id} has NO history! This will cause continuation detection to fail!")
    
    # Get intent BEFORE modifying history
    intent, intent_context = await get_user_intent(router_llm, query, chat_history)
    
    # NOW add the user query to history
    chat_history.append(("user", query))
    
    async def event_generator():
        # ALWAYS return the session_id so frontend can reuse it
        metadata = {
            'session_id': session_id,
            'intent': intent,
            'message': f"Session {session_id} - Intent: {intent}"
        }
        yield f"event: metadata\ndata: {json.dumps(metadata)}\n\n"
        logger.info(f"[Stream] Sent metadata with session_id: {session_id}")
        
        # Log intent decision
        logger.info(f"[Stream] Processing as {intent} for session {session_id}")
        
        if intent == "drug_discovery":
            # Collect the summary from the stream to save to history
            summary_content = None
            async for event_data in stream_pipeline_events(app_state, query, thread_config):
                yield event_data
                # Extract summary from the stream for history
                if "event: summary" in event_data:
                    try:
                        lines = event_data.split("\n")
                        for line in lines:
                            if line.startswith("data: "):
                                data = json.loads(line[6:])
                                if data.get("type") == "summary":
                                    summary_content = data.get("content", "")
                                    logger.info(f"[Stream] Captured summary: {summary_content[:100]}...")
                                    break
                    except Exception as e:
                        logger.error(f"[Stream] Failed to extract summary: {e}")
            
            # Save the AI's response to history for context in next query
            if summary_content:
                chat_history.append(("ai", summary_content))
                session_cache[session_id] = chat_history
                logger.info(f"[Stream] Saved AI response to session {session_id}")
            else:
                # If we didn't get a summary, save a generic response
                fallback_msg = "I've completed the analysis. What would you like to explore next?"
                chat_history.append(("ai", fallback_msg))
                session_cache[session_id] = chat_history
                logger.warning(f"[Stream] No summary captured, saved fallback to session {session_id}")
        
        elif intent == "continuation":
            # For continuations, we need to resume the pipeline with the existing state
            # The checkpointer should have the previous state saved
            logger.info(f"Handling continuation query: {query}")
            
            # Collect the summary from the stream to save to history
            summary_content = None
            async for event_data in stream_pipeline_events(app_state, query, thread_config, is_continuation=True):
                yield event_data
                # Extract summary from the stream for history
                if "event: summary" in event_data:
                    try:
                        lines = event_data.split("\n")
                        for line in lines:
                            if line.startswith("data: "):
                                data = json.loads(line[6:])
                                if data.get("type") == "summary":
                                    summary_content = data.get("content", "")
                                    logger.info(f"[Stream Continuation] Captured summary: {summary_content[:100]}...")
                                    break
                    except Exception as e:
                        logger.error(f"[Stream Continuation] Failed to extract summary: {e}")
            
            # Save the AI's response to history
            if summary_content:
                chat_history.append(("ai", summary_content))
                session_cache[session_id] = chat_history
                logger.info(f"[Stream Continuation] Saved AI response to session {session_id}")
            else:
                fallback_msg = "I've processed your request. What would you like to do next?"
                chat_history.append(("ai", fallback_msg))
                session_cache[session_id] = chat_history
                logger.warning(f"[Stream Continuation] No summary captured, saved fallback to session {session_id}")
        
        elif intent == "chitchat":
            response_text = await get_chitchat_response(router_llm, query, chat_history)
            chat_history.append(("ai", response_text.replace("F.A.D.E: ", "")))
            session_cache[session_id] = chat_history # Save back to cache
            
            yield f"event: summary\ndata: {json.dumps({'type': 'summary', 'content': response_text.replace('F.A.D.E: ', '')})}\n\n"
            yield "event: done\ndata: {}\n\n"

        elif intent == "about_system":
            response_text = await get_rag_response(router_llm, query, project_context, sample_queries, chat_history)
            chat_history.append(("ai", response_text.replace("F.A.D.E: ", "")))
            session_cache[session_id] = chat_history # Save back to cache
            
            yield f"event: summary\ndata: {json.dumps({'type': 'summary', 'content': response_text.replace('F.A.D.E: ', '')})}\n\n"
            yield "event: done\ndata: {}\n\n"
            
    return StreamingResponse(
        event_generator(), 
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"
        }
    )

# --- The /chat Endpoint ---
@app_fastapi.post("/chat", response_model=ChatResponse)
async def handle_chat_request(request: ChatRequest):
    """
    Main endpoint for the F.A.D.E chatbot.
    It classifies intent and routes to the correct workflow.
    """
    
    session_id = request.session_id or f"api_session_{uuid.uuid4()}"
    query = request.query
    
    router_llm = app_state["router_llm"]
    langgraph_app = app_state["langgraph_app"]
    project_context = app_state["project_context"]
    sample_queries = app_state["sample_queries"]

    chat_history = session_cache.get(session_id, [])
    
    logger.info(f"[/chat] Session: {session_id}")
    logger.info(f"[/chat] Query: {query}")
    logger.info(f"[/chat] History length: {len(chat_history)}")
    
    # Get intent BEFORE modifying history
    intent, intent_context = await get_user_intent(router_llm, query, chat_history)
    logger.info(f"[/chat] Intent classified as: {intent}")
    
    # NOW add query to history
    chat_history.append(("user", query))
    session_cache[session_id] = chat_history
    
    if intent == "drug_discovery" or intent == "continuation":
        if intent == "drug_discovery":
            logger.info(f"Intent: drug_discovery. Starting NEW pipeline for session: {session_id}")
        else:
            logger.info(f"Intent: continuation. Resuming pipeline for session: {session_id}")
        
        thread_config: RunnableConfig = {"configurable": {"thread_id": session_id}}
        
        # For continuation, load the existing state and pass it to the pipeline
        existing_state = None
        is_continuation = (intent == "continuation")
        
        if is_continuation:
            logger.info(f"Continuation query will use existing thread state: {session_id}")
            try:
                # Try to get the latest state from the checkpointer
                logger.info(f"[State] Attempting to load state for continuation")
                
                # Try using aget_state first
                graph_state = await langgraph_app.aget_state(thread_config)
                if graph_state and hasattr(graph_state, 'values'):
                    existing_state = graph_state.values
                    logger.info(f"[State] Got state from aget_state with keys: {list(existing_state.keys())}")
                    
                    # Check for target_info
                    if existing_state.get('target_info'):
                        target_info = existing_state['target_info']
                        logger.info(f"[State] Found target_info: {target_info.get('protein_name')} ({target_info.get('uniprot_id')})")
                    else:
                        logger.warning(f"[State] No target_info in state!")
                else:
                    # Fallback to checkpoint
                    checkpoint_tuple = await langgraph_app.checkpointer.aget_tuple(thread_config)
                    if checkpoint_tuple and checkpoint_tuple.checkpoint:
                        checkpoint = checkpoint_tuple.checkpoint
                        if "channel_values" in checkpoint:
                            existing_state = checkpoint["channel_values"]
                            logger.info(f"[State] Loaded from checkpoint channel_values")
            except Exception as e:
                logger.error(f"[State] Could not load existing state: {e}", exc_info=True)
        
        # Create appropriate initial state based on whether it's a continuation
        if existing_state and is_continuation:
            # For continuation, preserve existing data
            initial_state = existing_state.copy()
            initial_state["query"] = query
            initial_state["messages"] = initial_state.get("messages", []) + [
                HumanMessage(content=f"Continuation: {query}")
            ]
            initial_state["is_continuation"] = True
            initial_state["should_continue"] = True
            initial_state["error"] = None
            initial_state["routing_decision"] = None  # Clear any previous routing
            logger.info(f"[/chat] Using existing state with target: {initial_state.get('target_info', {}).get('protein_name') if initial_state.get('target_info') else 'None'}")
        else:
            # New query - start fresh
            run_id = str(uuid.uuid4())[:8]
            initial_state: DrugDiscoveryState = {
                "run_id": run_id,
                "job_id": None,
                "timestamp": datetime.now(),
                "query": query,
                "messages": [HumanMessage(content=f"Starting drug discovery for: {query}")],
                "current_step": "initialization",
                "should_continue": True,
                "user_id": None,
                "is_continuation": False
            }
        
        # Run the pipeline with the appropriate initial state
        final_state = await langgraph_app.ainvoke(initial_state, config=thread_config)

        response_text = await summarize_with_llm(router_llm, final_state)

        chat_history.append(("ai", response_text))
        session_cache[session_id] = chat_history
        
        serializable_state = final_state.copy()
        if "messages" in serializable_state and serializable_state.get("messages"):
            serializable_state["messages"] = [
                {"type": msg.type, "content": msg.content} 
                for msg in serializable_state["messages"]
                if hasattr(msg, 'type') and hasattr(msg, 'content')
            ]
        
        return ChatResponse(
            session_id=session_id,
            intent=intent,
            response_text=response_text,
            pipeline_state=serializable_state
        )

    elif intent == "chitchat":
        logger.info("Intent: chitchat. Generating simple response...")
        response_text = await get_chitchat_response(router_llm, query, chat_history)
        chat_history.append(("ai", response_text.replace("F.A.D.E: ", "")))
        session_cache[session_id] = chat_history
        return ChatResponse(
            session_id=session_id,
            intent=intent,
            response_text=response_text.replace("F.A.D.E: ", "")
        )

    elif intent == "about_system":
        logger.info("Intent: about_system. Generating RAG response...")
        response_text = await get_rag_response(
            router_llm, 
            query, 
            project_context, 
            sample_queries, 
            chat_history
        )
        chat_history.append(("ai", response_text.replace("F.A.D.E: ", "")))
        session_cache[session_id] = chat_history
        return ChatResponse(
            session_id=session_id,
            intent=intent,
            response_text=response_text.replace("F.A.D.E: ", "")
        )

# --- Use 'app_fastapi' ---
@app_fastapi.get("/")
async def read_root():
    return {"hello": "F.A.D.E API is running. POST to /chat to start."}

@app_fastapi.get("/debug/sessions")
async def debug_sessions():
    """Debug endpoint to check active sessions."""
    sessions_info = {}
    for session_id in session_cache._cache:
        history = session_cache._cache[session_id]
        sessions_info[session_id] = {
            "entries": len(history),
            "last_user": history[-1][1][:100] if history and history[-1][0] == "user" else None,
            "last_ai": None
        }
        # Find last AI message
        for role, msg in reversed(history):
            if role == "ai":
                sessions_info[session_id]["last_ai"] = msg[:100]
                break
    
    return {
        "active_sessions": len(session_cache._cache),
        "sessions": sessions_info
    }