import sys
import uuid
import json
import asyncio
from pathlib import Path
from pydantic import BaseModel, Field
from typing import Any, Dict, Optional, List, Literal, AsyncGenerator
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
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

# --- (Pydantic model) ---
class UserIntent(BaseModel):
    """Classifies the user's intent."""
    intent: Literal["chitchat", "drug_discovery", "about_system"] = Field(...)

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

async def get_user_intent(llm: Any, query: str, history: list) -> str:
    """
    Uses a fast LLM to classify the user's query intent.
    """
    llm_with_tool = llm.with_structured_output(UserIntent)
    
    history_messages: List[BaseMessage] = []
    for role, content in history[-4:]:
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))

    prompt: List[BaseMessage] = [
        SystemMessage(content="""Classify the user's intent based on their latest query and the chat history.
    - 'drug_discovery' is for specific scientific tasks like "Find KRAS G12C".
    - 'chitchat' is for simple greetings or conversation like "Hello".
    - 'about_system' is when the user asks what you are, for help, or "What did you just say?".
    """),
    ]
    prompt.extend(history_messages)
    prompt.append(HumanMessage(content=f"Classify this query: \"{query}\""))

    try:
        result = await llm_with_tool.ainvoke(prompt) # Use .ainvoke()
        return result.intent
    except Exception:
        return "drug_discovery"

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
    """
    history_messages: List[BaseMessage] = []
    for role, content in history[-4:]:
        if role == "user":
            history_messages.append(HumanMessage(content=content))
        else:
            history_messages.append(AIMessage(content=content))
    
    prompt: List[BaseMessage] = [
        SystemMessage(content="You are F.A.D.E... A user is making simple chitchat..."),
    ]
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

# --- 2. Setup Logging and Config ---
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

# --- 4. Define Request/Response Models ---
class ChatRequest(BaseModel):
    query: str
    session_id: Optional[str] = None 

class ChatResponse(BaseModel):
    session_id: str
    intent: str
    response_text: Optional[str] = None
    pipeline_state: Optional[Dict[str, Any]] = None

# --- 5. Async Generator for Streaming Pipeline Events ---
# --- Pass app_state to access LLM and Graph App ---
async def stream_pipeline_events(query: str, thread_config: RunnableConfig) -> AsyncGenerator[str, None]:
    """
    Calls the LangGraph app's .astream_events() method and yields 
    SSE-formatted data for messages in real-time.
    """
    
    # Get objects from the shared app_state
    langgraph_app = app_state["langgraph_app"]
    router_llm = app_state["router_llm"]
    
    # Create initial state
    run_id = str(uuid.uuid4())[:8]
    initial_state: DrugDiscoveryState = {
        "run_id": run_id,
        "job_id": None,
        "timestamp": datetime.now(),
        "query": query,
        "messages": [HumanMessage(content=f"Starting drug discovery for: {query}")],
        "current_step": "initialization",
        "should_continue": True,
        "user_id": None
    }

    # Track messages we've already sent
    sent_messages = set()
    accumulated_state = {}
    all_messages = []
    
    # Stream the initial message immediately
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
            # Process different event types
            event_kind = event.get("event", "")
            
            if event_kind == "on_chain_stream":
                # This is a state update from a node
                chunk = event.get("data", {}).get("chunk", {})
                if not chunk:
                    continue
                    
                # Debug logging
                logger.debug(f"[STREAM DEBUG] Chunk keys: {list(chunk.keys())}")
                
                # Process each node's output
                for node_name, node_output in chunk.items():
                    if not isinstance(node_output, dict):
                        continue
                    
                    # Store node output in accumulated state
                    accumulated_state[node_name] = node_output
                    
                    # Check for messages in the node output
                    if "messages" in node_output:
                        messages = node_output["messages"]
                        if not isinstance(messages, list):
                            messages = [messages]
                        
                        for msg in messages:
                            # Extract content from message
                            content = None
                            if hasattr(msg, 'content'):
                                content = msg.content
                            elif isinstance(msg, str):
                                # Sometimes messages are strings with content='...' format
                                import re
                                match = re.search(r"content='([^']+)'", msg)
                                if match:
                                    content = match.group(1)
                                else:
                                    content = msg
                            elif isinstance(msg, dict) and 'content' in msg:
                                content = msg['content']
                            
                            if content:
                                # Create a unique key to avoid duplicates
                                msg_key = f"{node_name}:{content[:50]}"
                                
                                if msg_key not in sent_messages:
                                    sent_messages.add(msg_key)
                                    all_messages.append(content)
                                    
                                    # Stream the message immediately
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
                # A node is starting
                node_name = event.get("name", "unknown")
                logger.info(f"[NODE START] {node_name}")
                # Could yield a node_start event if desired
                
            elif event_kind == "on_node_end":
                # A node completed
                node_name = event.get("name", "unknown")
                logger.info(f"[NODE END] {node_name}")
                # Could yield a node_end event if desired
            
            # Small delay to ensure SSE events are flushed
            await asyncio.sleep(0.01)

    except Exception as e:
        logger.error(f"Streaming error: {e}")
        error_payload = json.dumps({
            "type": "error",
            "content": str(e)
        })
        yield f"event: error\ndata: {error_payload}\n\n"
    
    # Prepare final state
    final_state = accumulated_state
    final_state["messages"] = all_messages
    
    if final_state:
        thread_id = thread_config.get("configurable", {}).get("thread_id", "unknown_thread")
        logger.info(f"Pipeline stream finished for session {thread_id}")

        # Generate summary
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
        
        summary_payload = json.dumps({
            "type": "summary",
            "content": summary,
            "state": clean_state
        }, default=str)
        yield f"event: summary\ndata: {summary_payload}\n\n"
    
    yield "event: done\ndata: {}\n\n"

# --- 6. The /chat-stream Endpoint ---
@app_fastapi.post("/chat-stream")
async def handle_chat_stream(request: ChatRequest):
    """
    Main endpoint for the F.A.D.E chatbot (streaming).
    It classifies intent and returns a StreamingResponse.
    """
    
    session_id = request.session_id or f"api_session_{uuid.uuid4()}"
    thread_config: RunnableConfig = {"configurable": {"thread_id": session_id}}
    query = request.query
    
    # --- : Access objects from app_state ---
    router_llm = app_state["router_llm"]
    project_context = app_state["project_context"]
    sample_queries = app_state["sample_queries"]
    
    logger.info(f"Classifying intent for stream query: {query}")
    intent = await get_user_intent(router_llm, query, history=[])
    
    async def event_generator():
        yield f"event: metadata\ndata: {json.dumps({'session_id': session_id, 'intent': intent})}\n\n"
        
        if intent == "drug_discovery":
            # --- : Pass app_state to the generator ---
            async for event_data in stream_pipeline_events(query, thread_config):
                yield event_data
        
        elif intent == "chitchat":
            response_text = await get_chitchat_response(router_llm, query, history=[])
            yield f"data: {json.dumps({'type': 'message', 'content': response_text.replace('F.A.D.E: ', '')})}\n\n"
            yield "event: done\ndata: {}\n\n"

        elif intent == "about_system":
            response_text = await get_rag_response(router_llm, query, project_context, sample_queries, history=[])
            yield f"data: {json.dumps({'type': 'message', 'content': response_text.replace('F.A.D.E: ', '')})}\n\n"
            yield "event: done\ndata: {}\n\n"
            
    return StreamingResponse(
        event_generator(), 
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"  # Disable nginx buffering if present
        }
    )

# --- 7. The /chat Endpoint ---
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
    
    logger.info(f"Classi_fying intent for query: {query}")
    
    intent = await get_user_intent(router_llm, query, history=[])
    
    if intent == "drug_discovery":
        logger.info(f"Intent: drug_discovery. Starting pipeline for session: {session_id}")
        
        thread_config: RunnableConfig = {"configurable": {"thread_id": session_id}}
        
        final_state = await run_drug_discovery_pipeline(langgraph_app, query, thread_config)

        response_text = await summarize_with_llm(router_llm, final_state)
        
        # -- Create a serializable copy of the state ---
        serializable_state = final_state.copy()
        if "messages" in serializable_state:
            serializable_state["messages"] = [
                {"type": msg.type, "content": msg.content} 
                for msg in serializable_state.get("messages", [])
            ]
        # --------------------------------------------------
        
        return ChatResponse(
            session_id=session_id,
            intent=intent,
            response_text=response_text,
            pipeline_state=serializable_state # <-- Use the clean copy
        )

    elif intent == "chitchat":
        logger.info("Intent: chitchat. Generating simple response...")
        response_text = await get_chitchat_response(router_llm, query, history=[])
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
            history=[]
        )
        return ChatResponse(
            session_id=session_id,
            intent=intent,
            response_text=response_text.replace("F.A.D.E: ", "")
        )

# --- Use 'app_fastapi' ---
@app_fastapi.get("/")
async def read_root():
    return {"hello": "F.A.D.E API is running. POST to /chat to start."}