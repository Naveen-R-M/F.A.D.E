"""
API routes for unified backend - NO FALLBACKS.
Fails immediately if any component is not available.
"""

from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime
import json
import os

# FastAPI imports
from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import JSONResponse

# Import models and intelligence components
from fade.api.models import (
    JobRequest, ChatResponse, JobResponse, 
    ConsentResponse, ErrorResponse, JobStatus,
    MemoryOperation, MemoryResponse, HealthStatus
)
from fade.intelligence.memory import MemoryStore
from fade.intelligence.rag import RAG
from fade.intelligence.classifier import QueryClassifier
from fade.intelligence.gemini_agent import GeminiAgent
from fade.api.pipeline_bridge import pipeline_bridge

# Create router
router = APIRouter(prefix="/api", tags=["unified"])

# Initialize intelligence components
memory_store = None
rag_system = None
classifier = None
gemini_agent = None


def init_intelligence():
    """
    Initialize all intelligence components.
    NO FALLBACKS - Raises exception if any component fails.
    """
    global memory_store, rag_system, classifier, gemini_agent
    
    print("[INIT] Loading Memory Store...", flush=True)
    memory_store = MemoryStore()
    print("[INIT] Memory Store loaded", flush=True)
    
    print("[INIT] Loading RAG System...", flush=True)
    rag_system = RAG()
    print(f"[INIT] RAG loaded with {len(rag_system.sample_queries)} queries, {len(rag_system.context_chunks)} context chunks", flush=True)
    
    print("[INIT] Loading Query Classifier...", flush=True)
    classifier = QueryClassifier()
    print("[INIT] Query Classifier loaded", flush=True)
    
    # Initialize LLM agent - NO FALLBACK
    print("[INIT] Loading LLM Agent...", flush=True)
    from fade.intelligence.gemini_agent import GeminiAgent
    gemini_agent = GeminiAgent()  # Will raise if not available
    print(f"[INIT] LLM Agent active: Gemini", flush=True)
    
    print("[INIT] All intelligence components initialized successfully!", flush=True)


@router.post("/query", response_model_exclude_unset=True)
async def handle_query(
    request: JobRequest,
    background_tasks: BackgroundTasks
):
    """
    Unified endpoint handling both chat and job queries.
    NO FALLBACKS - Raises exception on any error.
    """
    # Ensure components are initialized - NO FALLBACK
    if not classifier:
        print("[QUERY] Components not initialized, loading now...", flush=True)
        init_intelligence()  # Will raise if fails
    
    print(f"[QUERY] Received: '{request.query}'", flush=True)
    print(f"[QUERY] User: {request.user_id}, Consent: {request.consent_given}", flush=True)
    
    # Classify the query - NO ERROR HANDLING
    classification = classifier.classify(request.query)
    print(f"[CLASSIFY] Type: {classification['type']}, Confidence: {classification['confidence']:.2f}", flush=True)
    print(f"[CLASSIFY] Reasoning: {classification['reasoning']}", flush=True)
    
    if classification["type"] == "CHAT":
        # Handle as educational/informational query
        response = await handle_chat_query(request, classification)
        return response
    
    elif classification["type"] == "JOB":
        # Check for user consent
        if not request.consent_given:
            # Generate consent request
            response = await generate_consent_request(request, classification)
            return response
        else:
            # Submit job to LangGraph pipeline
            response = await submit_job_to_langgraph(request, background_tasks)
            return response
    
    # NO DEFAULT CASE - Raise error
    raise ValueError(f"Unknown classification type: {classification['type']}")


async def handle_chat_query(request: JobRequest, classification: Dict) -> ChatResponse:
    """
    Handle chat queries with LLM responses.
    NO FALLBACKS - Requires working LLM.
    """
    # Get RAG context
    context = rag_system.get_context_for_response(request.query)
    
    # Get recent memory
    recent_memory = memory_store.get_recent(request.user_id, k=5)
    
    # Generate response using LLM - NO FALLBACK
    if not gemini_agent or not gemini_agent.initialized:
        raise RuntimeError("LLM Agent not initialized")
    
    gemini_response = gemini_agent.generate_response(
        query=request.query,
        query_type="CHAT",
        context=context,
        memory=recent_memory
    )
    message = gemini_response["message"]
    
    # Store in memory
    memory_store.add(request.user_id, "user", request.query)
    memory_store.add(request.user_id, "assistant", message)
    
    return ChatResponse(
        message=message,
        reasoning=classification.get("reasoning", ""),
        confidence=classification.get("confidence", 0.0),
        context_used=bool(context)
    )


async def generate_consent_request(request: JobRequest, classification: Dict) -> ConsentResponse:
    """
    Generate consent request for LangGraph job submission.
    NO FALLBACKS - Requires working LLM.
    """
    # Use LLM to generate preview - NO FALLBACK
    if not gemini_agent or not gemini_agent.initialized:
        raise RuntimeError("LLM Agent not initialized")
    
    brief_answer = gemini_agent.generate_consent_preview(request.query)
    
    # Create job preview
    job_preview = {
        "query": request.query,
        "estimated_time": "6-8 hours",
        "computational_resources": "HPC cluster with GPU",
        "orchestration": "LangGraph state management",
        "pipeline_stages": [
            "Target Research (RCSB direct search)",
            "Structure Resolution (PDB/AlphaFold3)",
            "Pocket Detection (fpocket/P2Rank/Kalasanty)",
            "Molecule Generation (DiffSBDD)",
            "Screening (ADMET evaluation)",
            "Analysis (comprehensive results)"
        ]
    }
    
    return ConsentResponse(
        message=brief_answer,
        call_to_action="Shall I submit this to the LangGraph drug discovery pipeline?",
        job_preview=job_preview
    )


async def submit_job_to_langgraph(request: JobRequest, background_tasks: BackgroundTasks) -> JobResponse:
    """
    Submit job to LangGraph pipeline.
    NO FALLBACKS - Pipeline must be available.
    """
    if not pipeline_bridge.pipeline_available:
        raise RuntimeError("LangGraph pipeline not available")
    
    # Generate job ID
    import uuid
    job_id = f"job_{str(uuid.uuid4())[:8]}"
    
    # Submit to pipeline
    submission_result = await pipeline_bridge.submit_job(
        job_id=job_id,
        query=request.query,
        user_id=request.user_id
    )
    
    if not submission_result.get("success"):
        raise RuntimeError(f"Pipeline submission failed: {submission_result.get('message', 'Unknown error')}")
    
    # Store in memory
    memory_store.add(request.user_id, "user", request.query)
    memory_store.add(request.user_id, "assistant", f"Job {job_id} submitted to pipeline")
    
    return JobResponse(
        job_id=job_id,
        status="submitted",
        message=submission_result.get("message", "Job submitted"),
        tracking_url=f"/api/jobs/{job_id}/status",
        estimated_time="6-8 hours"
    )


@router.get("/jobs/{job_id}/status", response_model=JobStatus)
async def get_job_status(job_id: str):
    """
    Get status of a specific job.
    NO FALLBACKS - Job must exist.
    """
    job = pipeline_bridge.get_job_status(job_id)
    
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    
    return JobStatus(**job)


@router.get("/health", response_model=HealthStatus)
async def health_check():
    """
    Health check endpoint.
    NO FALLBACKS - All components must be healthy.
    """
    # Check each component
    if not memory_store:
        raise RuntimeError("Memory store not initialized")
    if not rag_system:
        raise RuntimeError("RAG system not initialized")
    if not classifier:
        raise RuntimeError("Classifier not initialized")
    if not gemini_agent or not gemini_agent.initialized:
        raise RuntimeError("LLM Agent not initialized")
    if not pipeline_bridge.pipeline_available:
        raise RuntimeError("LangGraph pipeline not available")
    
    return HealthStatus(
        status="healthy",
        components={
            "memory": "healthy",
            "rag": "healthy",
            "classifier": "healthy",
            "llm": "healthy",
            "pipeline": "healthy"
        },
        timestamp=datetime.utcnow()
    )


@router.post("/memory", response_model=MemoryResponse)
async def manage_memory(operation: MemoryOperation):
    """
    Manage user memory.
    NO FALLBACKS - Memory operations must succeed.
    """
    if operation.operation == "clear":
        memory_store.clear(operation.user_id)
        return MemoryResponse(
            success=True,
            message=f"Memory cleared for user {operation.user_id}"
        )
    
    elif operation.operation == "export":
        memories = memory_store.get_recent(operation.user_id, k=100)
        return MemoryResponse(
            success=True,
            message="Memory exported",
            data=memories
        )
    
    elif operation.operation == "search":
        if not operation.query:
            raise ValueError("Query required for search operation")
        results = memory_store.search(operation.user_id, operation.query)
        return MemoryResponse(
            success=True,
            message=f"Found {len(results)} memories",
            data=results
        )
    
    # NO DEFAULT CASE
    raise ValueError(f"Unknown operation: {operation.operation}")
