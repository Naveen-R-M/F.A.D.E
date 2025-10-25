"""
API routes for unified backend with Gemini AI integration.
Connected to LangGraph drug discovery pipeline.
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
from fade.intelligence.deepseek_agent import DeepSeekAgent as GeminiAgent
from fade.api.pipeline_bridge import pipeline_bridge

# Create router
router = APIRouter(prefix="/api", tags=["unified"])

# Initialize intelligence components
memory_store = None
rag_system = None
classifier = None
gemini_agent = None


def init_intelligence():
    """Initialize all intelligence components."""
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
    
    # Initialize LLM agent (DeepSeek or Gemini)
    print("[INIT] Loading LLM Agent...", flush=True)
    try:
        # Try DeepSeek first if available
        from fade.intelligence.deepseek_agent import DeepSeekAgent
        gemini_agent = DeepSeekAgent()
        agent_name = "DeepSeek via Ollama"
    except ImportError:
        # Fall back to Gemini
        from fade.intelligence.gemini_agent import GeminiAgent
        api_key = os.getenv("GEMINI_API_KEY")
        gemini_agent = GeminiAgent(api_key)
        agent_name = "Gemini AI" if api_key else "Fallback mode"
    
    if gemini_agent.initialized:
        print(f"[INIT] LLM Agent active: {agent_name}", flush=True)
    else:
        print(f"[INIT] LLM Agent not initialized - using fallback responses", flush=True)
    
    print("[INIT] All intelligence components initialized successfully!", flush=True)


@router.post("/query", response_model_exclude_unset=True)
async def handle_query(
    request: JobRequest,
    background_tasks: BackgroundTasks
):
    """
    Unified endpoint handling both chat and job queries.
    Uses LangGraph for job orchestration.
    """
    try:
        # Ensure components are initialized
        if not classifier:
            print("[QUERY] Components not initialized, loading now...", flush=True)
            init_intelligence()
        
        print(f"[QUERY] Received: '{request.query}'", flush=True)
        print(f"[QUERY] User: {request.user_id}, Consent: {request.consent_given}", flush=True)
        
        # Classify the query
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
        
    except Exception as e:
        print(f"[ERROR] Query handling failed: {e}")
        return ErrorResponse(
            error="Failed to process query",
            details=str(e)
        )


async def handle_chat_query(request: JobRequest, classification: Dict) -> ChatResponse:
    """Handle chat queries with Gemini-enhanced responses."""
    
    # Get RAG context
    context = rag_system.get_context_for_response(request.query)
    
    # Get recent memory
    recent_memory = memory_store.get_recent(request.user_id, k=5)
    
    # Generate response using Gemini
    if gemini_agent and gemini_agent.initialized:
        gemini_response = gemini_agent.generate_response(
            query=request.query,
            query_type="CHAT",
            context=context,
            memory=recent_memory
        )
        message = gemini_response["message"]
        ai_generated = True
    else:
        # Fallback to context-based response
        if context:
            message = f"According to F.A.D.E documentation: {context[:500]}..."
        else:
            message = "F.A.D.E is a LangGraph-powered drug discovery platform that designs molecules targeting specific proteins."
        ai_generated = False
    
    # Store in memory
    memory_store.add(request.user_id, "user", request.query)
    memory_store.add(request.user_id, "assistant", message)
    
    return ChatResponse(
        message=message,
        reasoning=classification.get("reasoning", "Informational query"),
        confidence=classification.get("confidence", 0.5),
        context_used=bool(context)
    )


async def generate_consent_request(request: JobRequest, classification: Dict) -> ConsentResponse:
    """Generate consent request for LangGraph job submission."""
    
    # Use Gemini to generate preview if available
    if gemini_agent and gemini_agent.initialized:
        brief_answer = gemini_agent.generate_consent_preview(request.query)
    else:
        # Fallback preview
        brief_answer = f"""The LangGraph F.A.D.E pipeline will process: "{request.query}"

• Target Research: Identify protein target
• Structure Resolution: Get 3D structure
• Pocket Detection: Find binding sites
• Molecule Generation: Create candidates
• Screening: Evaluate properties
• Analysis: Compile results

Uses HPC cluster (6-8 hours)."""
    
    # Create job preview
    job_preview = {
        "query": request.query,
        "estimated_time": "6-8 hours",
        "computational_resources": "HPC cluster with GPU",
        "orchestration": "LangGraph state management",
        "pipeline_stages": [
            "Target Research (UniProt integration)",
            "Structure Resolution (AlphaFold3/PDB)",
            "Pocket Detection (fpocket)",
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
    """Submit job to LangGraph drug discovery pipeline."""
    
    # Store query in memory
    memory_store.add(request.user_id, "user", request.query)
    memory_store.add(request.user_id, "assistant", f"Submitted job {request.id} to LangGraph pipeline")
    
    # Add background task to run LangGraph workflow
    background_tasks.add_task(
        run_langgraph_workflow,
        job_id=request.id,
        query=request.query,
        user_id=request.user_id
    )
    
    return JobResponse(
        job_id=request.id,
        status="submitted",
        message=f"Drug discovery job submitted to LangGraph pipeline. Estimated time: 6-8 hours.",
        tracking_url=f"/api/jobs/{request.id}/status",
        estimated_time="6-8 hours"
    )


async def run_langgraph_workflow(job_id: str, query: str, user_id: str):
    """Execute LangGraph drug discovery workflow via pipeline bridge."""
    print(f"[LANGGRAPH] Starting workflow for job {job_id}: {query}")
    
    try:
        # Submit to pipeline bridge
        result = await pipeline_bridge.submit_job(
            job_id=job_id,
            query=query,
            user_id=user_id
        )
        
        if result["success"]:
            print(f"[LANGGRAPH] Job {job_id} submitted successfully")
        else:
            print(f"[ERROR] Failed to submit job {job_id}")
            
    except Exception as e:
        print(f"[ERROR] LangGraph workflow submission failed: {e}")


@router.get("/jobs/{job_id}/status", response_model=JobStatus)
async def get_job_status(job_id: str):
    """Get status of LangGraph job from pipeline bridge."""
    
    job = pipeline_bridge.get_job_status(job_id)
    
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    
    return JobStatus(
        id=job["id"],
        status=job["status"],
        progress=job.get("progress", 0),
        message=job.get("message", ""),
        created_at=datetime.fromisoformat(job["created_at"]),
        updated_at=datetime.fromisoformat(job["updated_at"]),
        completed_at=datetime.fromisoformat(job["completed_at"]) if job.get("completed_at") else None,
        results=job.get("results")
    )


@router.post("/memory/{user_id}/clear", response_model=MemoryResponse)
async def clear_user_memory(user_id: str):
    """Clear conversation memory for user."""
    
    if not memory_store:
        init_intelligence()
    
    memory_store.clear(user_id)
    
    return MemoryResponse(
        status="ok",
        message=f"Cleared memory for user {user_id}"
    )


@router.get("/memory/{user_id}/recent", response_model=MemoryResponse)
async def get_recent_memory(user_id: str, k: int = 10):
    """Get recent conversation history."""
    
    if not memory_store:
        init_intelligence()
    
    recent = memory_store.get_recent(user_id, k)
    
    return MemoryResponse(
        status="ok",
        message=f"Retrieved {len(recent)} recent messages",
        data=recent
    )




@router.get("/jobs", response_model=list)
async def list_jobs(user_id: Optional[str] = None):
    """List all jobs or jobs for specific user."""
    return pipeline_bridge.list_jobs(user_id)

@router.get("/health", response_model=HealthStatus)
async def health_check():
    """Check health status including Gemini availability."""
    
    # Check component status
    intelligence_status = {
        "memory": memory_store is not None,
        "rag": rag_system is not None,
        "classifier": classifier is not None,
        "gemini": gemini_agent is not None and gemini_agent.initialized
    }
    
    # LangGraph pipeline availability (will be true in Phase 5)
    langgraph_available = pipeline_bridge.pipeline_available
    
    return HealthStatus(
        status="healthy" if all(intelligence_status.values()[:3]) else "degraded",
        service="fade-unified-backend",
        version="4.0-gemini",
        intelligence=intelligence_status,
        pipeline=langgraph_available
    )
