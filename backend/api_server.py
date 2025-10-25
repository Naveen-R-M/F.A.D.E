"""
Unified API server for F.A.D.E backend.
Combines intelligence layer (from gateway) with drug discovery pipeline.
"""

import sys
import logging
from fastapi import FastAPI, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
from pathlib import Path
import os

# Force print statements to appear immediately
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

# Import our unified routes
from fade.api.routes import router, init_intelligence

# Import existing pipeline (will be connected in Phase 5)
# from fade.workflows.drug_discovery import run_drug_discovery_pipeline

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize services on startup."""
    # Startup
    print("=" * 60, flush=True)
    print("F.A.D.E UNIFIED BACKEND STARTUP", flush=True)
    print("=" * 60, flush=True)
    
    print("[STARTUP] Initializing intelligence components...", flush=True)
    init_intelligence()
    
    print("[STARTUP] Checking component status...", flush=True)
    
    # Check which LLM is being used
    try:
        from fade.intelligence.deepseek_agent import DeepSeekAgent
        print("[STARTUP] Using DeepSeek LLM via Ollama", flush=True)
    except ImportError:
        try:
            from fade.intelligence.gemini_agent import GeminiAgent
            if os.getenv("GEMINI_API_KEY"):
                print("[STARTUP] Using Gemini AI", flush=True)
            else:
                print("[STARTUP] No API key - using fallback mode", flush=True)
        except ImportError:
            print("[STARTUP] No LLM agent available - fallback mode", flush=True)
    
    # Check pipeline status
    try:
        from fade.api.pipeline_bridge import pipeline_bridge
        if pipeline_bridge.pipeline_available:
            print("[STARTUP] LangGraph pipeline connected", flush=True)
        else:
            print("[STARTUP] LangGraph pipeline not available - simulation mode", flush=True)
    except ImportError:
        print("[STARTUP] Pipeline bridge not available", flush=True)
    
    print("[STARTUP] All components loaded successfully!", flush=True)
    print("-" * 60, flush=True)
    
    yield
    
    # Shutdown
    print("[SHUTDOWN] Cleaning up...", flush=True)
    print("=" * 60, flush=True)


# Create FastAPI app
app = FastAPI(
    title="F.A.D.E Unified Backend",
    description="Intelligent drug discovery platform with integrated gateway",
    version="4.0-unified",
    lifespan=lifespan
)

# Middleware to log all requests
@app.middleware("http")
async def log_requests(request, call_next):
    """Log all incoming requests."""
    # Skip logging for health checks to reduce noise
    if request.url.path != "/api/health":
        print(f"[REQUEST] {request.method} {request.url.path}", flush=True)
    response = await call_next(request)
    return response

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:3001", 
        "http://127.0.0.1:3000",
        "http://127.0.0.1:3001"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount the unified API routes
app.include_router(router)

# Root endpoint
@app.get("/")
async def root():
    """Root endpoint with system information."""
    print("[ACCESS] Root endpoint accessed", flush=True)
    return {
        "service": "F.A.D.E Unified Backend",
        "version": "4.0",
        "status": "operational",
        "endpoints": {
            "query": "/api/query",
            "health": "/api/health",
            "jobs": "/api/jobs/{job_id}/status",
            "memory": "/api/memory/{user_id}/clear",
            "docs": "/docs"
        }
    }

# Legacy endpoints (for backward compatibility)
@app.post("/jobs")
async def legacy_jobs_endpoint(request: dict, background_tasks: BackgroundTasks):
    """Legacy endpoint that redirects to new unified API."""
    print(f"[LEGACY] Job request received: {request.get('query', 'unknown')}", flush=True)
    
    from fade.api.models import JobRequest
    
    # Convert legacy request to new format
    job_request = JobRequest(
        query=request.get("query", ""),
        user_id=request.get("user_id", "default"),
        consent_given=request.get("consent_given", False)
    )
    
    # Use the new unified handler
    from fade.api.routes import handle_query
    return await handle_query(job_request, background_tasks)


if __name__ == "__main__":
    import uvicorn
    
    port = int(os.getenv("API_PORT", "8000"))
    host = os.getenv("API_HOST", "0.0.0.0")
    
    print(f"Starting F.A.D.E Unified Backend on {host}:{port}", flush=True)
    print(f"Documentation available at http://localhost:{port}/docs", flush=True)
    
    uvicorn.run(
        app,
        host=host,
        port=port,
        reload=True,
        log_level="info"
    )
