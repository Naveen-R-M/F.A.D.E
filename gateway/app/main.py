from pathlib import Path
from fastapi import FastAPI, BackgroundTasks, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from .core.models import JobRequest
from .core.gateway import Gateway
from .agents.open_router_agent import OpenRouterAgent
from .agents.gemini_agent import GeminiToolAgent
from .config import CORS_ORIGINS

app = FastAPI(
    title="F.A.D.E Gateway",
    description="OpenRouter (OpenAI-compatible) tool orchestration + RAG + Memory",
    version="3.1-openrouter"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

base_dir = Path(__file__).resolve().parent.parent
gateway = Gateway(base_dir)
# agent = OpenRouterAgent(gateway)
agent = GeminiToolAgent(gateway)

@app.get("/")
def root():
    return {
        "message": "F.A.D.E Gateway via OpenRouter",
        "status": "healthy",
        "features": [
            "OpenAI-compatible Chat Completions",
            "Function/tool calling with local dispatch",
            "Per-user memory + stdlib RAG",
            "Simulated job queue"
        ]
    }

@app.post("/jobs")
async def handle_request(job_request: JobRequest, background_tasks: BackgroundTasks):
    print(f"\n📥 Request (user {job_request.user_id}): {job_request.query}")

    result = agent.run(job_request.user_id, job_request.query, job_request.id)

    if result["type"] == "job":
        return {
            "type": "job",
            "job_id": result["job_id"],
            "status": result["status"],
            "message": result["message"]
        }
    else:
        return {
            "type": "chat",
            "message": result["message"],
            "reasoning": "OpenRouter tool orchestration",
            "confidence": 0.8
        }

@app.post("/memory/clear/{user_id}")
def clear_user_memory(user_id: str):
    if not user_id:
        raise HTTPException(status_code=400, detail="user_id required")
    gateway.mem_clear(user_id)
    return {"status": "ok", "message": f"Cleared memory for user_id={user_id}"}

@app.get("/health")
def health():
    return {"status": "healthy"}
