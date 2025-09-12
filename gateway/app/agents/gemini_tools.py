# app/agents/gemini_tools.py
from __future__ import annotations
import json
from typing import Dict, Any, Optional
from ..core.gateway import Gateway

# Optional HPC integration
try:
    from ..hpc.runner import submit_fade_job  # def submit_fade_job(job_id:str, prompt:str)->Dict[str,Any]
except Exception:
    submit_fade_job = None  # type: ignore[misc]


def build_tools_and_dispatch(gateway: Gateway, tool_state: Optional[Dict[str, Any]] = None):
    """
    Returns:
      gemini_tools_spec: list[dict] -> google.generativeai "tools" spec (function_declarations)
      dispatch:          dict[name -> callable] -> executes the tool
    """

    if tool_state is None:
        tool_state = {}
    # budgets & caches in shared state
    tool_state.setdefault("budget", {"memory_read": 1, "rag_search": 1, "job_submit": 1, "job_status": 2})
    tool_state.setdefault("calls", {})        # e.g., {"memory_read": 1}
    tool_state.setdefault("rag_cache", {})    # key: (query,k) -> hits

    def _use_budget(name: str) -> bool:
        calls = tool_state["calls"].get(name, 0)
        if calls >= tool_state["budget"].get(name, 0):
            return False
        tool_state["calls"][name] = calls + 1
        return True

    # ---- Tool implementations (server side) ----
    def memory_read(user_id: str, k: int) -> Dict[str, Any]:
        if not _use_budget("memory_read"):
            return {"turns": [], "note": "memory_read_budget_exhausted"}
        if not (1 <= k <= 12):
            raise ValueError("k must be between 1 and 12")
        turns = gateway.mem_recent(user_id, k)
        # keep content short to reduce tokens
        for t in turns:
            if len(t["content"]) > 400:
                t["content"] = t["content"][:400] + " …"
        return {"turns": turns}


    def memory_write(user_id: str, role: str, content: str) -> Dict[str, Any]:
        if role not in ("user", "assistant"):
            raise ValueError("role must be 'user' or 'assistant'")
        if not content:
            raise ValueError("content must be non-empty")
        gateway.mem_add(user_id, role, content)
        return {"ok": True}

    def rag_search(query: str, k: int) -> Dict[str, Any]:
        if not _use_budget("rag_search"):
            return {"hits": [], "note": "rag_search_budget_exhausted"}
        if not (1 <= k <= 10):
            raise ValueError("k must be between 1 and 10")

        key = (query.strip().lower(), int(k))
        cache = tool_state["rag_cache"]
        if key in cache:
            return cache[key]

        hits = gateway.rag_search(query, k)
        # trim & normalize payload
        slim = []
        for h in hits:
            snippet = h.get("snippet", "")
            if len(snippet) > 360:
                snippet = snippet[:360] + " …"
            slim.append({"pid": h.get("pid", ""), "score": float(h.get("score", 0.0)), "snippet": snippet})
        result = {"hits": slim}
        cache[key] = result
        return result


    def job_submit(job_id: str, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Submit a F.A.D.E job to HPC (Nextflow). Expect a natural-language 'prompt'.
        """
        prompt = (
            payload.get("prompt")
            or payload.get("query")
            or payload.get("natural_query")
            or payload.get("task")
        )
        if not prompt:
            return {"ok": False, "error": "payload.prompt (or query/natural_query/task) is required"}

        if submit_fade_job is None:
            return {"ok": False, "error": "HPC runner not configured (submit_fade_job unavailable)"}

        try:
            result = submit_fade_job(job_id, prompt)
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}

        ok = bool(result.get("ok", False))

        if ok and tool_state is not None:
            tool_state["last_job_id"] = result.get("job_id") or job_id

        try:
            gateway.job_upsert(job_id, "submitted" if result.get("ok") else "error",
                            message="Background submission started" if result.get("ok") else "Submission failed",
                            cwd=result.get("cwd", ""))
        except Exception:
            pass

        status = "submitted" if result.get("ok") else "error"
        logdir = result.get("cwd", "")
        
        return {
            "ok": bool(result.get("ok")),
            "job_id": result.get("job_id", job_id),
            "status": status,
            "message": f"Nextflow submission started (background). Logs: {logdir}/submit.log and {logdir}/run.log",
            "cwd": logdir,
            "command": result.get("command", ""),
            "exit_code": result.get("exit_code"),
            "stdout": result.get("stdout", "")[:4000],
            "stderr": result.get("stderr", "")[:4000],
        }

    def job_status(job_id: str) -> Dict[str, Any]:
        job = gateway.job_status(job_id)
        if job:
            try:
                return job.model_dump()   # pydantic v2
            except Exception:
                return job.dict()         # pydantic v1
        return {"error": "not_found"}

    # ---- Gemini tool spec (Function Declarations)
    # NOTE: keep schemas MINIMAL: {type, properties, required}. No 'additionalProperties'.
    gemini_tools_spec = [{
        "function_declarations": [
            {
                "name": "memory_read",
                "description": "Return the last k conversation turns for a user to resolve follow-ups.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "user_id": {"type": "string", "description": "Stable user identifier."},
                        "k": {"type": "integer"}
                    },
                    "required": ["user_id", "k"]
                }
            },
            {
                "name": "memory_write",
                "description": "Append one message to memory (side-effect).",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "user_id": {"type": "string"},
                        "role": {"type": "string"},
                        "content": {"type": "string"}
                    },
                    "required": ["user_id", "role", "content"]
                }
            },
            {
                "name": "rag_search",
                "description": "Retrieve top-k relevant snippets from local corpus. Use only if helpful; otherwise answer from your own knowledge.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "query": {"type": "string"},
                        "k": {"type": "integer"}
                    },
                    "required": ["query", "k"]
                }
            },
            {
                "name": "job_submit",
                "description": "Submit a pipeline job (idempotent if job_id reused).",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "job_id": {"type": "string"},
                        "payload": {"type": "object"}
                    },
                    "required": ["job_id", "payload"]
                }
            },
            {
                "name": "job_status",
                "description": "Fetch job status by job_id.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "job_id": {"type": "string"}
                    },
                    "required": ["job_id"]
                }
            }
        ]
    }]

    dispatch = {
        "memory_read": memory_read,
        "memory_write": memory_write,
        "rag_search": rag_search,
        "job_submit": job_submit,
        "job_status": job_status,
    }

    return gemini_tools_spec, dispatch
