# app/agents/openrouter_tools.py
from typing import Dict, Any, List, Literal
from ..core.gateway import Gateway
from ..hpc.runner import submit_fade_job

def build_tools_and_dispatch(gateway: Gateway):
    """
    Returns:
      tools_spec: list[dict] -> OpenAI Chat Completions "tools" spec (JSON schema)
      dispatch:   dict[name -> callable] -> executes the tool
    """

    # ---- Tool implementations (server side) ----
    def memory_read(user_id: str, k: int) -> Dict[str, Any]:
        if not (1 <= k <= 12):
            raise ValueError("k must be between 1 and 12")
        turns = gateway.mem_recent(user_id, k)
        return {"turns": turns}

    def memory_write(user_id: str, role: str, content: str) -> Dict[str, Any]:
        if role not in ("user", "assistant"):
            raise ValueError("role must be 'user' or 'assistant'")
        if not content:
            raise ValueError("content must be non-empty")
        gateway.mem_add(user_id, role, content)
        return {"ok": True}

    def rag_search(query: str, k: int) -> Dict[str, Any]:
        if not (1 <= k <= 10):
            raise ValueError("k must be between 1 and 10")
        hits = gateway.rag_search(query, k)
        return {"hits": hits}

    def job_submit(job_id: str, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Submit a F.A.D.E Nextflow job by SSH'ing to HPC and running the submission script
        from the scratch space. Payload should include the natural-language query.
        """
        # Accept several common keys; prefer 'prompt'
        prompt = (
            payload.get("prompt")
            or payload.get("query")
            or payload.get("natural_query")
            or payload.get("task")
        )
        if not prompt:
            return {"error": "payload.prompt (or query/natural_query) is required"}

        result = submit_fade_job(job_id, prompt)

        # Update in-memory job store status based on exit code
        if result.get("ok"):
            # Submitted successfully
            status = "submitted"
            message = "Nextflow submission succeeded"
        else:
            status = "error"
            message = f"Submission failed (exit {result.get('exit_code')}). See stderr."

        # Upsert into your job store if available
        try:
            # NOTE: gateway is captured from outer scope when build_tools_and_dispatch() is called
            # If you prefer, expose a setter on Gateway/JobStore to keep this clean.
            pass
        except Exception:
            pass

        # Return a concise summary to the model + full logs for debugging if needed
        return {
            "job_id": result.get("job_id"),
            "status": status,
            "message": message,
            "cwd": result.get("cwd"),
            "command": result.get("command"),
            "exit_code": result.get("exit_code"),
            "stdout": result.get("stdout", "")[:4000],  # trim for token budget
            "stderr": result.get("stderr", "")[:4000],
        }

    def job_status(job_id: str) -> Dict[str, Any]:
        job = gateway.job_status(job_id)
        return (job.model_dump() if hasattr(job, "model_dump")
                else (job.dict() if job else {"error": "not_found"}))

    # ---- OpenAI tool specs (what the model sees) ----
    tools_spec = [
        {
            "type": "function",
            "function": {
                "name": "memory_read",
                "description": "Return the last k conversation turns for a user to resolve follow-ups.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "user_id": {"type": "string", "description": "Stable user identifier."},
                        "k": {"type": "integer", "minimum": 1, "maximum": 12}
                    },
                    "required": ["user_id", "k"],
                    "additionalProperties": False
                }
            }
        },
        {
            "type": "function",
            "function": {
                "name": "memory_write",
                "description": "Append one message to memory (side-effect).",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "user_id": {"type": "string"},
                        "role": {"type": "string", "enum": ["user", "assistant"]},
                        "content": {"type": "string", "minLength": 1}
                    },
                    "required": ["user_id", "role", "content"],
                    "additionalProperties": False
                }
            }
        },
        {
            "type": "function",
            "function": {
                "name": "rag_search",
                "description": "Retrieve top-k relevant snippets from local corpus. Use only if helpful; otherwise answer from your own knowledge.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "query": {"type": "string"},
                        "k": {"type": "integer", "minimum": 1, "maximum": 10}
                    },
                    "required": ["query", "k"],
                    "additionalProperties": False
                }
            }
        },
        {
            "type": "function",
            "function": {
                "name": "job_submit",
                "description": "Submit a pipeline job (idempotent if job_id reused).",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "job_id": {"type": "string"},
                        "payload": {"type": "object"}
                    },
                    "required": ["job_id", "payload"],
                    "additionalProperties": False
                }
            }
        },
        {
            "type": "function",
            "function": {
                "name": "job_status",
                "description": "Fetch job status by job_id.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "job_id": {"type": "string"}
                    },
                    "required": ["job_id"],
                    "additionalProperties": False
                }
            }
        }
    ]

    dispatch = {
        "memory_read": memory_read,
        "memory_write": memory_write,
        "rag_search": rag_search,
        "job_submit": job_submit,
        "job_status": job_status,
    }

    return tools_spec, dispatch
