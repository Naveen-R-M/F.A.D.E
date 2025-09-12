# app/agents/openrouter_agent.py
import json
from typing import Dict, Any, Optional

from openai import OpenAI
from ..config import (
    OPENROUTER_API_KEY, OPENROUTER_BASE_URL, OPENROUTER_MODEL,
    OPENROUTER_REFERRER, OPENROUTER_TITLE
)
from ..core.gateway import Gateway
from .open_router_tools import build_tools_and_dispatch

SYSTEM_PROMPT = """
You are F.A.D.E, an AI drug discovery copilot.
POLICY:
  - Always answer the user's question directly in THIS message (self-contained).
  - Never say you already answered, “as above,” or “I’ve provided an overview.”
  - Use tools only when needed; otherwise answer from your own knowledge.
  - If RAG returns weak/irrelevant passages, ignore them and answer from your knowledge.
  - For workflow requests, propose the job payload and next steps clearly.
  - Keep answers concise, structured, and actionable.
"""

class OpenRouterAgent:
    def __init__(self, gateway: Gateway):
        if not OPENROUTER_API_KEY:
            raise RuntimeError("Missing OPENROUTER_API_KEY.")
        self.client = OpenAI(base_url=OPENROUTER_BASE_URL, api_key=OPENROUTER_API_KEY)
        self.model = OPENROUTER_MODEL
        self.gateway = gateway
        self.tools_spec, self.dispatch = build_tools_and_dispatch(gateway)

        # Optional headers for OpenRouter ranking/analytics
        self.extra_headers = {}
        if OPENROUTER_REFERRER:
            self.extra_headers["HTTP-Referer"] = OPENROUTER_REFERRER
        if OPENROUTER_TITLE:
            self.extra_headers["X-Title"] = OPENROUTER_TITLE

        # FIX 2: declare tools_supported so Pylance is happy
        self.tools_supported: Optional[bool] = None

    def _compose_server_context(self, user_id: str, user_query: str, request_id: str) -> str:
        """
        Build a concise, server-side context block (recent memory + RAG hits).
        Sent as a system message so the model doesn't think it already answered.
        """
        recent = self.gateway.mem_recent(user_id, 6)
        mem_lines = [f"{t['role']}: {t['content']}" for t in recent]
        memory_block = "\n".join(mem_lines) if mem_lines else "(none)"

        hits = self.gateway.rag_search(user_query, 6)
        rag_lines = [f"[{i+1}] ({h['pid']}, score={h.get('score', 0):.2f}) {h['snippet']}" for i, h in enumerate(hits)]
        rag_block = "\n".join(rag_lines) if rag_lines else "(none)"

        top_score = hits[0]["score"] if hits else 0.0
        guidance = "RAG looks informative." if top_score >= 0.40 else "RAG looks weak; prefer your own knowledge."

        return (
            "SERVER_CONTEXT\n"
            f"- Request: {json.dumps({'user_id': user_id, 'request_id': request_id})}\n\n"
            f"Recent conversation:\n{memory_block}\n\n"
            f"RAG candidates:\n{rag_block}\n\n"
            f"Guidance: {guidance}\n"
            "Instruction: Use the above only if helpful; otherwise answer from your own knowledge."
        )

    def _call(self, messages: list, tools: Optional[list] = None):
        kwargs = {
            "model": self.model,
            "messages": messages,
            "tool_choice": "auto" if tools else "none",
            "temperature": 0.3,
            "extra_headers": self.extra_headers or None,
        }
        if tools is not None:
            kwargs["tools"] = tools
        return self.client.chat.completions.create(**kwargs)

    def run(self, user_id: str, user_query: str, request_id: str) -> Dict[str, Any]:
        # seed: write the user message into memory on the server side
        self.gateway.mem_add(user_id, "user", user_query)

        # FIX 1: define server_ctx before using it
        server_ctx = self._compose_server_context(user_id, user_query, request_id)

        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "system", "content": server_ctx},
            {"role": "user", "content": json.dumps({"user_id": user_id, "request_id": request_id, "query": user_query})}
        ]

        job_submitted: Optional[str] = None

        # tool-calling loop
        for _ in range(8):
            resp = self._call(messages, tools=self.tools_spec)
            msg = resp.choices[0].message

            # If the model wants tools, execute them
            if msg.tool_calls:
                messages.append({
                    "role": msg.role,
                    "content": msg.content or "",
                    "tool_calls": [tc.model_dump() for tc in msg.tool_calls]
                })
                for tc in msg.tool_calls:
                    name = tc.function.name
                    args = json.loads(tc.function.arguments or "{}")
                    try:
                        result = self.dispatch[name](**args)
                        if name == "job_submit":
                            # track job id to report back
                            job_submitted = str(result.get("job_id"))
                    except Exception as e:
                        result = {"error": f"{type(e).__name__}: {e}"}

                    # append tool result
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tc.id,
                        "name": name,
                        "content": json.dumps(result)
                    })
                continue

            # No tool call → final answer
            final_text = msg.content or "(no response)"

            if self._looks_meta_reply(final_text):
                messages.append({"role": "system", "content": "Answer directly and self-contained. Do not say you've already answered."})
                messages.append({"role": "user", "content": user_query})
                resp2 = self._call(messages, tools=None if self.tools_supported is False else self.tools_spec)
                final_text = resp2.choices[0].message.content or final_text

            # persist assistant turn even if model forgot to call memory_write
            self.gateway.mem_add(user_id, "assistant", final_text)

            if job_submitted:
                job = self.gateway.job_status(job_submitted)
                status = job.status if job else "queued"
                return {"type": "job", "job_id": job_submitted, "status": status, "message": final_text}
            return {"type": "chat", "message": final_text}

        # Safety exit after many tool hops
        fallback = "I’ve applied available tools and context. Here’s my concise answer based on current knowledge."
        self.gateway.mem_add(user_id, "assistant", fallback)
        return {"type": "chat", "message": fallback}

    def _looks_meta_reply(self, text: str) -> bool:
        t = text.lower()
        bad = ["i've provided", "as mentioned above", "as stated earlier", "as i explained above"]
        return any(p in t for p in bad)
