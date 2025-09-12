# app/agents/gemini_agent.py
from __future__ import annotations
import json
import time
from typing import Dict, Any, Optional

import google.generativeai as genai

from ..config import GEMINI_API_KEY, GEMINI_MODEL
from ..core.gateway import Gateway
from .gemini_tools import build_tools_and_dispatch

SYSTEM_PROMPT = """You are F.A.D.E, an AI drug discovery copilot.

POLICY:
- Always answer the user's question directly in THIS message (self-contained). 
  Never say you already answered, “as above,” or “I’ve provided an overview.”
- Use tools only when needed; otherwise answer from your own knowledge.
- If RAG returns weak/irrelevant passages, ignore it and answer from your knowledge.
- Keep answers concise, structured, and actionable.

TOOL USE & BUDGETS:
- Server context may already include recent memory and RAG candidates. Do NOT re-fetch them unless it materially improves the answer.
- You may call tools to refine context, but budget is limited and non-retryable with the same args:
  - At most ONE memory_read and ONE rag_search per request.
  - Suppress repeated calls with identical arguments.
- Job tools:
  - Never submit a job without explicit consent.
  - After a job is submitted, briefly confirm next steps AND still answer the user’s question when appropriate (don’t reply with submission text only).

JOB INTENT & CONSENT:
- Treat requests to **find/design/generate molecules** or mentioning **inhibitors/targeting** a protein (e.g., “find molecules targeting KRAS G12D”) as potential JOB intent.
- For such requests:
  1) Provide a brief, high-level answer (≤5 bullets) and note key constraints (e.g., BBB, MW, logP).
  2) Ask a single yes/no question: “Shall I submit a F.A.D.E job now?” Do NOT submit without consent.
  3) Only if the user clearly consents (e.g., “yes”, “submit”, “run it”, “go ahead”, “proceed”, “start the job”), call:
     job_submit(job_id=<request.id>, payload={"prompt": <original user query>})
- Greetings/smalltalk (e.g., “hi”, “hello”, “thanks”) are NOT consent and should be answered normally with no job submission.
- If consent is ambiguous or the tool responds with consent_required/cta_missing, ask for explicit yes/no and DO NOT call job_submit again in the same turn.
- Don’t repeatedly ask to submit; if the user already declined, wait until they indicate new intent or add constraints.

OUTPUT STYLE:
- Be direct, precise, and avoid verbosity.
- When proposing a job, include a one-line payload summary and invite optional constraints (e.g., scaffold, MW, logP, pKa, off-targets).
"""

class GeminiToolAgent:
    """
    Manual tool-calling loop using google-generativeai Chat sessions:
      1) chat.send_message(initial_content)
      2) If function_call(s): dispatch locally
      3) Send results back via chat.send_message([...function_response parts...])
      4) Repeat until final text
      5) Log latency per step and total
    """

    def __init__(self, gateway: Gateway):
        if not GEMINI_API_KEY:
            raise RuntimeError("Missing GEMINI_API_KEY.")
        genai.configure(api_key=GEMINI_API_KEY)

        self.model_name = GEMINI_MODEL or "gemini-1.5-flash"
        self.gateway = gateway

        # Mutable state the tools can write into (e.g., last_job_id)
        self.tool_state: Dict[str, Any] = {}

        # Build tool spec + dispatch (already compatible with google-generativeai)
        self.tools_spec, self.dispatch = build_tools_and_dispatch(gateway, tool_state=self.tool_state)

        # Create a model bound to system prompt + our tools
        self.model = genai.GenerativeModel(
            model_name=self.model_name,
            system_instruction=SYSTEM_PROMPT,
            tools=self.tools_spec,
        )

    # ---------- Helpers ----------

    def _compose_server_context(self, user_id: str, user_query: str, request_id: str) -> str:
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

    def _fallback_manual(self, user_id: str, user_query: str, request_id: str) -> str:
        server_ctx = self._compose_server_context(user_id, user_query, request_id)
        prompt = (
            f"{SYSTEM_PROMPT}\n\n"
            f"{server_ctx}\n\n"
            f"User query:\n{user_query}\n\n"
            "Answer directly and self-contained."
        )
        t0 = time.perf_counter()
        resp = self.model.generate_content(prompt)
        t1 = time.perf_counter()
        print(f"[Gemini] fallback call latency: {int((t1 - t0)*1000)} ms")
        return resp.text or "(no response)"

    def _looks_meta_reply(self, text: str) -> bool:
        t = (text or "").lower()
        bad = ["i've provided", "as mentioned above", "as stated earlier", "as i explained above"]
        return any(p in t for p in bad)

    def _is_job_intent(self, text: str) -> bool:
        q = (text or "").lower()
        patterns = [
            "find molecules", "design molecules", "generate molecules",
            "find inhibitors", "design inhibitors",
            "targeting", "inhibitors for", "optimize molecules",
            "bbb permeability", "blood-brain barrier", "bbb"
        ]
        return any(p in q for p in patterns)
    
    def _looks_like_greeting(self, text: str) -> bool:
        q = (text or "").strip().lower()
        return q in {"hi", "hello", "hey", "howdy", "good morning", "good afternoon", "good evening"}

    def _has_affirmative_consent(self, text: str) -> bool:
        q = (text or "").lower()
        yes = ["yes", "yep", "yeah", "sure", "please do", "do it", "go ahead",
            "run it", "start it", "proceed", "submit", "ok submit", "okay submit", "start the job"]
        no  = ["no", "not now", "later", "stop", "don't", "do not", "cancel"]
        if any(n in q for n in no):
            return False
        return any(y in q for y in yes)

    def _prior_assistant_cta(self, user_id: str) -> bool:
        # Look back a few assistant turns for a CTA we generate
        turns = self.gateway.mem_recent(user_id, 6)
        for t in reversed(turns):
            if t["role"] != "assistant":
                continue
            msg = t.get("content", "").lower()
            if "would you like me to submit a f.a.d.e job" in msg:
                return True
        return False    

    # ---------- Function-call arg normalization ----------

    def _coerce_args(self, call) -> Any:
        """
        Gemini may return function_call.args as proto-plus MapComposite/ListComposite.
        Convert recursively to plain Python types for dispatch & dedupe.
        """
        def to_plain(x):
            # dict-like (including MapComposite)
            if hasattr(x, "items"):
                return {str(k): to_plain(v) for k, v in x.items()}
            # list/tuple-like (including ListComposite)
            if isinstance(x, (list, tuple)):
                return [to_plain(v) for v in x]
            # structs that expose .to_dict()
            if hasattr(x, "to_dict") and callable(getattr(x, "to_dict")):
                return to_plain(x.to_dict())
            return x

        try:
            if hasattr(call, "args") and call.args is not None:
                return to_plain(call.args)
            if hasattr(call, "args_json") and call.args_json:
                return json.loads(call.args_json)
        except Exception:
            pass
        return {}
    
    # ---------- Main entry ----------

    def run(self, user_id: str, user_query: str, request_id: str) -> Dict[str, Any]:
        MAX_TOOL_STEPS = 2            # keep it tight
        REQUEST_DEADLINE_MS = 45000   # total time budget
        seen_calls = set()            # (name, json_args)
        
        total_start = time.perf_counter()

        # Record user turn
        self.gateway.mem_add(user_id, "user", user_query)
        final_text = ""

        # Initial user content includes server context
        server_ctx = self._compose_server_context(user_id, user_query, request_id)
        payload = {"user_id": user_id, "request_id": request_id, "query": user_query}
        initial_content = f"{server_ctx}\n\nUSER_PAYLOAD_JSON:\n{json.dumps(payload)}\n\nUSER_QUERY:\n{user_query}"

        try:
            # Start a chat session (keeps history for tool loops)
            chat = self.model.start_chat()

            # First generation
            t0 = time.perf_counter()
            resp = chat.send_message(initial_content)
            print(f"[Gemini] gen step 1 latency: {int((time.perf_counter() - t0) * 1000)} ms")

            # Manual tool-calling loop (max 8 hops)
            for step in range(2, 2 + MAX_TOOL_STEPS):
                if (time.perf_counter() - total_start) * 1000 > REQUEST_DEADLINE_MS:
                    print("[Gemini] deadline exceeded, finalizing.")
                    break

                function_calls = []
                for cand in getattr(resp, "candidates", []) or []:
                    content = getattr(cand, "content", None)
                    if not content or not getattr(content, "parts", None):
                        continue
                    for part in content.parts:
                        if hasattr(part, "function_call") and part.function_call:
                            function_calls.append(part.function_call)

                if not function_calls:
                    final_text = resp.text or "(no response)"
                    break

                tool_parts = []
                for call in function_calls:
                    name = call.name
                    args = self._coerce_args(call)  # <-- convert to plain dict
                    # robust dedupe signature; default=str avoids proto dumps
                    sig = (name, json.dumps(args, sort_keys=True, default=str))
                    if sig in seen_calls:
                        print(f"[Gemini] duplicate tool call suppressed: {sig}")
                        continue
                    seen_calls.add(sig)

                    try:
                        # --- HARD CONSENT GATE FOR JOB SUBMISSION ---
                        if name == "job_submit":
                            # Only allow if user clearly consented in THIS turn AND (either prior CTA exists or text itself is job-intent)
                            if not self._has_affirmative_consent(user_query) or self._looks_like_greeting(user_query):
                                result = {
                                    "ok": False,
                                    "error": "consent_required",
                                    "reason": "Explicit consent required. Ask the user to confirm before submitting.",
                                }
                                tool_parts.append({
                                    "function_response": {
                                        "name": name,
                                        "response": {
                                            "name": name,
                                            "content": [{"text": json.dumps(result, ensure_ascii=False)}]
                                        }
                                    }
                                })
                                continue
                            # Optional: require that we recently asked a CTA
                            if not self._prior_assistant_cta(user_id) and not self._is_job_intent(user_query):
                                result = {
                                    "ok": False,
                                    "error": "cta_missing",
                                    "reason": "No prior CTA detected. Ask for confirmation first, then submit.",
                                }
                                tool_parts.append({
                                    "function_response": {
                                        "name": name,
                                        "response": {
                                            "name": name,
                                            "content": [{"text": json.dumps(result, ensure_ascii=False)}]
                                        }
                                    }
                                })
                                continue
                        start_tool = time.perf_counter()
                        result = self.dispatch[name](**args)
                        elapsed_ms = int((time.perf_counter() - start_tool) * 1000)
                        print(f"[Gemini] tool {name} exec: {elapsed_ms} ms")

                        if name == "job_submit" and isinstance(result, dict) and result.get("job_id"):
                            self.tool_state["last_job_id"] = result["job_id"]
                    except Exception as e:
                        elapsed_ms = int((time.perf_counter() - start_tool) * 1000)
                        print(f"[Gemini] tool {name} exec error after {elapsed_ms} ms: {type(e).__name__}: {e}")
                        result = {"error": f"{type(e).__name__}: {e}"}


                    # IMPORTANT: wrap the result as text content inside function_response
                    tool_parts.append({
                        "function_response": {
                            "name": name,
                            "response": {
                                "name": name,
                                "content": [{"text": json.dumps(result, ensure_ascii=False)}]
                            }
                        }
                    })
                
                if not tool_parts:
                    final_text = resp.text or "Proceeding with available information."
                    break

                t_tool = time.perf_counter()
                resp = chat.send_message(tool_parts)
                print(f"[Gemini] gen step {step} latency: {int((time.perf_counter() - t_tool) * 1000)} ms")
            else:
                # Loop exhausted
                final_text = resp.text or "Here’s a concise answer based on available context."

        except Exception as e:
            from google.api_core.exceptions import ResourceExhausted
            print(f"[Gemini error] {type(e).__name__}: {e}")
            if isinstance(e, ResourceExhausted):
                final_text = (
                    "I’m temporarily rate-limited by the Gemini free tier. "
                    "Please retry in ~30 seconds or switch to a lighter route. "
                    "I’ve kept your request in context."
                )
            else:
                final_text = self._fallback_manual(user_id, user_query, request_id)

        # Guard against meta reply once
        if self._looks_meta_reply(final_text):
            fix_prompt = (
                f"{SYSTEM_PROMPT}\n\n{server_ctx}\n\nUser query:\n{user_query}\n\n"
                "Answer now, directly and self-contained. Do NOT refer to a previous answer."
            )
            t_fix = time.perf_counter()
            try:
                resp2 = self.model.generate_content(fix_prompt)
                print(f"[Gemini] meta-guard latency: {int((time.perf_counter() - t_fix)*1000)} ms")
                final_text = resp2.text or final_text
            except Exception:
                pass

        # Persist assistant turn
        self.gateway.mem_add(user_id, "assistant", final_text)
        # If no job was submitted and this looks like a JOB intent, append a CTA
        if not self.tool_state.get("last_job_id") and self._is_job_intent(user_query):
            final_text = final_text.rstrip() + (
                f"\n\n**Would you like me to submit a F.A.D.E job now "
                f"(job_id: {request_id}) with BBB permeability as a hard constraint?** "
                f"Reply **yes** to proceed, or add constraints (e.g., scaffold, MW, logP, pKa, off-targets)."
            )

        # Persist assistant turn
        self.gateway.mem_add(user_id, "assistant", final_text)

        total_ms = int((time.perf_counter() - total_start) * 1000)
        print(f"[Gemini] total request latency: {total_ms} ms")

        # If a job was submitted, surface that
        job_id = self.tool_state.get("last_job_id")
        if job_id:
            job = self.gateway.job_status(str(job_id))
            status = job.status if job else "submitted"
            return {"type": "job", "job_id": str(job_id), "status": status, "message": final_text, "latency_ms": total_ms}

        return {"type": "chat", "message": final_text, "latency_ms": total_ms}
