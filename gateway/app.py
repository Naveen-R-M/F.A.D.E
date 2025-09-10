#!/usr/bin/env python3
"""
F.A.D.E Gateway - Hybrid Version (RAG + LLM classification)
- Parses sample_queries.txt (markdown) and project_context.txt
- Lightweight, dependency-free RAG (token overlap + fuzzy)
- LLM-driven classification (Gemini) with robust fallback heuristics
- RAG-informed answer generation
"""

import os
import re
import json
import difflib
import textwrap
import subprocess
import asyncio
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime

from fastapi import FastAPI, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Try to import Gemini for response generation
try:
    import google.generativeai as genai
    GEMINI_AVAILABLE = True
except ImportError:
    GEMINI_AVAILABLE = False


# === Data Models ===

class JobRequest(BaseModel):
    id: str
    query: str
    model: str
    current_time: str
    user_id: str

class JobStatus(BaseModel):
    id: str
    status: str
    message: Optional[str] = None
    error: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None


# === Fixed Hybrid Gateway with RAG + LLM ===

class FixedHybridGateway:
    def __init__(self):
        self.jobs: Dict[str, JobStatus] = {}
        self.active_processes: Dict[str, subprocess.Popen] = {}

        # Load project context (legacy helper) + full corpus for RAG
        self.project_context = self._load_project_context()
        self.corpus = self._load_corpus()

        # Initialize Gemini for classification + answer generation
        self.gemini_model = None
        self.gemini_working = False
        self._setup_gemini()

        # Backend path (if you later hook real pipelines)
        self.backend_path = Path(__file__).parent.parent / "backend"

        print("✅ Fixed Hybrid Gateway initialized")
        print(f"🤖 Gemini status: {'Working' if self.gemini_working else 'Not available - using fallback'}")
        print(f"📁 Backend path: {self.backend_path}")

    # -------------------------
    # Setup & file I/O helpers
    # -------------------------

    def _setup_gemini(self):
        """Setup and test Gemini"""
        if not GEMINI_AVAILABLE:
            print("⚠️ google-generativeai not installed")
            return

        api_key = os.getenv("GEMINI_API_KEY")
        if not api_key:
            print("⚠️ GEMINI_API_KEY not found in environment")
            return

        try:
            genai.configure(api_key=api_key)
            self.gemini_model = genai.GenerativeModel("gemini-2.5-pro")

            # Test Gemini with a simple query
            test_response = self.gemini_model.generate_content("What is 2+2?")
            if test_response and getattr(test_response, "text", None):
                self.gemini_working = True
                print("✅ Gemini initialized and tested successfully")
            else:
                print("❌ Gemini test failed - using fallback responses")

        except Exception as e:
            print(f"❌ Gemini setup failed: {e}")
            print("🔄 Will use fallback responses")

    def _load_project_context(self) -> str:
        """Legacy: Load project context (used as default if corpus missing)."""
        try:
            context_file = Path(__file__).parent / "project_context.txt"
            if context_file.exists():
                return context_file.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pass
        return "F.A.D.E is an AI-powered drug discovery platform."

    def _read_text_file(self, path: Path) -> str:
        try:
            if path.exists():
                return path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pass
        return ""

    # -------------------------
    # RAG utilities (stdlib)
    # -------------------------

    def _tokenize(self, text: str):
        return [t for t in re.findall(r"[A-Za-z0-9]+", text.lower()) if len(t) > 2]

    def _chunk(self, text: str, size: int = 700, overlap: int = 120):
        if not text:
            return []
        words = text.split()
        chunks, i = [], 0
        while i < len(words):
            chunk_words = words[i:i+size]
            chunks.append(" ".join(chunk_words))
            i += max(1, size - overlap)
        return chunks

    def _score(self, query: str, passage: str) -> float:
        # hybrid: token overlap + fuzzy character similarity
        qtok = set(self._tokenize(query))
        ptok = set(self._tokenize(passage))
        jacc = len(qtok & ptok) / (len(qtok | ptok) + 1e-9)
        fuzz = difflib.SequenceMatcher(None, query.lower(), passage.lower()).ratio()
        return 0.6 * jacc + 0.4 * fuzz

    def _parse_sample_queries_markdown(self, md: str):
        """
        Parse sample_queries.txt that uses Markdown sections:
          ## CHAT QUERIES ...
          ## JOB QUERIES ...
        Returns:
            labeled: list[(label, text)] where label in {"CHAT","JOB"}
            unlabeled: list[str]
        """
        labeled, unlabeled = [], []
        current_label = None
        for raw in md.splitlines():
            line = raw.strip()
            if not line:
                continue

            # Detect top-level group headers
            if re.match(r"^##\s*CHAT QUERIES", line, flags=re.I):
                current_label = "CHAT"
                continue
            if re.match(r"^##\s*JOB QUERIES", line, flags=re.I):
                current_label = "JOB"
                continue

            # Accept bullet/numbered items as examples
            if re.match(r"^[-*]\s+", line) or re.match(r"^\d+\.\s+", line):
                text = re.sub(r"^[-*]\s+|\d+\.\s+", "", line).strip()
                if text:
                    if current_label in {"CHAT", "JOB"}:
                        labeled.append((current_label, text))
                    else:
                        unlabeled.append(text)

        return labeled, unlabeled

    def _load_corpus(self):
        """Load RAG sources: project_context + parsed sample queries (CHAT/JOB)."""
        base = Path(__file__).parent
        project_path = base / "project_context.txt"
        samples_path = base / "sample_queries.txt"

        project_text = self._read_text_file(project_path)
        sample_text = self._read_text_file(samples_path)

        labeled, unlabeled = self._parse_sample_queries_markdown(sample_text)

        # Chunk project context for retrieval
        project_chunks = [(f"project_context.txt:{i}", c)
                          for i, c in enumerate(self._chunk(project_text))]

        # Treat each example line as a passage
        sample_chunks = [(f"sample_queries.txt:{i}", c)
                         for i, c in enumerate([t for _, t in labeled] + unlabeled)]

        return {
            "project_text": project_text,
            "project_chunks": project_chunks,
            "sample_labeled": labeled,     # list[(label, text)]
            "sample_unlabeled": unlabeled, # list[str]
            "sample_chunks": sample_chunks # list[(id, passage)]
        }

    def _retrieve(self, query: str, k: int = 6):
        passages = []
        passages.extend(self.corpus.get("project_chunks", []))
        passages.extend(self.corpus.get("sample_chunks", []))
        if not passages:
            return []
        scored = [(self._score(query, p), pid, p) for pid, p in passages]
        scored.sort(key=lambda x: x[0], reverse=True)
        return scored[:k]

    def _format_rag_block(self, query: str, k: int = 6) -> str:
        hits = self._retrieve(query, k=k)
        if not hits:
            return ""
        lines = []
        for rank, (score, pid, p) in enumerate(hits, 1):
            snippet = textwrap.shorten(p.replace("\n", " "), width=600, placeholder=" …")
            lines.append(f"[{rank}] ({pid}) {snippet}")
        return "\n".join(lines)

    def _extract_json(self, text: str) -> Optional[dict]:
        """Robust JSON extractor: try direct parse, then first {...} block."""
        try:
            return json.loads(text)
        except Exception:
            pass
        m = re.search(r"\{.*\}", text, flags=re.DOTALL)
        if not m:
            return None
        try:
            return json.loads(m.group(0))
        except Exception:
            return None

    # -------------------------
    # Classification (LLM + fallback)
    # -------------------------

    async def _classify_via_llm(self, query: str) -> Optional[Dict[str, Any]]:
        """LLM-driven classification using few-shots from sample_queries + RAG context."""
        if not self.gemini_working:
            return None

        jobs = [t for lab, t in self.corpus.get("sample_labeled", []) if lab == "JOB"][:3]
        chats = [t for lab, t in self.corpus.get("sample_labeled", []) if lab == "CHAT"][:3]

        examples = []
        for t in jobs:
            examples.append({"query": t, "decision": "JOB", "reason": "matches job example"})
        for t in chats:
            examples.append({"query": t, "decision": "CHAT", "reason": "matches chat example"})

        rag_block = self._format_rag_block(query, k=6)
        proj = self.corpus.get("project_text", "") or self.project_context

        prompt = f"""
You are the *Classifier Agent* for F.A.D.E. Decide whether a user query is:

- "JOB": user wants the platform to run a computational drug discovery workflow (e.g., design/generate/optimize molecules, docking runs, pipeline execution).
- "CHAT": user wants information or explanations (no workflow execution).

Return STRICT JSON only with keys: decision ("JOB"|"CHAT"), reason (string), confidence (0..1).

PROJECT CONTEXT (excerpt):
\"\"\"\n{textwrap.shorten(proj, width=2000, placeholder=" …")}\n\"\"\" 

RETRIEVED PASSAGES:
{rag_block if rag_block else "(none)"}

FEW-SHOT EXAMPLES:
{json.dumps(examples, ensure_ascii=False)}

Classify this query:
\"\"\"{query}\"\"\" 

JSON ONLY:
"""
        try:
            resp = self.gemini_model.generate_content(prompt)
            if not resp or not getattr(resp, "text", None):
                return None
            data = self._extract_json(resp.text.strip())
            if not data or data.get("decision") not in ("JOB", "CHAT"):
                return None
            return {
                "decision": data["decision"],
                "reasoning": data.get("reason", ""),
                "confidence": float(data.get("confidence", 0.7))
            }
        except Exception:
            return None

    def _classify_reliably(self, query: str) -> Dict[str, Any]:
        """Heuristic classification with small knowledge base + labeled nearest-neighbor boost."""
        query_lower = query.lower()
        print(f"🔍 Classifying: '{query}'")

        # Knowledge base (exact-ish matches)
        knowledge_responses = {
            "what is kras": """**KRAS Protein**

KRAS is a critical oncogene and GTPase enzyme that acts as a molecular switch controlling cell growth and proliferation. Key facts:

**Function & Role:**
• Controls cell growth, proliferation, and survival pathways
• Acts as molecular "on/off" switch between GTP-bound (active) and GDP-bound (inactive) states
• Essential for normal cell function, but dangerous when mutated

**Cancer Connection:**
• **Mutated in ~30% of all cancers** - most common oncogene
• **Pancreatic cancer**: 90% have KRAS mutations  
• **Lung adenocarcinoma**: 25-30% have KRAS mutations
• **Colorectal cancer**: 40-50% have KRAS mutations

**Key Mutations:**
• **G12D, G12C, G12V** - Most common and studied
• **Q61H, Q61R** - Also significant but less common

**Drug Development:**
• Historically "undruggable" due to protein structure
• **2021 breakthrough**: Sotorasib (AMG 510) - first KRAS G12C inhibitor  
• Active research continues for other mutations

F.A.D.E can help design inhibitors targeting specific KRAS mutations. Would you like to explore this?""",

            "tell me about protein targets": """**Protein Targets in Drug Discovery**

Protein targets are specific proteins that drugs are designed to interact with to treat diseases:

**Major Target Classes:**
• **Kinases** (35% of approved drugs) - EGFR, KRAS, PI3K, CDK4/6, ALK
• **G-Protein Coupled Receptors** (30%) - Dopamine, serotonin, adrenergic receptors
• **Ion Channels** (15%) - Calcium, sodium, potassium channels
• **Nuclear Receptors** (10%) - Estrogen, androgen, thyroid receptors  
• **Enzymes** (10%) - Proteases, phosphatases, metabolic enzymes

**Key Cancer Targets:**
• **EGFR** - Growth factor receptor (lung, colorectal cancer)
• **HER2** - Overexpressed in ~20% of breast cancers
• **KRAS** - GTPase mutated in 30% of cancers
• **TP53** - "Guardian of genome" (challenging to target)
**VEGF/VEGFR** - Angiogenesis pathway

**Target Validation Criteria:**
• **Druggability** - Can small molecules bind effectively?  
• **Disease relevance** - Is it causally involved in the disease?
• **Safety profile** - What are the side effects of inhibition?
• **Clinical accessibility** - Can drugs reach the target tissue?

F.A.D.E specializes in designing molecules for both established and emerging targets.""",

            "what is f.a.d.e": """**F.A.D.E: Fully Agentic Drug Engine**

F.A.D.E is an AI-powered drug discovery platform that transforms natural language queries into complete computational drug discovery workflows:

**Core Technology:**
• **Multi-Agent Architecture**
• **Natural Language Interface**
• **HPC Integration**
• **End-to-End Pipeline**

**Agents:** Target Selector, Structure Predictor, Molecule Generator, Evaluator, Docking Agent, Refiner

**What You Get:**
• **3-5 optimized drug candidates**
• **Binding affinity predictions** and interaction maps
• **Comprehensive property assessment** (ADMET, toxicity, bioavailability)
• **Technical data** (SMILES, SDF files) + natural language explanations
• **Complete workflow documentation** and result provenance

**Timeline:** Typically 6-8 hours for full analysis

Ready to start your drug discovery project?"""
        }

        # KB direct hits (CHAT)
        for key_phrase, response in knowledge_responses.items():
            if key_phrase in query_lower:
                print(f"✅ Knowledge base match: {key_phrase}")
                return {
                    "decision": "CHAT",
                    "response": response,
                    "confidence": 0.95,
                    "reasoning": "Knowledge base match"
                }

        # Strong job patterns
        job_patterns = [
            "find molecules", "design molecules", "find drugs", "design drugs",
            "find compounds", "design compounds", "find inhibitors", "design inhibitors",
            "molecules targeting", "drugs targeting", "compounds for", "inhibitors for",
            "create molecules", "develop drugs", "generate molecules",
            "target ", "inhibit ", "block ", "design "
        ]
        if any(p in query_lower for p in job_patterns):
            print("✅ Job pattern detected")
            return {
                "decision": "JOB",
                "response": "",
                "confidence": 0.9,
                "reasoning": "Drug discovery request (pattern)"
            }

        # Nearest-neighbor boost from labeled examples
        best_lab, best_score = None, 0.0
        for lab, txt in self.corpus.get("sample_labeled", []):
            s = self._score(query, txt)
            if s > best_score:
                best_score, best_lab = s, lab
        if best_lab and best_score >= 0.55:
            return {
                "decision": best_lab,
                "response": "" if best_lab == "JOB" else "I'll provide an intelligent response using AI...",
                "confidence": min(0.85 + (best_score - 0.55), 0.99),
                "reasoning": f"Nearest labeled sample ({best_lab}, score={best_score:.2f})"
            }

        # Default: CHAT needs intelligent response
        print("✅ General question - needs intelligent response")
        return {
            "decision": "CHAT",
            "response": "I'll provide an intelligent response using AI...",
            "confidence": 0.8,
            "reasoning": "General question"
        }

    # -------------------------
    # Generation (LLM + fallback)
    # -------------------------

    async def _generate_gemini_response(self, query: str) -> str:
        """Generate response using Gemini + lightweight RAG."""
        rag_block = self._format_rag_block(query, k=8)
        proj = self.corpus.get("project_text", "") or self.project_context

        prompt = f"""
You are F.A.D.E, an AI-powered drug discovery assistant.

USER QUESTION:
\"\"\"{query}\"\"\"

PROJECT CONTEXT (excerpt):
\"\"\"\n{textwrap.shorten(proj, width=2000, placeholder=" …")}\n\"\"\" 

RETRIEVED PASSAGES (most relevant first):
{rag_block if rag_block else "(no passages retrieved)"}

GUIDELINES:
- Use the retrieved passages and project context when relevant.
- Be accurate and concise; optionally cite passage numbers like [1], [2].
- If a detail isn't in context but is standard, state it plainly.
- If unsure, say so briefly and give best-practice guidance.

Answer:
"""
        response = self.gemini_model.generate_content(prompt)
        return response.text.strip() if response and getattr(response, "text", None) else "I couldn't generate a response."

    def _generate_fallback_response(self, query: str) -> str:
        """Generate fallback response when Gemini is not available."""
        q = query.lower()
        if any(k in q for k in ["spike protein", "covid", "sars-cov-2"]):
            return """**SARS-CoV-2 Spike Protein Information**

The SARS-CoV-2 spike protein is the primary target for COVID-19 vaccines and therapeutics:

**Key Components:**
• **S1 subunit** - Contains receptor binding domain (RBD)
• **S2 subunit** - Responsible for membrane fusion
• **Receptor Binding Domain (RBD)** - Binds to human ACE2 receptor
• **N-terminal Domain (NTD)** - Additional binding region

**Important Variants:**
• **Alpha, Beta, Gamma, Delta, Omicron** variants have mutations in spike protein
• **Key mutation sites**: E484K, N501Y, D614G, P681H/R

**Drug Discovery Targets:**
• RBD-ACE2 interaction inhibitors
• Fusion inhibitors targeting S2 subunit  
• Neutralizing antibodies against RBD

F.A.D.E can help design inhibitors targeting spike protein interactions. Would you like to explore this?"""

        elif "list" in q or "proteins" in q:
            return """I understand you're looking for protein information. As a drug discovery platform, F.A.D.E can help with:

**Common Drug Targets:**
• Cancer proteins (EGFR, KRAS, HER2, TP53)
• Infectious disease targets (viral proteins, bacterial enzymes)
• Neurological targets (neurotransmitter receptors, ion channels)
• Metabolic targets (kinases, nuclear receptors)

**F.A.D.E Capabilities:**
• Target identification and validation
• Structure prediction and analysis  
• Molecule generation and optimization
• Property prediction and filtering

If you specify a protein or disease area, I can dig in further — or F.A.D.E can design molecules targeting any protein you specify."""

        else:
            return f"""Happy to help with drug discovery.

F.A.D.E can:
• Answer questions about proteins, diseases, and drug mechanisms
• Design molecules targeting specific proteins
• Run complete computational workflows
• Provide detailed analysis of drug candidates

Regarding “{query}”, I can elaborate or help you launch a design pipeline if that’s your goal."""

    # -------------------------
    # Orchestration
    # -------------------------

    async def classify_and_respond(self, query: str) -> Dict[str, Any]:
        """RAG-aware classification + generation."""
        # 1) Try LLM classification using labeled examples + RAG
        llm_cls = await self._classify_via_llm(query)
        if llm_cls:
            if llm_cls["decision"] == "JOB":
                return {**llm_cls, "response": "", "decision": "JOB"}
            # decision == CHAT → generate answer
            if self.gemini_working:
                try:
                    answer = await self._generate_gemini_response(query)
                    return {
                        "decision": "CHAT",
                        "response": answer,
                        "confidence": min(llm_cls.get("confidence", 0.7) + 0.05, 1.0),
                        "reasoning": llm_cls.get("reasoning", "LLM classification")
                    }
                except Exception as e:
                    print(f"⚠️ Gemini failed during response: {e}")
                    return {
                        "decision": "CHAT",
                        "response": self._generate_fallback_response(query),
                        "confidence": 0.6,
                        "reasoning": "Gemini failed; fallback used"
                    }

        # 2) Heuristic fallback (uses sample labels too)
        heuristic = self._classify_reliably(query)
        if heuristic["decision"] == "JOB":
            return heuristic

        # 3) Generate answer (Gemini or fallback)
        if self.gemini_working:
            try:
                answer = await self._generate_gemini_response(query)
                heuristic["response"] = answer
                heuristic["confidence"] = min(heuristic.get("confidence", 0.7) + 0.05, 1.0)
                return heuristic
            except Exception as e:
                print(f"⚠️ Gemini failed during fallback generation: {e}")

        heuristic["response"] = self._generate_fallback_response(query)
        return heuristic

    # -------------------------
    # Job submission (simulated)
    # -------------------------

    async def submit_job(self, job_request: JobRequest, background_tasks: BackgroundTasks) -> JobStatus:
        """Submit job to backend pipeline (simulation for now)."""
        job = JobStatus(
            id=job_request.id,
            status="queued",
            message="Job queued for F.A.D.E pipeline",
            started_at=datetime.now().isoformat()
        )

        self.jobs[job.id] = job
        print(f"🚀 Job {job.id} created (simulation mode)")
        return job

    def get_job_status(self, job_id: str) -> Optional[JobStatus]:
        """Get job status"""
        return self.jobs.get(job_id)


# === FastAPI App ===

app = FastAPI(
    title="F.A.D.E Gateway - Fixed Hybrid",
    description="Reliable classification + RAG-informed responses",
    version="2.4"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize gateway
gateway = FixedHybridGateway()

# === Routes ===

@app.get("/")
def root():
    return {
        "message": "F.A.D.E Gateway v2.4 - Fixed Hybrid Intelligence (RAG + LLM)",
        "status": "healthy",
        "features": [
            "RAG from project_context.txt + sample_queries.txt",
            "LLM-driven classification (Gemini) with fallback",
            "Lightweight retrieval + robust JSON parsing"
        ]
    }

@app.post("/jobs")
async def handle_request(job_request: JobRequest, background_tasks: BackgroundTasks):
    print(f"\n📥 Request: {job_request.query}")

    # Classify and respond
    result = await gateway.classify_and_respond(job_request.query)

    if result["decision"] == "CHAT":
        print(f"💬 Chat response (confidence: {result['confidence']:.2f})")
        return {
            "type": "chat",
            "message": result["response"],
            "reasoning": result.get("reasoning", ""),
            "confidence": result.get("confidence", 0.5)
        }
    else:
        print(f"⚙️ Job submission (confidence: {result['confidence']:.2f})")
        job = await gateway.submit_job(job_request, background_tasks)
        return {
            "type": "job",
            "job_id": job.id,
            "status": job.status,
            "message": job.message,
            "confidence": result.get("confidence", 0.5)
        }

@app.get("/health")
def health():
    return {
        "status": "healthy",
        "classification": "llm_with_fallback",
        "response_generation": "gemini_with_rag" if gateway.gemini_working else "fallback_only",
        "gemini_working": gateway.gemini_working,
        "total_jobs": len(gateway.jobs)
    }


if __name__ == "__main__":
    import uvicorn

    print("🚀 Starting Fixed Hybrid F.A.D.E Gateway v2.4...")
    print("🔧 Gateway: http://localhost:8000")
    print("🏥 Health:  http://localhost:8000/health")

    uvicorn.run(app, host="0.0.0.0", port=8000)
