"""
DeepSeek LLM Agent for intelligent responses via Ollama.
"""

import os
import json
import requests
from typing import Dict, Any, Optional, List

class DeepSeekAgent:
    """Handles intelligent responses using DeepSeek via Ollama."""
    
    def __init__(self, api_key: Optional[str] = None):
        """Initialize DeepSeek Agent."""
        self.model = os.getenv("AI_MODEL", "deepseek-r1:8b")
        self.api_base = os.getenv("AI_API_BASE", "http://localhost:11434/v1")
        self.initialized = False
        
        try:
            response = requests.get(f"{self.api_base.replace('/v1', '')}/api/tags")
            if response.status_code == 200:
                self.initialized = True
                print(f"[INFO] DeepSeek Agent initialized with model: {self.model}")
        except Exception as e:
            print(f"[WARNING] Could not connect to Ollama: {e}")
    
    def _get_system_prompt(self) -> str:
        return """You are F.A.D.E, an AI drug discovery copilot powered by LangGraph.

F.A.D.E (Fully Agentic Drug Engine) uses LangGraph to orchestrate AI agents for drug discovery.
LangGraph is a framework for building stateful, multi-agent applications with graph-based orchestration.

For questions about F.A.D.E: Explain the drug discovery platform and its agents.
For questions about LangGraph: Explain it's a graph orchestration framework.
For drug queries: Explain the pipeline stages."""
    
    def generate_response(self, query: str, query_type: str, 
                         context: Optional[str] = None, 
                         memory: Optional[List[Dict]] = None) -> Dict[str, Any]:
        """Generate response using DeepSeek."""
        
        if not self.initialized:
            return self._fallback_response(query, query_type, context)
        
        try:
            messages = [{"role": "system", "content": self._get_system_prompt()}]
            
            prompt = f"User query: {query}"
            if context:
                prompt = f"Context: {context}\n\n{prompt}"
            
            messages.append({"role": "user", "content": prompt})
            
            response = requests.post(
                f"{self.api_base}/chat/completions",
                json={"model": self.model, "messages": messages, "max_tokens": 500}
            )
            
            if response.status_code == 200:
                data = response.json()
                return {
                    "message": data['choices'][0]['message']['content'],
                    "ai_generated": True,
                    "model": self.model
                }
        except Exception as e:
            print(f"[ERROR] DeepSeek failed: {e}")
        
        return self._fallback_response(query, query_type, context)
    
    def _fallback_response(self, query: str, query_type: str, 
                           context: Optional[str]) -> Dict[str, Any]:
        """Improved fallback responses."""
        query_lower = query.lower()
        
        if "langgraph" in query_lower:
            message = "LangGraph is a framework for building stateful, multi-agent applications with graph-based orchestration."
        elif "fade" in query_lower or "f.a.d.e" in query_lower:
            message = "F.A.D.E (Fully Agentic Drug Engine) is an AI-powered drug discovery platform using LangGraph orchestration."
        else:
            message = "F.A.D.E uses LangGraph to orchestrate drug discovery workflows."
        
        return {"message": message, "ai_generated": False, "fallback": True}
    
    def generate_consent_preview(self, query: str) -> str:
        """Generate job preview."""
        return f"""The LangGraph F.A.D.E pipeline will:
• Target Research: Identify protein via UniProt
• Structure Resolution: Get 3D structure
• Pocket Detection: Find binding sites  
• Molecule Generation: Create candidates
• Screening: Evaluate properties

Process requires HPC cluster (6-8 hours)."""

# Compatibility alias
GeminiAgent = DeepSeekAgent
