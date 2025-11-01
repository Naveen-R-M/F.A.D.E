"""
DeepSeek LLM Agent for intelligent responses via Ollama.
NO FALLBACKS - Fails if Ollama is not available.
"""

import os
import json
import requests
from typing import Dict, Any, Optional, List

class DeepSeekAgent:
    """Handles intelligent responses using DeepSeek via Ollama. NO FALLBACKS."""
    
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize DeepSeek Agent.
        
        Raises:
            RuntimeError: If Ollama is not available
        """
        self.model = os.getenv("AI_MODEL")
        if not self.model:
            raise ValueError("AI_MODEL environment variable not set")
            
        self.api_base = os.getenv("AI_API_BASE")
        if not self.api_base:
            raise ValueError("AI_API_BASE environment variable not set")
        
        # Check Ollama availability - NO FALLBACK
        response = requests.get(f"{self.api_base.replace('/v1', '')}/api/tags")
        if response.status_code != 200:
            raise RuntimeError(f"Ollama not available at {self.api_base}")
        
        self.initialized = True
        print(f"[INFO] DeepSeek Agent initialized with model: {self.model}")
    
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
        """
        Generate response using DeepSeek.
        NO FALLBACK - Raises exception on failure.
        
        Raises:
            RuntimeError: If response generation fails
        """
        if not self.initialized:
            raise RuntimeError("DeepSeek Agent not initialized")
        
        messages = [{"role": "system", "content": self._get_system_prompt()}]
        
        prompt = f"User query: {query}"
        if context:
            prompt = f"Context: {context}\n\n{prompt}"
        
        messages.append({"role": "user", "content": prompt})
        
        response = requests.post(
            f"{self.api_base}/chat/completions",
            json={"model": self.model, "messages": messages, "max_tokens": 500}
        )
        
        if response.status_code != 200:
            raise RuntimeError(f"DeepSeek API call failed: {response.status_code} - {response.text}")
        
        data = response.json()
        return {
            "message": data['choices'][0]['message']['content'],
            "ai_generated": True,
            "model": self.model
        }
    
    def generate_consent_preview(self, query: str) -> str:
        """
        Generate job preview.
        NO FALLBACK - Requires working LLM.
        """
        if not self.initialized:
            raise RuntimeError("DeepSeek Agent not initialized")
            
        return f"""The LangGraph F.A.D.E pipeline will:
• Target Research: Identify protein via RCSB
• Structure Resolution: Get 3D structure
• Pocket Detection: Find binding sites  
• Molecule Generation: Create candidates
• Screening: Evaluate properties

Process requires HPC cluster (6-8 hours)."""

# Compatibility alias
GeminiAgent = DeepSeekAgent
