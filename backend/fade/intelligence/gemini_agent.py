"""
Phase 4.1: Create Gemini Agent for intelligent responses.
This integrates Google's Gemini AI for enhanced query handling.
"""

import os
import json
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime

# Gemini imports
try:
    import google.generativeai as genai
    GEMINI_AVAILABLE = True
except ImportError:
    GEMINI_AVAILABLE = False
    print("[WARNING] google-generativeai not installed. Install with: pip install google-generativeai")


class GeminiAgent:
    """
    Handles intelligent responses using Gemini AI with RAG and Memory.
    Integrated with LangGraph workflow for drug discovery.
    """
    
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize Gemini Agent.
        
        Args:
            api_key: Gemini API key (or from environment)
        """
        self.api_key = api_key or os.getenv("GEMINI_API_KEY")
        self.model = None
        self.initialized = False
        
        if not self.api_key:
            print("[WARNING] No GEMINI_API_KEY found. Agent will use fallback responses.")
            return
        
        if not GEMINI_AVAILABLE:
            print("[ERROR] google-generativeai package not installed")
            return
        
        try:
            # Configure Gemini
            genai.configure(api_key=self.api_key)
            
            # Initialize model with system instruction
            self.model = genai.GenerativeModel(
                'gemini-1.5-flash',
                system_instruction=self._get_system_prompt()
            )
            
            self.initialized = True
            print("[INFO] Gemini Agent initialized successfully")
            
        except Exception as e:
            print(f"[ERROR] Failed to initialize Gemini: {e}")
            self.initialized = False
    
    def _get_system_prompt(self) -> str:
        """Get system prompt for Gemini."""
        return """You are F.A.D.E, an AI drug discovery copilot powered by LangGraph workflows.

ABOUT F.A.D.E:
F.A.D.E (Fully Agentic Drug Engine) is a comprehensive drug discovery platform that uses LangGraph to orchestrate multiple AI agents. The system transforms natural language queries into complete drug discovery workflows.

KEY COMPONENTS:
1. LangGraph Workflow: Orchestrates the entire pipeline with state management
2. Target Research Agent: Identifies and validates protein targets
3. Structure Resolution Agent: Uses AlphaFold3 or RCSB PDB for 3D structures
4. Pocket Detection Agent: Identifies binding sites using fpocket
5. Molecule Generation Agent: Creates novel drug candidates using DiffSBDD
6. Screening Agent: Evaluates ADMET properties
7. Analysis Agent: Provides comprehensive results

RESPONSE GUIDELINES:
- For CHAT queries: Provide educational, accurate information about drug discovery
- For JOB queries: Explain what the LangGraph pipeline will do (but don't submit without consent)
- Use provided RAG context when relevant
- Be specific about the computational approach and timeline (6-8 hours on HPC)
- Mention LangGraph orchestration when discussing the pipeline
- Always be helpful and technically accurate

IMPORTANT:
- The system uses real computational tools, not simulations
- HPC cluster resources are required for jobs
- LangGraph ensures reliable state management across the workflow
- Each agent is specialized for its specific task
"""
    
    def generate_response(
        self,
        query: str,
        query_type: str,
        context: Optional[str] = None,
        memory: Optional[List[Dict]] = None
    ) -> Dict[str, Any]:
        """
        Generate intelligent response using Gemini.
        
        Args:
            query: User query
            query_type: "CHAT" or "JOB"
            context: RAG context
            memory: Conversation history
            
        Returns:
            Response with message and metadata
        """
        if not self.initialized:
            return self._fallback_response(query, query_type, context)
        
        try:
            # Build prompt
            prompt = self._build_prompt(query, query_type, context, memory)
            
            # Generate response
            response = self.model.generate_content(prompt)
            
            return {
                "message": response.text,
                "ai_generated": True,
                "model": "gemini-1.5-flash"
            }
            
        except Exception as e:
            print(f"[ERROR] Gemini generation failed: {e}")
            return self._fallback_response(query, query_type, context)
    
    def _build_prompt(
        self,
        query: str,
        query_type: str,
        context: Optional[str],
        memory: Optional[List[Dict]]
    ) -> str:
        """Build prompt for Gemini."""
        
        prompt_parts = []
        
        # Add memory context if available
        if memory and len(memory) > 0:
            memory_text = "\n".join([
                f"{m['role']}: {m['content']}" 
                for m in memory[-5:]  # Last 5 messages
            ])
            prompt_parts.append(f"Recent conversation:\n{memory_text}\n")
        
        # Add RAG context if available
        if context:
            prompt_parts.append(f"Relevant context:\n{context}\n")
        
        # Add query type hint
        if query_type == "JOB":
            prompt_parts.append("This is a JOB query requesting drug discovery computation.\n")
        else:
            prompt_parts.append("This is a CHAT query seeking information.\n")
        
        # Add the actual query
        prompt_parts.append(f"User query: {query}\n")
        
        # Add response instruction
        if query_type == "JOB":
            prompt_parts.append(
                "Provide a brief explanation of what the LangGraph pipeline will do for this request. "
                "Include the specific agents that will be used and mention it takes 6-8 hours on HPC."
            )
        else:
            prompt_parts.append(
                "Provide a helpful, accurate response. Use the context if relevant."
            )
        
        return "\n".join(prompt_parts)
    
    def _fallback_response(
        self,
        query: str,
        query_type: str,
        context: Optional[str]
    ) -> Dict[str, Any]:
        """Generate fallback response without AI."""
        
        if query_type == "JOB":
            message = f"""I'll process your request: "{query}"

The F.A.D.E LangGraph pipeline will:
â€¢ Target Research: Identify and validate the protein target
â€¢ Structure Resolution: Obtain 3D structure via AlphaFold3 or PDB
â€¢ Pocket Detection: Find binding sites using fpocket
â€¢ Molecule Generation: Create candidates with DiffSBDD
â€¢ Screening: Evaluate ADMET properties
â€¢ Analysis: Compile comprehensive results

This LangGraph-orchestrated workflow uses HPC resources and takes 6-8 hours."""
        
        else:
            # For CHAT, try to use context
            if context and "F.A.D.E" in context:
                # Extract first relevant sentence
                sentences = context.split('.')
                relevant = [s for s in sentences if 'F.A.D.E' in s or 'drug' in s]
                if relevant:
                    message = relevant[0].strip() + "."
                else:
                    message = "F.A.D.E is a LangGraph-powered drug discovery platform that designs molecules for specific protein targets."
            else:
                message = "F.A.D.E is a LangGraph-powered drug discovery platform that designs molecules for specific protein targets."
        
        return {
            "message": message,
            "ai_generated": False,
            "fallback": True
        }
    
    def generate_consent_preview(self, query: str) -> str:
        """
        Generate preview of what job will do.
        
        Args:
            query: Job query
            
        Returns:
            Preview message
        """
        if self.initialized:
            try:
                prompt = f"""
                Brief explanation (5 bullet points max) of what F.A.D.E's LangGraph pipeline will do for:
                {query}
                
                Include key constraints and mention it uses HPC resources.
                Be specific about which agents will be involved.
                """
                
                response = self.model.generate_content(prompt)
                return response.text
                
            except Exception as e:
                print(f"[ERROR] Failed to generate preview: {e}")
        
        # Fallback preview
        return f"""The LangGraph-orchestrated F.A.D.E pipeline will:
â€¢ Use Target Research agent to identify the protein
â€¢ Run Structure Resolution to get 3D structure  
â€¢ Apply Pocket Detection to find binding sites
â€¢ Generate molecules using DiffSBDD
â€¢ Screen candidates for drug-like properties

Process requires HPC cluster (6-8 hours)."""
