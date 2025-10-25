"""
Query classifier for determining CHAT vs JOB queries.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional
from fade.intelligence.rag import RAG


class QueryClassifier:
    """
    Classifies queries as CHAT (informational) or JOB (computational).
    Uses pattern matching and RAG for classification.
    """
    
    def __init__(self, data_dir: Optional[Path] = None):
        """
        Initialize classifier with patterns and RAG.
        
        Args:
            data_dir: Directory containing data files
        """
        self.data_dir = data_dir or Path("data")
        self.rag = RAG(self.data_dir)
        
        # Define job-related keywords and patterns
        self.job_keywords = [
            # Action verbs
            "find", "design", "generate", "create", "discover",
            "develop", "synthesize", "optimize", "screen", "identify",
            
            # Target objects
            "molecules", "compounds", "inhibitors", "drugs", "ligands",
            "candidates", "modulators", "antagonists", "agonists",
            
            # Target specifications
            "targeting", "binding", "against", "for", "that bind",
            "that inhibit", "that target", "that modulate",
            
            # Protein/target mentions
            "kras", "egfr", "braf", "her2", "pd-1", "pd-l1",
            "bcl-2", "vegf", "alk", "met", "ros1", "jak2"
        ]
        
        # Chat-related patterns
        self.chat_patterns = [
            r"^what (is|are)",
            r"^how (does|do)",
            r"^explain",
            r"^tell me about",
            r"^describe",
            r"^why",
            r"^when",
            r"^can you explain",
            r"^help me understand"
        ]
        
        # Greeting patterns
        self.greeting_patterns = [
            r"^(hi|hello|hey)",
            r"^good (morning|afternoon|evening)",
            r"^thank",
            r"^thanks",
            r"^bye",
            r"^goodbye"
        ]
    
    def classify(self, query: str) -> Dict:
        """
        Classify a query as CHAT or JOB.
        
        Args:
            query: User query to classify
            
        Returns:
            Classification result with type, confidence, and reasoning
        """
        query_lower = query.lower().strip()
        
        # Check for greetings first
        if self._is_greeting(query_lower):
            return {
                'type': 'CHAT',
                'confidence': 0.99,
                'reasoning': 'Greeting or casual conversation'
            }
        
        # Check for clear chat patterns
        if self._matches_chat_patterns(query_lower):
            return {
                'type': 'CHAT',
                'confidence': 0.90,
                'reasoning': 'Informational or educational query pattern'
            }
        
        # Calculate job score based on keywords
        job_score = self._calculate_job_score(query_lower)
        
        # Strong job signal
        if job_score >= 3:
            confidence = min(0.85 + (job_score - 3) * 0.03, 0.99)
            return {
                'type': 'JOB',
                'confidence': confidence,
                'reasoning': f'Contains {job_score} drug discovery keywords'
            }
        
        # No job signal
        if job_score == 0:
            # Use RAG to check against sample queries
            rag_result = self.rag.classify_query(query)
            
            # If RAG is confident, use its result
            if rag_result['confidence'] > 0.7:
                return rag_result
            
            # Default to CHAT
            return {
                'type': 'CHAT',
                'confidence': 0.85,
                'reasoning': 'No drug discovery keywords detected'
            }
        
        # Ambiguous (1-2 job keywords)
        # Use RAG classification with pattern score as tiebreaker
        rag_result = self.rag.classify_query(query)
        
        if rag_result['type'] == 'JOB' and job_score >= 2:
            return {
                'type': 'JOB',
                'confidence': min(rag_result['confidence'] * 1.1, 0.90),
                'reasoning': f'RAG similarity + {job_score} job keywords'
            }
        elif rag_result['type'] == 'CHAT':
            return {
                'type': 'CHAT',
                'confidence': rag_result['confidence'],
                'reasoning': rag_result['reasoning']
            }
        else:
            # When uncertain, prefer CHAT to avoid unnecessary jobs
            return {
                'type': 'CHAT',
                'confidence': 0.60,
                'reasoning': 'Ambiguous query, defaulting to chat'
            }
    
    def _is_greeting(self, query: str) -> bool:
        """Check if query is a greeting."""
        for pattern in self.greeting_patterns:
            if re.match(pattern, query):
                return True
        return False
    
    def _matches_chat_patterns(self, query: str) -> bool:
        """Check if query matches chat patterns."""
        for pattern in self.chat_patterns:
            if re.search(pattern, query):
                return True
        return False
    
    def _calculate_job_score(self, query: str) -> int:
        """Calculate job score based on keyword presence."""
        score = 0
        words = query.split()
        
        for keyword in self.job_keywords:
            if keyword in query:
                score += 1
                # Extra weight for action verbs at the beginning
                if keyword in ["find", "design", "generate", "create"] and words and words[0] == keyword:
                    score += 1
        
        return score
