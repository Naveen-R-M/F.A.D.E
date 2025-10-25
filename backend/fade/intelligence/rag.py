"""
RAG (Retrieval-Augmented Generation) system for context-aware responses.
"""

import json
from pathlib import Path
from typing import List, Dict, Optional
import re
from difflib import SequenceMatcher

class RAG:
    """
    Simple RAG system for retrieving relevant context.
    Ported from gateway/app/core/rag.py
    """
    
    def __init__(self, data_dir: Optional[Path] = None):
        """
        Initialize RAG system.
        
        Args:
            data_dir: Directory containing sample_queries.txt and project_context.txt
        """
        self.data_dir = data_dir or Path("data")
        self.sample_queries = []
        self.project_context = ""
        self.context_chunks = []
        
        # Load data files
        self._load_sample_queries()
        self._load_project_context()
    
    def _load_sample_queries(self) -> None:
        """Load sample queries from file."""
        queries_file = self.data_dir / "sample_queries.txt"
        
        if queries_file.exists():
            with open(queries_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            current_category = None
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Check if this is a category header
                if line.startswith('[') and line.endswith(']'):
                    current_category = line[1:-1]
                    continue
                
                # Parse query type and text
                if ':' in line:
                    query_type, query_text = line.split(':', 1)
                    self.sample_queries.append({
                        'type': query_type.strip(),
                        'text': query_text.strip(),
                        'category': current_category
                    })
        
        print(f"[INFO] Loaded {len(self.sample_queries)} sample queries")
    
    def _load_project_context(self) -> None:
        """Load project context from file."""
        context_file = self.data_dir / "project_context.txt"
        
        if context_file.exists():
            with open(context_file, 'r', encoding='utf-8') as f:
                self.project_context = f.read()
            
            # Split into chunks for retrieval
            # Split on double newlines or section headers
            chunks = re.split(r'\n\n+|(?=##\s)', self.project_context)
            
            for chunk in chunks:
                chunk = chunk.strip()
                if chunk and len(chunk) > 50:  # Ignore very small chunks
                    self.context_chunks.append({
                        'content': chunk,
                        'snippet': chunk[:200] + '...' if len(chunk) > 200 else chunk
                    })
        
        print(f"[INFO] Loaded {len(self.context_chunks)} context chunks")
    
    def search(self, query: str, k: int = 5) -> List[Dict]:
        """
        Search for relevant context based on query.
        
        Args:
            query: Search query
            k: Number of results to return
            
        Returns:
            List of relevant context snippets with scores
        """
        results = []
        query_lower = query.lower()
        
        # Search through sample queries
        for sample in self.sample_queries:
            score = self._calculate_similarity(query_lower, sample['text'].lower())
            if score > 0.3:  # Threshold for relevance
                results.append({
                    'type': 'sample_query',
                    'content': sample['text'],
                    'query_type': sample['type'],
                    'category': sample['category'],
                    'score': score,
                    'snippet': f"[{sample['type']}] {sample['text']}"
                })
        
        # Search through context chunks
        for chunk in self.context_chunks:
            score = self._calculate_similarity(query_lower, chunk['content'].lower())
            if score > 0.2:  # Lower threshold for context
                results.append({
                    'type': 'context',
                    'content': chunk['content'],
                    'score': score,
                    'snippet': chunk['snippet']
                })
        
        # Sort by score and return top k
        results.sort(key=lambda x: x['score'], reverse=True)
        return results[:k]
    
    def _calculate_similarity(self, text1: str, text2: str) -> float:
        """
        Calculate similarity between two texts.
        Uses both exact matching and fuzzy matching.
        """
        # Check for exact phrase matches
        if text1 in text2 or text2 in text1:
            return 0.9
        
        # Check for keyword overlap
        words1 = set(text1.split())
        words2 = set(text2.split())
        
        if not words1 or not words2:
            return 0.0
        
        # Jaccard similarity
        intersection = words1.intersection(words2)
        union = words1.union(words2)
        jaccard = len(intersection) / len(union) if union else 0
        
        # Sequence matching for fuzzy similarity
        sequence = SequenceMatcher(None, text1, text2).ratio()
        
        # Weighted combination
        return (jaccard * 0.6) + (sequence * 0.4)
    
    def classify_query(self, query: str) -> Dict:
        """
        Classify a query as CHAT or JOB based on sample queries.
        
        Args:
            query: Query to classify
            
        Returns:
            Classification result with type and confidence
        """
        # Search for similar queries
        similar = self.search(query, k=3)
        
        if not similar:
            return {
                'type': 'CHAT',
                'confidence': 0.5,
                'reasoning': 'No similar queries found'
            }
        
        # Count query types in similar results
        job_count = 0
        chat_count = 0
        
        for result in similar:
            if result.get('query_type') == 'JOB':
                job_count += result['score']
            elif result.get('query_type') == 'CHAT':
                chat_count += result['score']
        
        # Determine classification
        if job_count > chat_count:
            confidence = min(job_count / (job_count + chat_count), 0.95)
            return {
                'type': 'JOB',
                'confidence': confidence,
                'reasoning': f'Similar to {len([r for r in similar if r.get("query_type") == "JOB"])} job queries'
            }
        else:
            confidence = min(chat_count / (job_count + chat_count) if (job_count + chat_count) > 0 else 0.5, 0.95)
            return {
                'type': 'CHAT',
                'confidence': confidence,
                'reasoning': f'Similar to {len([r for r in similar if r.get("query_type") == "CHAT"])} chat queries'
            }
    
    def get_context_for_response(self, query: str) -> str:
        """
        Get relevant context to help answer a query.
        
        Args:
            query: User query
            
        Returns:
            Relevant context text
        """
        results = self.search(query, k=3)
        
        if not results:
            return ""
        
        context_parts = []
        for result in results:
            if result['type'] == 'context' and result['score'] > 0.3:
                context_parts.append(result['content'])
        
        return "\n\n".join(context_parts)
