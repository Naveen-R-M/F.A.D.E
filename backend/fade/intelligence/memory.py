"""
Memory management for conversation history.
"""

import json
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime

class MemoryStore:
    """
    Manages conversation history for users.
    Ported from gateway/app/core/memory.py
    """
    
    def __init__(self, storage_path: Optional[Path] = None):
        """
        Initialize memory store.
        
        Args:
            storage_path: Optional path to store persistent memory (for future use)
        """
        # In-memory storage for now (like gateway)
        self._memory: Dict[str, List[Dict]] = {}
        self.storage_path = storage_path or Path("data/memory")
        
        # Create storage directory if using persistent storage
        if self.storage_path:
            self.storage_path.mkdir(parents=True, exist_ok=True)
    
    def add(self, user_id: str, role: str, content: str) -> None:
        """
        Add a message to user's conversation history.
        
        Args:
            user_id: User identifier
            role: Message role ('user' or 'assistant')
            content: Message content
        """
        if user_id not in self._memory:
            self._memory[user_id] = []
        
        self._memory[user_id].append({
            "role": role,
            "content": content,
            "timestamp": datetime.utcnow().isoformat()
        })
        
        # Keep only last 100 messages per user to manage memory
        if len(self._memory[user_id]) > 100:
            self._memory[user_id] = self._memory[user_id][-100:]
    
    def get_recent(self, user_id: str, k: int = 10) -> List[Dict]:
        """
        Get k most recent messages for a user.
        
        Args:
            user_id: User identifier
            k: Number of recent messages to return
            
        Returns:
            List of recent messages
        """
        if user_id not in self._memory:
            return []
        
        return self._memory[user_id][-k:]
    
    def get_all(self, user_id: str) -> List[Dict]:
        """
        Get all messages for a user.
        
        Args:
            user_id: User identifier
            
        Returns:
            All messages for the user
        """
        return self._memory.get(user_id, [])
    
    def clear(self, user_id: str) -> None:
        """
        Clear all messages for a user.
        
        Args:
            user_id: User identifier
        """
        if user_id in self._memory:
            del self._memory[user_id]
    
    def clear_all(self) -> None:
        """Clear all memory for all users."""
        self._memory.clear()
    
    def get_user_count(self) -> int:
        """Get number of users with stored conversations."""
        return len(self._memory)
    
    def get_message_count(self, user_id: str) -> int:
        """Get number of messages for a specific user."""
        return len(self._memory.get(user_id, []))
    
    def export_to_file(self, user_id: str, filepath: Path) -> None:
        """
        Export a user's conversation to a JSON file.
        
        Args:
            user_id: User identifier
            filepath: Path to save the JSON file
        """
        if user_id in self._memory:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(self._memory[user_id], f, indent=2)
    
    def import_from_file(self, user_id: str, filepath: Path) -> None:
        """
        Import a user's conversation from a JSON file.
        
        Args:
            user_id: User identifier
            filepath: Path to the JSON file
        """
        if filepath.exists():
            with open(filepath, 'r', encoding='utf-8') as f:
                self._memory[user_id] = json.load(f)
    
    def search_messages(self, user_id: str, query: str) -> List[Dict]:
        """
        Search through a user's messages for a query string.
        
        Args:
            user_id: User identifier
            query: Search query
            
        Returns:
            Messages containing the query
        """
        if user_id not in self._memory:
            return []
        
        query_lower = query.lower()
        results = []
        
        for message in self._memory[user_id]:
            if query_lower in message["content"].lower():
                results.append(message)
        
        return results
