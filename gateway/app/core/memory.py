from datetime import datetime
from typing import Dict, List, Dict as TDict
from . import models
from ..config import MAX_TURNS_PER_USER, MAX_CHARS_PER_USER

class MemoryStore:
    """Per-user short transcript with char+turn caps."""
    def __init__(self):
        self._mem: Dict[str, List[TDict[str, str]]] = {}

    def add(self, user_id: str, role: str, content: str) -> None:
        if not user_id: return
        entry = {"role": role, "content": content, "ts": datetime.now().isoformat()}
        self._mem.setdefault(user_id, []).append(entry)
        # turn cap
        max_items = MAX_TURNS_PER_USER * 2
        if len(self._mem[user_id]) > max_items:
            self._mem[user_id] = self._mem[user_id][-max_items:]
        # char cap
        total = sum(len(e["content"]) for e in self._mem[user_id])
        while self._mem[user_id] and total > MAX_CHARS_PER_USER:
            dropped = self._mem[user_id].pop(0)
            total -= len(dropped["content"])

    def recent(self, user_id: str, max_turns: int = 6):
        convo = self._mem.get(user_id, [])
        return convo[-max_turns*2:] if convo else []

    def clear(self, user_id: str) -> None:
        if user_id in self._mem:
            del self._mem[user_id]
