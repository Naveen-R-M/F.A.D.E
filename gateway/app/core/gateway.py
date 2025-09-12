from pathlib import Path
from typing import Optional
from .memory import MemoryStore
from .rag import RAG
from .jobs import JobStore
from .models import JobRequest, JobStatus
from dataclasses import dataclass, asdict
from typing import Dict, Optional
from datetime import datetime

@dataclass
class JobRec:
    id: str
    status: str
    message: str = ""
    cwd: str = ""
    created_at: str = ""
    updated_at: str = ""

class Gateway:
    """Holds shared services: memory, rag, jobs."""
    def __init__(self, base_dir: Path):
        self._jobs: Dict[str, JobRec] = {}
        self.memory = MemoryStore()
        self.rag = RAG(base_dir)
        self.jobs = JobStore()

    # Convenience wrappers
    def mem_add(self, user_id: str, role: str, content: str): self.memory.add(user_id, role, content)
    def mem_recent(self, user_id: str, k: int): return self.memory.recent(user_id, k)
    def mem_clear(self, user_id: str): self.memory.clear(user_id)

    def rag_search(self, query: str, k: int): return self.rag.retrieve(query, k)

    async def submit_job(self, jr: JobRequest) -> JobStatus:
        job = JobStatus(id=jr.id, status="queued", message="Job queued", started_at=jr.current_time)
        # You can expand to real async submission here
        return job

    def submit_job_sync(self, job_id: str, payload: dict):
        return self.jobs.submit_sync(job_id, payload)

    def job_upsert(self, job_id: str, status: str, message: str = "", cwd: str = "") -> JobRec:
        now = datetime.utcnow().isoformat()
        rec = self._jobs.get(job_id)
        if rec:
            rec.status = status
            rec.message = message or rec.message
            rec.cwd = cwd or rec.cwd
            rec.updated_at = now
        else:
            rec = JobRec(id=job_id, status=status, message=message, cwd=cwd, created_at=now, updated_at=now)
            self._jobs[job_id] = rec
        return rec

    def job_status(self, job_id: str) -> Optional[JobRec]:
        return self._jobs.get(job_id)
