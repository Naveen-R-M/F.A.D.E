from typing import Optional, Dict, Any
from datetime import datetime
from .models import JobStatus

class JobStore:
    def __init__(self):
        self._jobs: Dict[str, JobStatus] = {}

    def submit_sync(self, job_id: str, payload: Dict[str, Any]):
        job = JobStatus(
            id=job_id,
            status="queued",
            message="Job queued for F.A.D.E pipeline",
            started_at=datetime.now().isoformat()
        )
        self._jobs[job.id] = job
        return {"job_id": job_id, "status": "queued", "message": job.message}

    def get(self, job_id: str) -> Optional[JobStatus]:
        return self._jobs.get(job_id)
    
    def set_status(self, job_id: str, status: str, message: str = ""):
        job = self._jobs.get(job_id)
        if job:
            job.status = status
            job.message = message
        else:
            self._jobs[job_id] = JobStatus(id=job_id, status=status, message=message)
            