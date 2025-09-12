from typing import Optional, Dict, Any
from pydantic import BaseModel

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

JsonObj = Dict[str, Any]
