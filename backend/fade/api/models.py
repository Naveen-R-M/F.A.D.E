"""
Phase 3.1: Create API request/response models for unified backend.
These Pydantic models define the structure for API communication.
"""

from pathlib import Path
from typing import Optional, Literal, Dict, Any, List
from datetime import datetime
from pydantic import BaseModel, Field
import uuid


# Request Models
class JobRequest(BaseModel):
    """Unified request model for both chat and job queries."""
    query: str = Field(..., description="User's natural language query")
    user_id: str = Field(default="default", description="User identifier")
    id: str = Field(default_factory=lambda: str(uuid.uuid4()), description="Request ID")
    consent_given: bool = Field(default=False, description="Whether user has consented to job submission")
    current_time: datetime = Field(default_factory=datetime.utcnow)
    
    class Config:
        json_schema_extra = {
            "example": {
                "query": "Find inhibitors for KRAS G12C",
                "user_id": "user123",
                "consent_given": False
            }
        }


# Response Models
class ChatResponse(BaseModel):
    """Response for chat/informational queries."""
    type: Literal["chat"] = "chat"
    message: str = Field(..., description="Response message")
    reasoning: Optional[str] = Field(None, description="Reasoning for response")
    confidence: Optional[float] = Field(None, description="Confidence score")
    context_used: Optional[bool] = Field(False, description="Whether RAG context was used")
    
    class Config:
        json_schema_extra = {
            "example": {
                "type": "chat",
                "message": "F.A.D.E is a drug discovery engine...",
                "reasoning": "Informational query",
                "confidence": 0.95
            }
        }


class JobResponse(BaseModel):
    """Response for job submissions."""
    type: Literal["job"] = "job"
    job_id: str = Field(..., description="Job identifier")
    status: str = Field(..., description="Job status")
    message: str = Field(..., description="Status message")
    tracking_url: Optional[str] = Field(None, description="URL to track job progress")
    estimated_time: Optional[str] = Field(None, description="Estimated completion time")
    
    class Config:
        json_schema_extra = {
            "example": {
                "type": "job",
                "job_id": "job_12345",
                "status": "submitted",
                "message": "Drug discovery job submitted",
                "tracking_url": "/api/jobs/job_12345/status",
                "estimated_time": "6-8 hours"
            }
        }


class ConsentResponse(BaseModel):
    """Response requiring user consent before job submission."""
    type: Literal["consent_required"] = "consent_required"
    message: str = Field(..., description="Brief answer about what will be done")
    call_to_action: str = Field(..., description="Question asking for consent")
    job_preview: Dict[str, Any] = Field(..., description="Preview of job details")
    
    class Config:
        json_schema_extra = {
            "example": {
                "type": "consent_required",
                "message": "I can design KRAS G12C inhibitors using...",
                "call_to_action": "Shall I submit a F.A.D.E job to design these molecules?",
                "job_preview": {
                    "query": "Find inhibitors for KRAS G12C",
                    "estimated_time": "6-8 hours",
                    "computational_resources": "HPC cluster with GPU"
                }
            }
        }


class ErrorResponse(BaseModel):
    """Error response model."""
    type: Literal["error"] = "error"
    error: str = Field(..., description="Error message")
    details: Optional[str] = Field(None, description="Additional error details")
    
    class Config:
        json_schema_extra = {
            "example": {
                "type": "error",
                "error": "Invalid request",
                "details": "Query cannot be empty"
            }
        }


# Job Status Models
class JobStatus(BaseModel):
    """Job status information."""
    id: str = Field(..., description="Job identifier")
    status: str = Field(..., description="Current job status")
    progress: Optional[int] = Field(None, description="Progress percentage")
    message: Optional[str] = Field(None, description="Status message")
    created_at: datetime = Field(..., description="Job creation time")
    updated_at: datetime = Field(..., description="Last update time")
    completed_at: Optional[datetime] = Field(None, description="Completion time")
    results: Optional[Dict[str, Any]] = Field(None, description="Job results if completed")
    
    class Config:
        json_schema_extra = {
            "example": {
                "id": "job_12345",
                "status": "running",
                "progress": 45,
                "message": "Running structure prediction",
                "created_at": "2024-01-01T12:00:00Z",
                "updated_at": "2024-01-01T13:00:00Z"
            }
        }


# Memory Models
class MemoryOperation(BaseModel):
    """Memory operation request."""
    operation: Literal["clear", "export", "search"] = Field(..., description="Operation to perform")
    user_id: str = Field(..., description="User identifier")
    query: Optional[str] = Field(None, description="Search query for search operation")
    
    class Config:
        json_schema_extra = {
            "example": {
                "operation": "clear",
                "user_id": "user123"
            }
        }


class MemoryResponse(BaseModel):
    """Memory operation response."""
    status: str = Field(..., description="Operation status")
    message: str = Field(..., description="Result message")
    data: Optional[List[Dict[str, Any]]] = Field(None, description="Memory data if applicable")
    
    class Config:
        json_schema_extra = {
            "example": {
                "status": "ok",
                "message": "Memory cleared for user123"
            }
        }


# Health Check Models
class HealthStatus(BaseModel):
    """Health check response."""
    status: str = Field(..., description="Service status")
    service: str = Field(..., description="Service name")
    version: str = Field(..., description="API version")
    intelligence: Dict[str, bool] = Field(..., description="Intelligence component status")
    pipeline: bool = Field(..., description="Pipeline availability")
    
    class Config:
        json_schema_extra = {
            "example": {
                "status": "healthy",
                "service": "fade-unified-backend",
                "version": "4.0",
                "intelligence": {
                    "memory": True,
                    "rag": True,
                    "classifier": True
                },
                "pipeline": True
            }
        }
