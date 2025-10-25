"""
Pipeline bridge connecting unified API to LangGraph drug discovery workflow.
Handles job submission, tracking, and result retrieval.
"""

import asyncio
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime
import json
import uuid

# LangGraph workflow import
try:
    from fade.workflows.drug_discovery import run_drug_discovery_pipeline
    PIPELINE_AVAILABLE = True
except ImportError:
    PIPELINE_AVAILABLE = False
    print("[WARNING] LangGraph pipeline not available. Check fade/workflows/drug_discovery.py")


class PipelineBridge:
    """
    Bridge between unified API and LangGraph drug discovery pipeline.
    Manages job lifecycle and state tracking.
    """
    
    def __init__(self):
        """Initialize pipeline bridge."""
        self.jobs: Dict[str, Dict] = {}  # In-memory job tracking
        self.pipeline_available = PIPELINE_AVAILABLE
        
        # Create output directory for jobs
        self.output_dir = Path("output/jobs")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"[INFO] Pipeline bridge initialized. LangGraph available: {self.pipeline_available}")
    
    async def submit_job(
        self,
        job_id: str,
        query: str,
        user_id: str,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Submit job to LangGraph pipeline.
        
        Args:
            job_id: Unique job identifier
            query: Natural language drug discovery query
            user_id: User identifier
            **kwargs: Additional parameters
            
        Returns:
            Job submission status
        """
        # Create job record
        job = {
            "id": job_id,
            "query": query,
            "user_id": user_id,
            "status": "queued",
            "progress": 0,
            "message": "Job queued for LangGraph processing",
            "created_at": datetime.utcnow().isoformat(),
            "updated_at": datetime.utcnow().isoformat(),
            "completed_at": None,
            "results": None,
            "error": None
        }
        
        self.jobs[job_id] = job
        
        if self.pipeline_available:
            # Run pipeline in background
            asyncio.create_task(self._run_pipeline(job_id, query, user_id))
            
            return {
                "success": True,
                "message": "Job submitted to LangGraph pipeline",
                "job_id": job_id
            }
        else:
            # Simulate pipeline execution for testing
            asyncio.create_task(self._simulate_pipeline(job_id, query, user_id))
            
            return {
                "success": True,
                "message": "Job submitted (simulation mode - LangGraph not connected)",
                "job_id": job_id,
                "warning": "Running in simulation mode"
            }
    
    async def _run_pipeline(self, job_id: str, query: str, user_id: str):
        """
        Execute actual LangGraph pipeline.
        
        Args:
            job_id: Job identifier
            query: Drug discovery query
            user_id: User identifier
        """
        try:
            # Update status
            self._update_job_status(job_id, "running", 10, "Starting LangGraph workflow")
            
            # Run the actual pipeline
            result = await run_drug_discovery_pipeline(
                query=query,
                user_id=user_id,
                job_id=job_id,
                output_dir=self.output_dir / job_id
            )
            
            # Update with results
            self._update_job_status(
                job_id, 
                "completed", 
                100, 
                "LangGraph workflow completed successfully",
                results=result
            )
            
        except Exception as e:
            # Handle errors
            self._update_job_status(
                job_id,
                "failed",
                self.jobs[job_id]["progress"],
                f"Pipeline error: {str(e)}",
                error=str(e)
            )
            print(f"[ERROR] LangGraph pipeline failed for job {job_id}: {e}")
    
    async def _simulate_pipeline(self, job_id: str, query: str, user_id: str):
        """
        Simulate pipeline execution for testing.
        
        Args:
            job_id: Job identifier
            query: Drug discovery query
            user_id: User identifier
        """
        # Simulate pipeline stages
        stages = [
            (20, "Target Research: Identifying protein target"),
            (35, "Structure Resolution: Obtaining 3D structure"),
            (50, "Pocket Detection: Finding binding sites"),
            (65, "Molecule Generation: Creating candidates"),
            (80, "Screening: Evaluating ADMET properties"),
            (95, "Analysis: Compiling results"),
            (100, "Pipeline completed (simulation)")
        ]
        
        for progress, message in stages:
            await asyncio.sleep(2)  # Simulate processing time
            self._update_job_status(job_id, "running", progress, message)
        
        # Simulate results
        mock_results = {
            "target": "KRAS G12C" if "KRAS" in query else "Target protein",
            "molecules_generated": 5,
            "top_candidate": {
                "smiles": "CC(C)C1=CC=C(C=C1)C(=O)N",
                "binding_affinity": -8.5,
                "properties": {
                    "MW": 315.4,
                    "logP": 2.3,
                    "TPSA": 45.2
                }
            },
            "simulation_note": "This is simulated data - connect LangGraph for real results"
        }
        
        self._update_job_status(
            job_id,
            "completed",
            100,
            "Simulation completed",
            results=mock_results
        )
    
    def _update_job_status(
        self,
        job_id: str,
        status: str,
        progress: int,
        message: str,
        results: Optional[Dict] = None,
        error: Optional[str] = None
    ):
        """Update job status."""
        if job_id in self.jobs:
            self.jobs[job_id].update({
                "status": status,
                "progress": progress,
                "message": message,
                "updated_at": datetime.utcnow().isoformat(),
                "results": results,
                "error": error
            })
            
            if status == "completed" or status == "failed":
                self.jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
            
            # Save to file for persistence
            job_file = self.output_dir / f"{job_id}.json"
            with open(job_file, 'w') as f:
                json.dump(self.jobs[job_id], f, indent=2)
    
    def get_job_status(self, job_id: str) -> Optional[Dict]:
        """
        Get current job status.
        
        Args:
            job_id: Job identifier
            
        Returns:
            Job status dictionary or None
        """
        # Check memory first
        if job_id in self.jobs:
            return self.jobs[job_id]
        
        # Check file system
        job_file = self.output_dir / f"{job_id}.json"
        if job_file.exists():
            with open(job_file, 'r') as f:
                return json.load(f)
        
        return None
    
    def list_jobs(self, user_id: Optional[str] = None) -> list:
        """
        List all jobs or jobs for specific user.
        
        Args:
            user_id: Optional user filter
            
        Returns:
            List of job summaries
        """
        jobs = []
        
        # From memory
        for job in self.jobs.values():
            if user_id is None or job["user_id"] == user_id:
                jobs.append({
                    "id": job["id"],
                    "query": job["query"],
                    "status": job["status"],
                    "progress": job["progress"],
                    "created_at": job["created_at"]
                })
        
        return sorted(jobs, key=lambda x: x["created_at"], reverse=True)


# Global bridge instance
pipeline_bridge = PipelineBridge()
