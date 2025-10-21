"""
Base Agent class for all drug discovery agents.

This provides common functionality like logging, error handling,
state validation, and checkpoint/resume capabilities.
"""

import logging
import time
import traceback
from abc import ABC, abstractmethod
from datetime import datetime
from typing import Dict, Any, Optional, List
from pathlib import Path
import json

from fade.state.drug_discovery_state import DrugDiscoveryState, ExecutionLogEntry


class BaseAgent(ABC):
    """
    Abstract base class for all agents in the drug discovery pipeline.
    
    Each agent should:
    1. Take the current state as input
    2. Perform its specific task
    3. Update the state with results
    4. Return the updated state
    """
    
    def __init__(self, 
                 name: str,
                 description: str,
                 checkpoint_dir: Optional[Path] = None,
                 max_retries: int = 3,
                 timeout: Optional[float] = None):
        """
        Initialize the base agent.
        
        Args:
            name: Agent identifier (e.g., "target_research", "structure_resolver")
            description: Human-readable description of agent's purpose
            checkpoint_dir: Directory for saving checkpoints
            max_retries: Maximum number of retries on failure
            timeout: Maximum execution time in seconds
        """
        self.name = name
        self.description = description
        self.checkpoint_dir = Path(checkpoint_dir) if checkpoint_dir else None
        self.max_retries = max_retries
        self.timeout = timeout
        
        # Set up logging
        self.logger = logging.getLogger(f"fade.agent.{name}")
        
    def __call__(self, state: DrugDiscoveryState) -> DrugDiscoveryState:
        """
        Main execution method with error handling and logging.
        
        Args:
            state: Current pipeline state
            
        Returns:
            Updated pipeline state
        """
        start_time = time.time()
        
        # Log execution start
        self.logger.info(f"Starting {self.name} agent")
        state = self._add_execution_log(state, "started")
        
        # Validate input state
        validation_result = self.validate_input(state)
        if not validation_result["is_valid"]:
            error_msg = f"Invalid input state: {', '.join(validation_result['errors'])}"
            self.logger.error(error_msg)
            state = self._add_error(state, error_msg)
            state["should_continue"] = False
            return state
            
        # Handle warnings
        if validation_result.get("warnings"):
            for warning in validation_result["warnings"]:
                self.logger.warning(warning)
                state = self._add_warning(state, warning)
        
        # Execute with retries
        retry_count = 0
        last_error = None
        
        while retry_count < self.max_retries:
            try:
                # Check for existing checkpoint
                if self.checkpoint_dir:
                    checkpoint_state = self._load_checkpoint(state["run_id"])
                    if checkpoint_state:
                        self.logger.info(f"Resuming from checkpoint for run {state['run_id']}")
                        state = checkpoint_state
                
                # Execute the agent's main logic
                state = self.execute(state)
                
                # Save checkpoint on successful execution
                if self.checkpoint_dir:
                    self._save_checkpoint(state)
                
                # Log successful completion
                execution_time = time.time() - start_time
                self.logger.info(f"Completed {self.name} agent in {execution_time:.2f}s")
                state = self._add_execution_log(state, "completed", 
                                               data={"execution_time": execution_time})
                
                # Update execution times
                if state.get("execution_times") is None:
                    state["execution_times"] = {}
                state["execution_times"][self.name] = execution_time
                
                return state
                
            except Exception as e:
                retry_count += 1
                last_error = e
                error_trace = traceback.format_exc()
                
                self.logger.error(f"Error in {self.name} agent (attempt {retry_count}/{self.max_retries}): {str(e)}")
                self.logger.debug(f"Full traceback: {error_trace}")
                
                if retry_count < self.max_retries:
                    wait_time = 2 ** retry_count  # Exponential backoff
                    self.logger.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                    
                    # Allow agent to handle error and potentially modify state
                    state = self.handle_error(state, e, retry_count)
        
        # All retries exhausted
        error_msg = f"Agent {self.name} failed after {self.max_retries} attempts: {str(last_error)}"
        self.logger.error(error_msg)
        state = self._add_error(state, error_msg)
        state = self._add_execution_log(state, "failed", message=error_msg)
        state["should_continue"] = False
        
        return state
    
    @abstractmethod
    def execute(self, state: DrugDiscoveryState) -> DrugDiscoveryState:
        """
        Execute the agent's main logic.
        
        This method should be implemented by each specific agent.
        
        Args:
            state: Current pipeline state
            
        Returns:
            Updated pipeline state
        """
        pass
    
    @abstractmethod
    def validate_input(self, state: DrugDiscoveryState) -> Dict[str, Any]:
        """
        Validate that the input state contains required information.
        
        This method should be implemented by each specific agent.
        
        Args:
            state: Current pipeline state
            
        Returns:
            Validation result with keys:
            - is_valid: bool
            - errors: List[str]
            - warnings: List[str]
        """
        pass
    
    def handle_error(self, state: DrugDiscoveryState, error: Exception, 
                    retry_count: int) -> DrugDiscoveryState:
        """
        Handle errors and potentially modify state before retry.
        
        Can be overridden by specific agents for custom error handling.
        
        Args:
            state: Current pipeline state
            error: The exception that occurred
            retry_count: Current retry attempt number
            
        Returns:
            Potentially modified state
        """
        # Default implementation just returns state unchanged
        return state
    
    def _add_execution_log(self, state: DrugDiscoveryState, status: str, 
                          message: Optional[str] = None, 
                          data: Optional[Dict[str, Any]] = None) -> DrugDiscoveryState:
        """Add an entry to the execution log."""
        if state.get("execution_log") is None:
            state["execution_log"] = []
            
        log_entry: ExecutionLogEntry = {
            "timestamp": datetime.now(),
            "agent": self.name,
            "action": "execute",
            "status": status,
            "message": message,
            "data": data
        }
        
        state["execution_log"].append(log_entry)
        return state
    
    def _add_error(self, state: DrugDiscoveryState, error_message: str) -> DrugDiscoveryState:
        """Add an error message to the state."""
        if state.get("error_messages") is None:
            state["error_messages"] = []
        state["error_messages"].append(f"[{self.name}] {error_message}")
        return state
    
    def _add_warning(self, state: DrugDiscoveryState, warning_message: str) -> DrugDiscoveryState:
        """Add a warning message to the state."""
        if state.get("warnings") is None:
            state["warnings"] = []
        state["warnings"].append(f"[{self.name}] {warning_message}")
        return state
    
    def _save_checkpoint(self, state: DrugDiscoveryState) -> None:
        """Save a checkpoint of the current state."""
        if not self.checkpoint_dir:
            return
            
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        checkpoint_path = self.checkpoint_dir / f"{state['run_id']}_{self.name}.json"
        
        # Convert state to JSON-serializable format
        # Note: This is simplified - in production you'd need proper serialization
        serializable_state = self._make_serializable(state)
        
        with open(checkpoint_path, 'w') as f:
            json.dump(serializable_state, f, indent=2)
            
        self.logger.debug(f"Saved checkpoint to {checkpoint_path}")
    
    def _load_checkpoint(self, run_id: str) -> Optional[DrugDiscoveryState]:
        """Load a checkpoint for a given run ID."""
        if not self.checkpoint_dir:
            return None
            
        checkpoint_path = self.checkpoint_dir / f"{run_id}_{self.name}.json"
        
        if not checkpoint_path.exists():
            return None
            
        with open(checkpoint_path, 'r') as f:
            state_dict = json.load(f)
            
        # Convert back to proper types
        # Note: This is simplified - in production you'd need proper deserialization
        state = self._deserialize_state(state_dict)
        
        return state
    
    def _make_serializable(self, obj: Any) -> Any:
        """Convert objects to JSON-serializable format."""
        if isinstance(obj, datetime):
            return obj.isoformat()
        elif isinstance(obj, Path):
            return str(obj)
        elif isinstance(obj, dict):
            return {k: self._make_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_serializable(item) for item in obj]
        else:
            return obj
    
    def _deserialize_state(self, state_dict: Dict[str, Any]) -> DrugDiscoveryState:
        """Convert JSON data back to proper types."""
        # Convert ISO strings back to datetime objects
        if state_dict.get("timestamp"):
            state_dict["timestamp"] = datetime.fromisoformat(state_dict["timestamp"])
            
        if state_dict.get("execution_log"):
            for log_entry in state_dict["execution_log"]:
                if log_entry.get("timestamp"):
                    log_entry["timestamp"] = datetime.fromisoformat(log_entry["timestamp"])
                    
        return state_dict


class AgentRegistry:
    """Registry for managing all agents in the pipeline."""
    
    _agents: Dict[str, BaseAgent] = {}
    
    @classmethod
    def register(cls, agent: BaseAgent) -> None:
        """Register an agent."""
        cls._agents[agent.name] = agent
        
    @classmethod
    def get(cls, name: str) -> Optional[BaseAgent]:
        """Get an agent by name."""
        return cls._agents.get(name)
    
    @classmethod
    def list_agents(cls) -> List[str]:
        """List all registered agent names."""
        return list(cls._agents.keys())
