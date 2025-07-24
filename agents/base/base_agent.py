"""
Base agent class for F.A.D.E framework.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional


class BaseAgent(ABC):
    """
    Abstract base class for all agents in the F.A.D.E framework.
    
    This class defines the common interface and functionality that all agents
    should implement, ensuring consistent behavior across the framework.
    """
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None) -> None:
        """
        Initialize the base agent.
        
        Args:
            name: A unique identifier for the agent.
            config: Optional configuration parameters for the agent.
        """
        self.name = name
        self.config = config or {}
        self.logger = self._setup_logger()
        
    def _setup_logger(self):
        """Set up a logger for the agent."""
        import logging
        logger = logging.getLogger(f"fade.agent.{self.name}")
        return logger
    
    @abstractmethod
    def process(self, input_data: Any) -> Any:
        """
        Process the input data and return the results.
        
        This is the main method that all agents must implement to perform their
        specific tasks.
        
        Args:
            input_data: The data to be processed by the agent.
            
        Returns:
            The processed results.
        """
        pass
    
    def validate_input(self, input_data: Any) -> bool:
        """
        Validate the input data before processing.
        
        Args:
            input_data: The data to be validated.
            
        Returns:
            True if the input is valid, False otherwise.
        """
        # Default implementation always returns True
        # Subclasses should override this method to implement specific validation
        return True
    
    def validate_output(self, output_data: Any) -> bool:
        """
        Validate the output data after processing.
        
        Args:
            output_data: The data to be validated.
            
        Returns:
            True if the output is valid, False otherwise.
        """
        # Default implementation always returns True
        # Subclasses should override this method to implement specific validation
        return True
