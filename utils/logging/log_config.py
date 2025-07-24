"""
Logging configuration module for F.A.D.E

This module provides a structured logging setup with different log files for each agent
and component of the system. It ensures that logs are organized and easy to find.
"""

import os
import sys
import logging
import logging.handlers
from typing import Dict, Optional, Union
from datetime import datetime

# Define agent types
AGENT_TYPES = [
    "target_selector",
    "structure_predictor",
    "molecule_generator",
    "evaluator",
    "docking",
    "refiner"
]

# Global storage for loggers to prevent duplicate configuration
_loggers: Dict[str, logging.Logger] = {}

def _ensure_logs_directory(logs_dir: str) -> None:
    """
    Ensure that the logs directory structure exists.
    
    Args:
        logs_dir: Base logs directory
    """
    # Create main logs directory
    os.makedirs(logs_dir, exist_ok=True)
    
    # Create directories for each agent type
    for agent_type in AGENT_TYPES:
        agent_dir = os.path.join(logs_dir, agent_type)
        os.makedirs(agent_dir, exist_ok=True)
    
    # Create directory for main logs
    os.makedirs(os.path.join(logs_dir, "main"), exist_ok=True)
    
    # Create directory for utility logs
    os.makedirs(os.path.join(logs_dir, "utils"), exist_ok=True)

def _get_run_id() -> str:
    """
    Get a unique run ID based on the current timestamp.
    
    Returns:
        String timestamp in the format YYYYMMDD_HHMMSS
    """
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def _get_log_file_path(logger_name: str, logs_dir: str, run_id: str) -> str:
    """
    Get the appropriate log file path based on the logger name.
    
    Args:
        logger_name: Name of the logger
        logs_dir: Base logs directory
        run_id: Unique run identifier
        
    Returns:
        Path to the log file
    """
    # Split the logger name to determine the appropriate directory
    parts = logger_name.split('.')
    
    if len(parts) >= 3 and parts[0] == "fade" and parts[1] == "agent":
        # This is an agent logger (e.g., fade.agent.target_selector)
        agent_type = parts[2]
        if agent_type in AGENT_TYPES:
            return os.path.join(logs_dir, agent_type, f"{agent_type}_{run_id}.log")
        else:
            # For custom agents not in the predefined list
            return os.path.join(logs_dir, "agents", f"{agent_type}_{run_id}.log")
    elif len(parts) >= 2 and parts[0] == "fade" and parts[1] == "utils":
        # This is a utility logger
        return os.path.join(logs_dir, "utils", f"utils_{run_id}.log")
    else:
        # Default to main logs
        return os.path.join(logs_dir, "main", f"fade_{run_id}.log")

def setup_logging(
    log_level: Union[str, int] = "INFO", 
    logs_dir: str = "logs",
    run_id: Optional[str] = None,
    console_output: bool = True
) -> None:
    """
    Set up logging for the application with structured organization.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        logs_dir: Directory where logs should be stored
        run_id: Optional unique run identifier, will be generated if not provided
        console_output: Whether to also output logs to the console
    """
    # Convert string log level to numeric value if needed
    if isinstance(log_level, str):
        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {log_level}")
        log_level = numeric_level
    
    # Generate run ID if not provided
    if run_id is None:
        run_id = _get_run_id()
    
    # Ensure logs directory structure exists
    _ensure_logs_directory(logs_dir)
    
    # Configure the root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Remove existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Add console handler if requested
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)
    
    # Add file handler for the main log
    main_log_path = os.path.join(logs_dir, "main", f"fade_{run_id}.log")
    file_handler = logging.FileHandler(main_log_path)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)
    
    # Log setup information
    root_logger.info(f"Logging initialized with level {logging.getLevelName(log_level)}")
    root_logger.info(f"Run ID: {run_id}")
    root_logger.info(f"Main log file: {main_log_path}")
    
    # Store run ID for future reference
    global _current_run_id, _current_logs_dir
    _current_run_id = run_id
    _current_logs_dir = logs_dir

# Store current run ID and logs directory
_current_run_id: Optional[str] = None
_current_logs_dir: str = "logs"

def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name, configured for the current run.
    
    Args:
        name: Name of the logger
        
    Returns:
        Configured logger instance
    """
    # Check if this logger already exists
    if name in _loggers:
        return _loggers[name]
    
    # Get the logger
    logger = logging.getLogger(name)
    
    # If we have a run ID, add a specific file handler for this logger
    if _current_run_id is not None:
        log_file = _get_log_file_path(name, _current_logs_dir, _current_run_id)
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Add file handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
        # Log the creation of this logger
        logger.info(f"Logger initialized with specific log file: {log_file}")
    
    # Store the logger for future reference
    _loggers[name] = logger
    
    return logger
