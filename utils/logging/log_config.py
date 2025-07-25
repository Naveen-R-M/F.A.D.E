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

# Define log levels
LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL
}

# Define log format
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

# Set up console formatter
_console_formatter = logging.Formatter(LOG_FORMAT, DATE_FORMAT)

# Set up file formatter with more details
_file_formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s",
    DATE_FORMAT
)

# Global variables for current run ID and logs directory
_current_run_id = None
_current_logs_dir = None


def _configure_root_logger(log_level: str = "INFO") -> None:
    """
    Configure the root logger.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    """
    # Get the root logger
    root_logger = logging.getLogger()
    
    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Set the log level
    root_logger.setLevel(LOG_LEVELS.get(log_level, logging.INFO))
    
    # Create a console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(_console_formatter)
    root_logger.addHandler(console_handler)


def _configure_fade_logger(
    logs_dir: str,
    run_id: str,
    log_level: str = "INFO"
) -> logging.Logger:
    """
    Configure the main F.A.D.E logger.
    
    Args:
        logs_dir: Directory where logs will be stored.
        run_id: Unique identifier for the current run.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        
    Returns:
        The configured logger.
    """
    # Create the main F.A.D.E logger
    logger = logging.getLogger("fade")
    
    # Set the log level
    logger.setLevel(LOG_LEVELS.get(log_level, logging.INFO))
    
    # Create a file handler for the main log file
    main_log_file = os.path.join(logs_dir, f"{run_id}.log")
    file_handler = logging.FileHandler(main_log_file)
    file_handler.setFormatter(_file_formatter)
    logger.addHandler(file_handler)
    
    # Create a file handler for errors only
    error_log_file = os.path.join(logs_dir, f"{run_id}_errors.log")
    error_handler = logging.FileHandler(error_log_file)
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(_file_formatter)
    logger.addHandler(error_handler)
    
    return logger


def _get_log_file_path(name: str, logs_dir: str, run_id: str) -> str:
    """
    Get the path for a specific log file.
    
    Args:
        name: Name of the logger.
        logs_dir: Directory where logs will be stored.
        run_id: Unique identifier for the current run.
        
    Returns:
        Path to the log file.
    """
    # Extract the component name from the logger name
    parts = name.split(".")
    
    if len(parts) > 1:
        # For agents and components (e.g., fade.agent.target_selector)
        component_type = parts[1]
        component_name = ".".join(parts[2:]) if len(parts) > 2 else "main"
        
        if component_type == "agent" and component_name in AGENT_TYPES:
            # Create a subdirectory for each agent
            agent_dir = os.path.join(logs_dir, "agents", component_name)
            os.makedirs(agent_dir, exist_ok=True)
            return os.path.join(agent_dir, f"{run_id}.log")
        else:
            # Create a subdirectory for each component type
            component_dir = os.path.join(logs_dir, component_type)
            os.makedirs(component_dir, exist_ok=True)
            return os.path.join(component_dir, f"{component_name}_{run_id}.log")
    else:
        # For other loggers
        return os.path.join(logs_dir, f"{name}_{run_id}.log")


def set_run_id(run_id: str) -> None:
    """
    Set the current run ID for logging.
    
    Args:
        run_id: Unique identifier for the current run.
    """
    global _current_run_id, _current_logs_dir
    _current_run_id = run_id


def setup_logging(log_level: str = "INFO", logs_dir: str = "logs") -> None:
    """
    Set up logging for the application.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        logs_dir: Directory where logs will be stored.
    """
    # Create a unique run ID based on the current time
    run_id = datetime.now().strftime("%Y%m%d%H%M%S")
    
    # Create the logs directory if it doesn't exist
    os.makedirs(logs_dir, exist_ok=True)
    
    # Configure the root logger
    _configure_root_logger(log_level)
    
    # Configure the main F.A.D.E logger
    _configure_fade_logger(logs_dir, run_id, log_level)
    
    # Set the current run ID and logs directory
    global _current_run_id, _current_logs_dir
    _current_run_id = run_id
    _current_logs_dir = logs_dir
    
    # Log the start of the application
    logger = logging.getLogger("fade")
    logger.info(f"Logging initialized with run ID: {run_id}")
    logger.info(f"Logs directory: {os.path.abspath(logs_dir)}")


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.
    
    Args:
        name: Name of the logger.
        
    Returns:
        The logger.
    """
    # Get the logger
    logger = logging.getLogger(name)
    
    # Add a file handler if run ID is set
    if _current_run_id is not None:
        log_file = _get_log_file_path(name, _current_logs_dir, _current_run_id)
        
        # Check if the logger already has a file handler
        has_file_handler = any(
            isinstance(h, logging.FileHandler) for h in logger.handlers
        )
        
        if not has_file_handler:
            # Create a file handler
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(_file_formatter)
            logger.addHandler(file_handler)
    
    return logger
