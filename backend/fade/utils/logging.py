"""
Logging configuration for F.A.D.E backend
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional

from rich.logging import RichHandler
from rich.console import Console

from fade.config import config


def setup_logging(
    name: Optional[str] = None,
    level: Optional[str] = None,
    log_file: Optional[Path] = None,
    use_rich: bool = True
) -> logging.Logger:
    """
    Set up logging configuration.
    
    Args:
        name: Logger name (defaults to "fade")
        level: Logging level (defaults to config.LOG_LEVEL)
        log_file: Optional file to log to
        use_rich: Whether to use Rich for pretty console output
        
    Returns:
        Configured logger instance
    """
    logger_name = name or "fade"
    log_level = level or config.LOG_LEVEL
    
    # Get or create logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Remove existing handlers to avoid duplicates
    logger.handlers.clear()
    
    # Console handler with Rich formatting
    if use_rich:
        console = Console(stderr=True)
        console_handler = RichHandler(
            console=console,
            show_time=True,
            show_path=False,
            markup=True,
            rich_tracebacks=True,
            tracebacks_show_locals=True
        )
    else:
        console_handler = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(config.LOG_FORMAT)
        console_handler.setFormatter(formatter)
    
    console_handler.setLevel(getattr(logging, log_level.upper()))
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        # Create log directory if it doesn't exist
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_formatter = logging.Formatter(config.LOG_FORMAT)
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(getattr(logging, log_level.upper()))
        logger.addHandler(file_handler)
        
    # Also create a default log file in the logs directory
    if not log_file and logger_name == "fade":
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_log_file = config.LOGS_DIR / f"fade_{timestamp}.log"
        
        file_handler = logging.FileHandler(default_log_file, encoding='utf-8')
        file_formatter = logging.Formatter(config.LOG_FORMAT)
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(getattr(logging, log_level.upper()))
        logger.addHandler(file_handler)
        
        logger.info(f"Logging to {default_log_file}")
    
    return logger


# Create default logger
default_logger = setup_logging()


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance with the specified name.
    
    Args:
        name: Logger name (will be prefixed with "fade.")
        
    Returns:
        Logger instance
    """
    if not name.startswith("fade."):
        name = f"fade.{name}"
    return logging.getLogger(name)
