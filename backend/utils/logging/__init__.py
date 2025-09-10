"""
Logging configuration for F.A.D.E
"""

# Import package modules to make them available through the package
from .log_config import setup_logging, get_logger

__all__ = ["setup_logging", "get_logger"]
