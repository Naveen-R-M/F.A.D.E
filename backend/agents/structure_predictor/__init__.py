"""
Structure Predictor Package - Enhanced with RCSB Integration
"""

from agents.structure_predictor.structure_predictor import StructurePredictor

# NEW: Export RCSB client for direct use if needed
try:
    from agents.structure_predictor.rcsb_client import RCSBClient
    __all__ = ["StructurePredictor", "RCSBClient"]
except ImportError:
    # RCSB client not available (missing dependencies)
    __all__ = ["StructurePredictor"]