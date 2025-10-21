"""
Configuration management for F.A.D.E backend
"""

import os
from pathlib import Path
from typing import Optional, Dict, Any
from dotenv import load_dotenv

# Load environment variables
load_dotenv()


class Config:
    """Central configuration for the drug discovery pipeline."""
    
    # Project paths
    PROJECT_ROOT = Path(__file__).parent.parent
    DATA_DIR = PROJECT_ROOT / "data"
    OUTPUT_DIR = PROJECT_ROOT / "output"
    CHECKPOINT_DIR = PROJECT_ROOT / "checkpoints"
    LOGS_DIR = PROJECT_ROOT / "logs"
    
    # Ensure directories exist
    for dir_path in [DATA_DIR, OUTPUT_DIR, CHECKPOINT_DIR, LOGS_DIR]:
        dir_path.mkdir(parents=True, exist_ok=True)
    
    # API Keys and Endpoints
    # LLM APIs
    OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
    ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY")
    GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
    AI_MODEL = os.getenv("AI_MODEL", "gpt-4")
    AI_API_BASE = os.getenv("AI_API_BASE")  # For local models like DeepSeek
    
    # Protein Database APIs
    UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb"
    RCSB_API_URL = "https://data.rcsb.org/graphql"
    ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api"
    CHEMBL_API_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    # Boltz2 Configuration
    BOLTZ2_SERVER_URL = os.getenv("BOLTZ2_SERVER_URL", "http://localhost:8080")
    BOLTZ2_API_KEY = os.getenv("BOLTZ2_API_KEY")
    BOLTZ2_MAX_SEQUENCE_LENGTH = int(os.getenv("BOLTZ2_MAX_SEQUENCE_LENGTH", "2048"))
    
    # SILCS Configuration
    SILCS_EXECUTABLE = os.getenv("SILCS_EXECUTABLE", "silcs")
    SILCS_LICENSE = os.getenv("SILCS_LICENSE")
    
    # DiffSBDD Configuration
    DIFFSBDD_MODEL_PATH = os.getenv("DIFFSBDD_MODEL_PATH", str(DATA_DIR / "models" / "diffsbdd"))
    DIFFSBDD_DEVICE = os.getenv("DIFFSBDD_DEVICE", "cuda" if os.getenv("CUDA_VISIBLE_DEVICES") else "cpu")
    
    # Pipeline Configuration
    MAX_MOLECULES_TO_GENERATE = int(os.getenv("MAX_MOLECULES_TO_GENERATE", "1000"))
    MAX_MOLECULES_TO_SCREEN = int(os.getenv("MAX_MOLECULES_TO_SCREEN", "100"))
    MAX_POCKETS_TO_ANALYZE = int(os.getenv("MAX_POCKETS_TO_ANALYZE", "5"))
    
    # Filtering Thresholds
    MAX_MOLECULAR_WEIGHT = float(os.getenv("MAX_MOLECULAR_WEIGHT", "500"))
    MIN_LOGP = float(os.getenv("MIN_LOGP", "-0.4"))
    MAX_LOGP = float(os.getenv("MAX_LOGP", "5.6"))
    MAX_ROTATABLE_BONDS = int(os.getenv("MAX_ROTATABLE_BONDS", "10"))
    MAX_HBD = int(os.getenv("MAX_HBD", "5"))  # Hydrogen Bond Donors
    MAX_HBA = int(os.getenv("MAX_HBA", "10"))  # Hydrogen Bond Acceptors
    MIN_QED_SCORE = float(os.getenv("MIN_QED_SCORE", "0.3"))
    MAX_SA_SCORE = float(os.getenv("MAX_SA_SCORE", "6.0"))  # Synthetic Accessibility
    
    # Performance Configuration
    ENABLE_CACHING = os.getenv("ENABLE_CACHING", "true").lower() == "true"
    CACHE_TTL = int(os.getenv("CACHE_TTL", "3600"))  # seconds
    MAX_PARALLEL_JOBS = int(os.getenv("MAX_PARALLEL_JOBS", "4"))
    REQUEST_TIMEOUT = int(os.getenv("REQUEST_TIMEOUT", "300"))  # seconds
    
    # Logging Configuration
    LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO")
    LOG_FORMAT = os.getenv("LOG_FORMAT", "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # HPC Configuration (if using cluster)
    USE_HPC = os.getenv("USE_HPC", "false").lower() == "true"
    SLURM_PARTITION = os.getenv("SLURM_PARTITION", "default")
    SLURM_ACCOUNT = os.getenv("SLURM_ACCOUNT")
    
    @classmethod
    def get_api_headers(cls, api_name: str) -> Dict[str, str]:
        """Get headers for API requests."""
        headers = {
            "Content-Type": "application/json",
            "User-Agent": "F.A.D.E/2.0 (Drug Discovery Pipeline)"
        }
        
        if api_name == "boltz2" and cls.BOLTZ2_API_KEY:
            headers["Authorization"] = f"Bearer {cls.BOLTZ2_API_KEY}"
            
        return headers
    
    @classmethod
    def get_llm_config(cls) -> Dict[str, Any]:
        """Get configuration for LLM initialization."""
        model = cls.AI_MODEL.lower()
        
        # Handle local models like DeepSeek
        if cls.AI_API_BASE:
            return {
                "model": cls.AI_MODEL,
                "api_key": cls.OPENAI_API_KEY or "dummy",  # Some local models need a dummy key
                "base_url": cls.AI_API_BASE,
                "temperature": 0.7,
                "max_tokens": 4000
            }
        
        # Handle cloud models
        if "gpt" in model:
            return {
                "model": cls.AI_MODEL,
                "api_key": cls.OPENAI_API_KEY,
                "temperature": 0.7,
                "max_tokens": 4000
            }
        elif "claude" in model:
            return {
                "model": cls.AI_MODEL,
                "api_key": cls.ANTHROPIC_API_KEY,
                "temperature": 0.7,
                "max_tokens": 4000
            }
        elif "gemini" in model:
            return {
                "model": cls.AI_MODEL,
                "api_key": cls.GOOGLE_API_KEY,
                "temperature": 0.7,
                "max_tokens": 4000
            }
        else:
            # Default to OpenAI-compatible endpoint
            return {
                "model": cls.AI_MODEL,
                "api_key": cls.OPENAI_API_KEY or "dummy",
                "base_url": cls.AI_API_BASE,
                "temperature": 0.7,
                "max_tokens": 4000
            }
    
    @classmethod
    def validate(cls) -> None:
        """Validate that required configuration is present."""
        errors = []
        
        # Check required API keys or base URL for LLM
        if not cls.AI_API_BASE and not any([cls.OPENAI_API_KEY, cls.ANTHROPIC_API_KEY, cls.GOOGLE_API_KEY]):
            errors.append("Either AI_API_BASE or at least one LLM API key must be set")
            
        # Check Boltz2 configuration (can be optional for now)
        # if not cls.BOLTZ2_SERVER_URL:
        #     errors.append("BOLTZ2_SERVER_URL must be set")
            
        if errors:
            raise ValueError(f"Configuration errors: {'; '.join(errors)}")


# Create a singleton config instance
config = Config()
