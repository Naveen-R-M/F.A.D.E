import os
from dotenv import load_dotenv

load_dotenv()

import os
from dotenv import load_dotenv

load_dotenv()

GEMINI_API_KEY = os.getenv("GEMINI_API_KEY") or ""
GEMINI_MODEL = os.getenv("GEMINI_MODEL") or "gemini-2.5-pro"

OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY", "")
OPENROUTER_BASE_URL = os.getenv("OPENROUTER_BASE_URL", "https://openrouter.ai/api/v1")
OPENROUTER_MODEL = os.getenv("OPENROUTER_MODEL", "deepseek/deepseek-chat-v3.1:free")

OPENROUTER_REFERRER = os.getenv("OPENROUTER_REFERRER", "")
OPENROUTER_TITLE = os.getenv("OPENROUTER_TITLE", "")

# CORS
CORS_ORIGINS = os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")

# Memory caps
MAX_TURNS_PER_USER = int(os.getenv("MAX_TURNS_PER_USER", "12"))
MAX_CHARS_PER_USER = int(os.getenv("MAX_CHARS_PER_USER", "6000"))

# CORS (adjust for your UI)
CORS_ORIGINS = os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")

# --- HPC / SSH ---
HPC_HOST = os.getenv("HPC_HOST", "")
HPC_PORT = int(os.getenv("HPC_PORT", "22"))
HPC_USER = os.getenv("HPC_USER", "")

SSH_PRIVATE_KEY_PATH = os.getenv("SSH_PRIVATE_KEY_PATH") or None
SSH_KEY_PASSPHRASE = os.getenv("SSH_KEY_PASSPHRASE") or None
HPC_PASSWORD = os.getenv("HPC_PASSWORD") or None  # only if not using keys

# Optional bastion/jump host
BASTION_HOST = os.getenv("BASTION_HOST") or None
BASTION_USER = os.getenv("BASTION_USER") or None
BASTION_PORT = int(os.getenv("BASTION_PORT", "22"))

# Where the script must run from:
HPC_SCRATCH_DIR = os.getenv("HPC_SCRATCH_DIR", "/scratch/rajagopalmohanraj.n")
HPC_SCRIPT_PATH = os.getenv("HPC_SCRIPT_PATH", os.path.join(HPC_SCRATCH_DIR, "F.A.D.E/run_fade_nextflow.sh"))

# Connection behavior
STRICT_HOST_KEY_CHECKING = os.getenv("STRICT_HOST_KEY_CHECKING", "true").lower() == "true"
KNOWN_HOSTS_PATH = os.getenv("KNOWN_HOSTS_PATH") or None
HPC_CONNECT_TIMEOUT = float(os.getenv("HPC_CONNECT_TIMEOUT", "15"))
HPC_KEEPALIVE_SECS = int(os.getenv("HPC_KEEPALIVE_SECS", "30"))

# Command timeout (seconds) — script should return quickly (e.g., after sbatch submission)
HPC_SUBMIT_TIMEOUT = float(os.getenv("HPC_SUBMIT_TIMEOUT", "120"))
