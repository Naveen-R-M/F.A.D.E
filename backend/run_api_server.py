#!/usr/bin/env python3
"""
Run the F.A.D.E Unified API Server with visible logging.
"""

import os
import sys
import logging
from pathlib import Path

# Add backend to path
backend_dir = Path(__file__).parent
sys.path.insert(0, str(backend_dir))

# Set environment variables if not set
os.environ.setdefault("API_PORT", "8000")
os.environ.setdefault("API_HOST", "0.0.0.0")

# Configure Python logging to show all output
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

# Force stdout to be unbuffered so prints appear immediately
sys.stdout.reconfigure(line_buffering=True)

if __name__ == "__main__":
    import uvicorn
    
    port = int(os.getenv("API_PORT", "8000"))
    host = os.getenv("API_HOST", "0.0.0.0")
    
    print("=" * 60)
    print("F.A.D.E Unified Backend")
    print("=" * 60)
    print(f"Starting server on {host}:{port}")
    print(f"Documentation: http://localhost:{port}/docs")
    print(f"Health check: http://localhost:{port}/api/health")
    print("-" * 60)
    
    # Configure uvicorn logging to be more verbose
    log_config = uvicorn.config.LOGGING_CONFIG
    log_config["formatters"]["default"]["fmt"] = "%(asctime)s - %(levelprefix)s %(message)s"
    log_config["formatters"]["access"]["fmt"] = '%(asctime)s - %(levelprefix)s %(client_addr)s - "%(request_line)s" %(status_code)s'
    
    # Use string import for proper reload functionality
    uvicorn.run(
        "api_server:app",
        host=host,
        port=port,
        reload=True,
        log_level="info",
        log_config=log_config,
        access_log=True,
        use_colors=True
    )
