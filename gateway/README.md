# F.A.D.E Gateway

🧬 **Intelligent gateway for the F.A.D.E drug discovery platform**

## Overview

The F.A.D.E Gateway is an AI-powered interface that:
- **Classifies queries** using Gemini AI (CHAT vs JOB)
- **Provides direct answers** for general questions
- **Routes drug discovery requests** to the backend pipeline
- **Manages job execution** with real-time status tracking

## Files

```
gateway/
├── app.py                 # Main FastAPI application
├── requirements.txt       # Python dependencies
├── sample_queries.txt     # Training examples for classification
├── project_context.txt    # Platform knowledge base
├── .env.example          # Environment template
└── README.md             # This file
```

## Quick Start

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Set up environment:**
   ```bash
   cp .env.example .env
   # Add your GEMINI_API_KEY to .env
   ```

3. **Run the gateway:**
   ```bash
   python app.py
   ```

4. **Update frontend:**
   ```bash
   # In frontend/.env.local
   NEXT_PUBLIC_BACKEND_URL=http://localhost:8000
   ```

## API Endpoints

- `POST /jobs` - Submit queries (auto-classified as chat or job)
- `GET /jobs/{id}/status` - Check job status
- `GET /jobs` - List all jobs  
- `GET /health` - System health check
- `GET /docs` - Interactive API documentation

## Examples

**Chat Query:**
```json
{"query": "What is F.A.D.E?"}
→ Returns direct explanation
```

**Job Query:**
```json
{"query": "Find EGFR inhibitors for lung cancer"}
→ Submits to pipeline, returns job ID
```

## Features

✅ **Gemini AI Classification** - Intelligent query routing  
✅ **RAG Enhancement** - Uses sample queries and context  
✅ **Backend Integration** - Connects to F.A.D.E pipeline  
✅ **Job Management** - Real-time status and monitoring  
✅ **Fallback Support** - Works without Gemini API  
✅ **CORS Enabled** - Ready for frontend integration  

## Configuration

The gateway automatically:
- Loads 120+ sample queries for better classification
- Uses project context for informed chat responses
- Falls back to keyword matching if Gemini is unavailable
- Connects to `../backend/main.py` for pipeline execution

## Monitoring

- **Gateway:** http://localhost:8000
- **Health:** http://localhost:8000/health  
- **Docs:** http://localhost:8000/docs
- **Jobs:** http://localhost:8000/jobs

Ready to power your drug discovery conversations! 🚀