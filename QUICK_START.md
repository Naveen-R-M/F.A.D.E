# F.A.D.E Quick Start Guide

## 🚀 Getting Started (30 seconds)

### 1. Setup Environment (One-time)
```bash
./setup_env.sh              # Setup Nextflow environment
./scripts/setup/fix_conda_env.sh  # Fix conda dependencies
```

### 2. Configure API Key
```bash
# Edit .env file and add your Gemini API key
echo "GEMINI_API_KEY=your_api_key_here" >> .env
```

### 3. Run Drug Discovery
```bash
# Use the main execution script
./run_fade_nextflow.sh "Find molecules targeting KRAS G12D with good BBB permeability"
```

## 📖 Example Queries
```bash
# EGFR targeting for lung cancer
./run_fade_nextflow.sh "Design EGFR inhibitors for lung cancer with oral bioavailability"

# BRAF targeting for melanoma
./run_fade_nextflow.sh "Find BRAF V600E inhibitors that cross blood-brain barrier"

# KRAS targeting for pancreatic cancer
./run_fade_nextflow.sh "Target KRAS G12D mutant with selective small molecules"
```

## 🔧 Maintenance
```bash
# Monitor running jobs
./monitor_jobs.py

# Clean up old files
./scripts/maintenance/clear_logs.py --days 7
./scripts/maintenance/clear_outputs.py --keep-best
```

## 📁 Project Structure
See `PROJECT_STRUCTURE.md` for detailed organization overview.

---
**F.A.D.E: Fully Agentic Drug Engine - Ready for AI-powered drug discovery!**
