# F.A.D.E Project Structure

## 📁 Directory Organization

```
F.A.D.E/
├── 📂 agents/                          # AI Agents (Core Components)
│   ├── base/                           # Base agent classes
│   ├── target_selector/                # Target identification agent
│   ├── structure_predictor/            # Protein structure prediction agent
│   └── molecule_generator/             # Molecule generation agent
│
├── 📂 nextflow/                        # Nextflow Workflow Implementation
│   ├── main.nf                         # Main workflow definition
│   ├── nextflow.config                 # Workflow configuration
│   ├── modules/                        # Individual workflow modules
│   └── bin/                            # Python wrapper scripts
│
├── 📂 utils/                           # Utility Libraries
│   ├── gemini_client.py               # Gemini API client
│   ├── uniprot_client.py              # UniProt API client
│   ├── alphafold_client.py            # AlphaFold3 integration
│   ├── slurm_client.py                # SLURM job management
│   ├── config_generator.py            # Configuration file generator
│   └── logging/                       # Logging utilities
│
├── 📂 data/                           # Data Storage
│   ├── inputs/                        # Input configurations and sequences
│   ├── outputs/                       # Generated results (runtime)
│   └── examples/                      # Example data and configurations
│
├── 📂 scripts/                        # Execution and Utility Scripts
│   ├── execution/                     # Main execution scripts
│   │   ├── run_nextflow_integrated.sh # Primary Nextflow execution
│   │   ├── run_nextflow_with_api.sh   # API-aware execution
│   │   └── run_fade.sh                # Legacy Python execution
│   ├── setup/                         # Environment setup scripts
│   │   ├── setup_nextflow_env.sh      # Nextflow environment setup
│   │   └── fix_conda_env.sh           # Conda environment fixes
│   └── maintenance/                   # Maintenance utilities
│       ├── check_disk_usage.py        # Disk usage monitoring
│       ├── clear_logs.py              # Log cleanup
│       ├── clear_outputs.py           # Output cleanup
│       └── monitor_all_jobs.py        # Job monitoring
│
├── 📂 docs/                           # Documentation
│   ├── integration/                   # Nextflow integration docs
│   ├── architecture/                  # System architecture docs
│   └── api/                          # API documentation
│
├── 📂 tests/                          # Testing Framework
│   ├── unit/                         # Unit tests for individual components
│   └── integration/                  # Integration tests for workflows
│
├── 📂 examples/                       # Usage Examples
│   └── queries/                      # Example drug discovery queries
│
├── 📂 logs/                          # Runtime Logs
│   └── (generated during execution)
│
├── 📂 results/                       # Results (Runtime)
│   └── (generated during execution)
│
├── 📄 main.py                        # Primary Python execution interface
├── 📄 README.md                      # Main project documentation
├── 📄 requirements.txt               # Python dependencies
├── 📄 .env                          # Environment variables (API keys)
├── 📄 .env.example                  # Environment template
├── 📄 TodoList.md                   # Development roadmap
└── 📄 Workflow.md                   # Workflow documentation
```

## 🚀 Quick Start Commands

### Nextflow Execution (Recommended)
```bash
# Setup environment (one-time)
./scripts/setup/setup_nextflow_env.sh
./scripts/setup/fix_conda_env.sh

# Run drug discovery workflow
./scripts/execution/run_nextflow_with_api.sh "Find molecules targeting KRAS G12D"
```

### Python Execution (Legacy)
```bash
python main.py --query "Your drug discovery query"
```

### Maintenance
```bash
# Monitor jobs
./scripts/maintenance/monitor_all_jobs.py

# Clean up old files
./scripts/maintenance/clear_logs.py
./scripts/maintenance/clear_outputs.py
```

## 📚 Documentation Index

- **Main README**: `README.md` - Primary project documentation
- **Workflow Guide**: `Workflow.md` - Detailed workflow explanation
- **Integration Guide**: `docs/integration/` - Nextflow integration documentation
- **API Documentation**: `docs/api/` - API and agent documentation
- **Architecture**: `docs/architecture/` - System architecture documentation

## 🎯 Key Files for Users

- **Primary Execution**: `./scripts/execution/run_nextflow_with_api.sh`
- **Environment Setup**: `./scripts/setup/setup_nextflow_env.sh`
- **Configuration**: `.env` (set your GEMINI_API_KEY here)
- **Examples**: `examples/queries/` 
- **Results**: `results/` (generated during execution)

This organization provides a clean, professional structure that separates:
- **Core components** (agents, utils)
- **Execution interfaces** (nextflow, scripts)
- **Documentation** (organized by topic)
- **Data and results** (isolated from code)
- **Testing and examples** (clear separation)
