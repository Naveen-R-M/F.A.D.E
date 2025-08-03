### Running on HPC Clusters

For better compatibility with HPC environments, you can use the provided shell script:

```bash
chmod +x run_fade.sh
./run_fade.sh "Find molecules targeting KRAS G12D with good BBB permeability"
```

To skip the structure prediction step:

```bash
./run_fade.sh "Find molecules targeting KRAS G12D with good BBB permeability" --skip-structure-prediction
```

This script uses `nohup` to ensure the job continues running even if your SSH session disconnects.# F.A.D.E: Fully Agentic Drug Engine

![F.A.D.E Logo](https://via.placeholder.com/800x200?text=F.A.D.E)

## Overview

F.A.D.E is an agentic AI pipeline for discovering drug-like molecules targeting specific proteins. The system transforms natural language queries into a computational drug discovery workflow, leveraging state-of-the-art tools for protein structure prediction, molecule generation, property evaluation, and docking simulations.

## Features

- **Natural Language Interface**: Describe your drug discovery goals in plain English
- **Multi-Agent Architecture**: Specialized agents handle different aspects of the discovery process
- **HPC Integration**: Designed for high-performance computing environments
- **State-of-the-Art Tools**: Utilizes AlphaFold3, RDKit, and Schrodinger Glide
- **LLM-Guided Molecule Generation**: Generates novel molecules guided by large language models and RDKit
- **Iterative Refinement**: Molecules are refined through multiple cycles of evaluation and improvement

## System Requirements

- **HPC Environment**: Access to a cluster with SLURM scheduler
- **GPU Resources**: For protein structure prediction
- **Software Dependencies**:
  - AlphaFold3 container
  - Schrodinger Suite with Glide
  - Python 3.9+
  - RDKit
  - LangGraph

## Installation

1. Clone this repository to your HPC environment:
```bash
git clone https://github.com/Naveen-R-M/F.A.D.E.git
cd F.A.D.E.
```

2. Create a conda environment:

   **For HPC (SLURM) environments:**
   ```bash
   module load miniconda3/24.11.1
   conda create -p $SCRATCH/conda-envs/fade python=3.9 -y
   conda activate $SCRATCH/conda-envs/fade
   pip install -r requirements.txt
   ```

   **For local development (general Conda):**
   ```bash
   conda create -n fade python=3.9 -y
   conda activate fade
   pip install -r requirements.txt
   ```

3. Configure your environment:
```bash
cp .env.example .env
# Edit .env with your API keys and paths
```

4. Verify access to required containers and modules:
```bash
module load schrodinger/2024-4
singularity exec --nv /shared/container_repository/AlphaFold/alphafold3.sif python -c "import sys; print(sys.version)"
```

## Usage

Run the main interface (defaults to background execution):

```bash
python main.py --query "Find molecules targeting KRAS G12D with good BBB permeability"
```

This will automatically run the process in the background and provide you with monitoring options. The results will be saved to a timestamped directory.

To force running in the foreground (not recommended for long-running jobs):

```bash
python main.py --query "Find molecules targeting KRAS G12D with good BBB permeability" --foreground
```

For batch processing:

```bash
python main.py --batch-file queries.txt --output-dir results/
```

### Monitoring Background Jobs

To monitor all running F.A.D.E jobs:

```bash
python monitor_jobs.py
```

To monitor a specific job:

```bash
python monitor_jobs.py results_20250728_123456
```

### Running on HPC Clusters

For better compatibility with HPC environments, you can use the provided shell script:

```bash
chmod +x run_fade.sh
./run_fade.sh "Find molecules targeting KRAS G12D with good BBB permeability"
```

To skip the structure prediction step:

```bash
./run_fade.sh "Find molecules targeting KRAS G12D with good BBB permeability" --skip-structure-prediction
```

This script uses `nohup` to ensure the job continues running even if your SSH session disconnects.

## Project Structure

```
fade/
├── agents/               # Specialized agents for different tasks
├── configs/              # Configuration files (e.g., default parameters)
├── data/                 # Input/output data storage
│   ├── inputs/           # User inputs and configurations
│   ├── outputs/          # Generated results
│   └── examples/         # Example configurations
├── logs/                 # Log files
├── main.py               # Main entry point
├── models/               # ML models for property prediction (future use)
├── requirements.txt      # Python dependencies
├── tests/                # Unit and integration tests
├── .env                  # Environment variables (not in version control)
├── utils/                # Utility functions and helpers
└── workflows/            # Pipeline orchestration and job templates
```

## Workflow

See [Workflow.md](Workflow.md) for a detailed description of the pipeline process.

## Configuration

The system is configured through:

1. **Environment Variables**: API keys, paths to external tools (.env file)
2. **Configuration Files**: Default parameters and settings (configs directory)
3. **Command-Line Arguments**: Runtime options and queries

## Examples

Example query:
```
"I'm looking for potential drug candidates that could target the KRAS G12D mutant protein, which is implicated in pancreatic cancer. The molecules should be able to cross the blood-brain barrier, have low toxicity, and follow Lipinski's rule of five for drug-likeness."
```

See the `data/examples/` directory for more sample configurations.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use F.A.D.E in your research, please cite:

```
@software{fade2025,
  author = {Your Name},
  title = {F.A.D.E: Fully Agentic Drug Engine},
  year = {2025},
  url = {https://github.com/Naveen-R-M/F.A.D.E.git}
}
```

## Acknowledgments

- This project was developed on the Northeastern University HPC cluster
- Special thanks to the contributors of AlphaFold, RDKit, and Schrodinger
