# F.A.D.E: Fully Agentic Drug Engine

![F.A.D.E Logo](https://via.placeholder.com/800x200?text=F.A.D.E)

## Overview

F.A.D.E is an agentic AI pipeline for discovering drug-like molecules targeting specific proteins. The system transforms natural language queries into a computational drug discovery workflow, leveraging state-of-the-art tools for protein structure prediction, molecule generation, property evaluation, and docking simulations.

## Features

- **Natural Language Interface**: Describe your drug discovery goals in plain English
- **Multi-Agent Architecture**: Specialized agents handle different aspects of the discovery process
- **HPC Integration**: Designed for high-performance computing environments
- **State-of-the-Art Tools**: Utilizes AlphaFold3, RDKit, and Schrodinger Glide
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
```bash
module load miniconda3/24.11.1
conda create -p $SCRATCH/conda-envs/fade python=3.9 -y
conda activate $SCRATCH/conda-envs/fade
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

Run the main interface:

```bash
python main.py
```

Or run via command line with a specific query:

```bash
python main.py --query "Find molecules targeting KRAS G12D with good BBB permeability"
```

For batch processing:

```bash
python main.py --batch-file queries.txt --output-dir results/
```

## Project Structure

```
fade/
├── agents/               # Specialized agents for different tasks
├── workflows/            # Pipeline orchestration and job templates
├── data/                 # Input/output data storage
│   ├── inputs/           # User inputs and configurations
│   ├── outputs/          # Generated results
│   └── examples/         # Example configurations
├── models/               # ML models for property prediction
├── utils/                # Utility functions and helpers
├── configs/              # Configuration files
├── logs/                 # Log files
├── requirements.txt      # Python dependencies
├── setup.py              # Package setup
├── .env                  # Environment variables (not in version control)
└── main.py               # Main entry point
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
  url = {https://github.com/yourusername/fade}
}
```

## Acknowledgments

- This project was developed on the Northeastern University HPC cluster
- Special thanks to the contributors of AlphaFold, RDKit, and Schrodinger
