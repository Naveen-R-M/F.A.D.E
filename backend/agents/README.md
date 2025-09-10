# Agents Directory

## Core Architecture

### Multi-Agent System Structure:

- **Target Selector Agent**: Processes natural language queries to identify protein targets and extract requirements.
- **Structure Predictor Agent**: Uses AlphaFold3 to generate 3D protein structures.
- **Molecule Generator Agent**: Creates potential drug candidates using AI and chemical libraries.
- **Evaluator Agent**: Assesses drug-like properties (ADMET, Lipinski's Rule of Five, etc.).
- **Docking Agent**: Simulates protein-ligand binding using Schrodinger Glide.
- **Refiner Agent**: Iteratively improves molecules based on evaluation results.
- **Memory Manager**: Tracks history and learns from previous