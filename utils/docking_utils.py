"""
Docking Utility Functions - Helper functions for molecular docking.
"""

import os
import logging
import tempfile
import subprocess
from typing import Dict, Any, List, Optional, Tuple

logger = logging.getLogger(__name__)

class DockingUtils:
    """
    Utility class for molecular docking.
    
    This class provides utilities for:
    1. Preparing proteins and ligands for docking
    2. Running docking simulations
    3. Analyzing docking results
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the DockingUtils class.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.tool = config.get("pipeline", {}).get("docking", {}).get("tool", "glide")
        self.schrodinger_path = config.get("paths", {}).get("schrodinger", "/shared/EL9/explorer/schrodinger/2024-4")
        self.precision = config.get("pipeline", {}).get("docking", {}).get("precision", "SP")
        self.exhaustiveness = config.get("pipeline", {}).get("docking", {}).get("exhaustiveness", 8)
        
        logger.info(f"DockingUtils initialized with {self.tool}")
    
    def prepare_receptor(self, structure_file: str, binding_site: Optional[Dict[str, Any]] = None) -> str:
        """
        Prepare a protein structure for docking.
        
        Args:
            structure_file: Path to PDB file
            binding_site: Optional binding site information
            
        Returns:
            Path to prepared receptor file
        """
        if not os.path.exists(structure_file):
            logger.error(f"Structure file {structure_file} not found")
            raise FileNotFoundError(f"Structure file {structure_file} not found")
        
        output_dir = os.path.dirname(structure_file)
        base_name = os.path.splitext(os.path.basename(structure_file))[0]
        
        if self.tool == "glide":
            return self._prepare_receptor_glide(structure_file, output_dir, base_name, binding_site)
        elif self.tool == "autodock":
            return self._prepare_receptor_autodock(structure_file, output_dir, base_name, binding_site)
        else:
            logger.error(f"Unsupported docking tool: {self.tool}")
            raise ValueError(f"Unsupported docking tool: {self.tool}")
    
    def _prepare_receptor_glide(self, structure_file: str, output_dir: str, base_name: str, binding_site: Optional[Dict[str, Any]]) -> str:
        """Prepare a receptor for Glide docking."""
        # In a real implementation, this would run Protein Preparation Wizard and Receptor Grid Generation
        # For this mock, we'll just create a fake grid file
        
        logger.info(f"Preparing receptor for Glide: {structure_file}")
        
        # Mock grid file
        grid_file = os.path.join(output_dir, f"{base_name}_grid.zip")
        
        # In a real implementation, we would run something like:
        # cmd = f"{self.schrodinger_path}/prepwizard -WAIT {structure_file} {output_dir}/{base_name}_prepared.mae"
        # subprocess.run(cmd, shell=True, check=True)
        # 
        # grid_cmd = f"{self.schrodinger_path}/glide_grid -WAIT -rec {output_dir}/{base_name}_prepared.mae -out {grid_file}"
        # subprocess.run(grid_cmd, shell=True, check=True)
        
        # Mock implementation
        with open(grid_file, 'w') as f:
            f.write("MOCK GLIDE GRID FILE\n")
            f.write(f"Original PDB: {structure_file}\n")
            if binding_site:
                f.write(f"Binding site: {binding_site}\n")
        
        logger.info(f"Receptor prepared, grid file: {grid_file}")
        return grid_file
    
    def _prepare_receptor_autodock(self, structure_file: str, output_dir: str, base_name: str, binding_site: Optional[Dict[str, Any]]) -> str:
        """Prepare a receptor for AutoDock Vina."""
        # In a real implementation, this would prepare PDBQT files and config
        # For this mock, we'll just create a fake PDBQT file
        
        logger.info(f"Preparing receptor for AutoDock: {structure_file}")
        
        # Mock PDBQT file
        pdbqt_file = os.path.join(output_dir, f"{base_name}.pdbqt")
        
        # In a real implementation, we would run something like:
        # cmd = f"prepare_receptor -r {structure_file} -o {pdbqt_file}"
        # subprocess.run(cmd, shell=True, check=True)
        
        # Mock implementation
        with open(pdbqt_file, 'w') as f:
            f.write("MOCK AUTODOCK PDBQT FILE\n")
            f.write(f"Original PDB: {structure_file}\n")
            if binding_site:
                f.write(f"Binding site: {binding_site}\n")
        
        logger.info(f"Receptor prepared, PDBQT file: {pdbqt_file}")
        return pdbqt_file
    
    def prepare_ligand(self, mol: Any, output_dir: str, name: str) -> str:
        """
        Prepare a ligand for docking.
        
        Args:
            mol: Molecule object (RDKit or similar)
            output_dir: Output directory
            name: Base name for output files
            
        Returns:
            Path to prepared ligand file
        """
        os.makedirs(output_dir, exist_ok=True)
        
        if self.tool == "glide":
            return self._prepare_ligand_glide(mol, output_dir, name)
        elif self.tool == "autodock":
            return self._prepare_ligand_autodock(mol, output_dir, name)
        else:
            logger.error(f"Unsupported docking tool: {self.tool}")
            raise ValueError(f"Unsupported docking tool: {self.tool}")
    
    def _prepare_ligand_glide(self, mol: Any, output_dir: str, name: str) -> str:
        """Prepare a ligand for Glide docking."""
        # In a real implementation, this would use LigPrep
        # For this mock, we'll just create a fake MAE file
        
        logger.info(f"Preparing ligand for Glide: {name}")
        
        # Mock MAE file
        mae_file = os.path.join(output_dir, f"{name}.mae")
        
        # In a real implementation, we would run something like:
        # cmd = f"{self.schrodinger_path}/ligprep -WAIT -i {sdf_file} -o {mae_file}"
        # subprocess.run(cmd, shell=True, check=True)
        
        # Mock implementation
        with open(mae_file, 'w') as f:
            f.write("MOCK GLIDE MAE FILE\n")
            f.write(f"Molecule: {mol.get('smiles', '')}\n")
        
        logger.info(f"Ligand prepared, MAE file: {mae_file}")
        return mae_file
    
    def _prepare_ligand_autodock(self, mol: Any, output_dir: str, name: str) -> str:
        """Prepare a ligand for AutoDock Vina."""
        # In a real implementation, this would prepare PDBQT files
        # For this mock, we'll just create a fake PDBQT file
        
        logger.info(f"Preparing ligand for AutoDock: {name}")
        
        # Mock PDBQT file
        pdbqt_file = os.path.join(output_dir, f"{name}.pdbqt")
        
        # In a real implementation, we would run something like:
        # cmd = f"prepare_ligand -l {sdf_file} -o {pdbqt_file}"
        # subprocess.run(cmd, shell=True, check=True)
        
        # Mock implementation
        with open(pdbqt_file, 'w') as f:
            f.write("MOCK AUTODOCK PDBQT FILE\n")
            f.write(f"Molecule: {mol.get('smiles', '')}\n")
        
        logger.info(f"Ligand prepared, PDBQT file: {pdbqt_file}")
        return pdbqt_file
    
    def run_docking(self, receptor_file: str, ligand_file: str, output_dir: str, name: str) -> Dict[str, Any]:
        """
        Run docking simulation.
        
        Args:
            receptor_file: Path to prepared receptor file
            ligand_file: Path to prepared ligand file
            output_dir: Output directory
            name: Base name for output files
            
        Returns:
            Dictionary with docking results
        """
        os.makedirs(output_dir, exist_ok=True)
        
        if self.tool == "glide":
            return self._run_docking_glide(receptor_file, ligand_file, output_dir, name)
        elif self.tool == "autodock":
            return self._run_docking_autodock(receptor_file, ligand_file, output_dir, name)
        else:
            logger.error(f"Unsupported docking tool: {self.tool}")
            raise ValueError(f"Unsupported docking tool: {self.tool}")
    
    def _run_docking_glide(self, grid_file: str, ligand_file: str, output_dir: str, name: str) -> Dict[str, Any]:
        """Run Glide docking."""
        # In a real implementation, this would run Glide
        # For this mock, we'll just create fake output files
        
        logger.info(f"Running Glide docking: {name}")
        
        # Create a Glide input file
        input_file = os.path.join(output_dir, f"{name}_glide.in")
        with open(input_file, 'w') as f:
            f.write(f"GRIDFILE {grid_file}\n")
            f.write(f"LIGANDFILE {ligand_file}\n")
            f.write(f"PRECISION {self.precision}\n")
            f.write(f"WRITE_CSV 1\n")
            f.write(f"POSES_PER_LIG 10\n")
            
        # Output files
        pose_file = os.path.join(output_dir, f"{name}_poses.mae")
        score_file = os.path.join(output_dir, f"{name}_scores.csv")
        
        # In a real implementation, we would run something like:
        # cmd = f"{self.schrodinger_path}/glide -WAIT -i {input_file} -o {pose_file}"
        # subprocess.run(cmd, shell=True, check=True)
        
        # Mock implementation
        with open(pose_file, 'w') as f:
            f.write("MOCK GLIDE POSE FILE\n")
        
        with open(score_file, 'w') as f:
            f.write("MOCK GLIDE SCORE FILE\n")
            f.write("Title,Binding Score,Pose Number\n")
            f.write(f"{name},-9.5,1\n")
        
        logger.info(f"Docking completed, pose file: {pose_file}")
        
        # Mock docking results
        return {
            "score": -9.5,
            "pose_file": pose_file,
            "score_file": score_file,
            "binding_site": "GTP binding pocket",
            "interactions": [
                {"type": "hydrogen_bond", "residue": "Asp57"},
                {"type": "pi_stacking", "residue": "Tyr32"}
            ]
        }
    
    def _run_docking_autodock(self, receptor_file: str, ligand_file: str, output_dir: str, name: str) -> Dict[str, Any]:
        """Run AutoDock Vina docking."""
        # In a real implementation, this would run AutoDock Vina
        # For this mock, we'll just create fake output files
        
        logger.info(f"Running AutoDock docking: {name}")
        
        # Create a Vina config file
        config_file = os.path.join(output_dir, f"{name}_vina.cfg")
        with open(config_file, 'w') as f:
            f.write(f"receptor = {receptor_file}\n")
            f.write(f"ligand = {ligand_file}\n")
            f.write(f"center_x = 0.0\n")
            f.write(f"center_y = 0.0\n")
            f.write(f"center_z = 0.0\n")
            f.write(f"size_x = 20.0\n")
            f.write(f"size_y = 20.0\n")
            f.write(f"size_z = 20.0\n")
            f.write(f"exhaustiveness = {self.exhaustiveness}\n")
            
        # Output files
        output_file = os.path.join(output_dir, f"{name}_out.pdbqt")
        log_file = os.path.join(output_dir, f"{name}_log.txt")
        
        # In a real implementation, we would run something like:
        # cmd = f"vina --config {config_file} --out {output_file} --log {log_file}"
        # subprocess.run(cmd, shell=True, check=True)
        
        # Mock implementation
        with open(output_file, 'w') as f:
            f.write("MOCK AUTODOCK OUTPUT FILE\n")
        
        with open(log_file, 'w') as f:
            f.write("MOCK AUTODOCK LOG FILE\n")
            f.write("Affinity: -8.7 (kcal/mol)\n")
        
        logger.info(f"Docking completed, output file: {output_file}")
        
        # Mock docking results
        return {
            "score": -8.7,
            "pose_file": output_file,
            "log_file": log_file,
            "binding_site": "GTP binding pocket",
            "interactions": [
                {"type": "hydrogen_bond", "residue": "Gln61"},
                {"type": "hydrophobic", "residue": "Val12"}
            ]
        }
    
    def analyze_docking_results(self, docking_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze docking results.
        
        Args:
            docking_results: Docking results from run_docking
            
        Returns:
            Dictionary with analysis results
        """
        # In a real implementation, this would analyze poses, interactions, etc.
        # For this mock, we'll just return the results
        
        logger.info("Analyzing docking results")
        
        # Mock analysis
        docking_results["analysis"] = {
            "quality": "good" if docking_results.get("score", 0) < -8.0 else "moderate",
            "key_interactions": [i.get("residue") for i in docking_results.get("interactions", [])],
            "potential_improvements": [
                "Add hydrogen bond donor near residue X",
                "Extend hydrophobic group into pocket Y"
            ]
        }
        
        return docking_results
    
    def create_docking_job_script(self, receptor_file: str, ligand_file: str, output_dir: str, name: str) -> str:
        """
        Create a SLURM job script for docking.
        
        Args:
            receptor_file: Path to prepared receptor file
            ligand_file: Path to prepared ligand file
            output_dir: Output directory
            name: Base name for output files
            
        Returns:
            Path to job script
        """
        slurm_config = self.config.get("slurm", {}).get("docking", {})
        
        script_content = ""
        
        if self.tool == "glide":
            # Create a Glide input file
            input_file = os.path.join(output_dir, f"{name}_glide.in")
            with open(input_file, 'w') as f:
                f.write(f"GRIDFILE {receptor_file}\n")
                f.write(f"LIGANDFILE {ligand_file}\n")
                f.write(f"PRECISION {self.precision}\n")
                f.write(f"WRITE_CSV 1\n")
                f.write(f"POSES_PER_LIG 10\n")
            
            script_content = f"""#!/bin/bash
#SBATCH --partition={slurm_config.get("partition", "compute")}
#SBATCH --cpus-per-task={slurm_config.get("cpus_per_task", 8)}
#SBATCH --mem={slurm_config.get("mem", "16G")}
#SBATCH --time={slurm_config.get("time", "02:00:00")}
#SBATCH --job-name=glide_{name}
#SBATCH --output={output_dir}/glide_{name}_%j.out
#SBATCH --error={output_dir}/glide_{name}_%j.err

# Load modules
module load schrodinger/2024-4

# Run Glide docking
$SCHRODINGER/glide \\
  -JOBNAME {name} \\
  -HOST localhost \\
  -WAIT \\
  -LOCAL \\
  -NOJOBID \\
  {input_file}

# Exit with the command's exit code
exit $?
"""
        elif self.tool == "autodock":
            # Create a Vina config file
            config_file = os.path.join(output_dir, f"{name}_vina.cfg")
            with open(config_file, 'w') as f:
                f.write(f"receptor = {receptor_file}\n")
                f.write(f"ligand = {ligand_file}\n")
                f.write(f"center_x = 0.0\n")
                f.write(f"center_y = 0.0\n")
                f.write(f"center_z = 0.0\n")
                f.write(f"size_x = 20.0\n")
                f.write(f"size_y = 20.0\n")
                f.write(f"size_z = 20.0\n")
                f.write(f"exhaustiveness = {self.exhaustiveness}\n")
            
            script_content = f"""#!/bin/bash
#SBATCH --partition={slurm_config.get("partition", "compute")}
#SBATCH --cpus-per-task={slurm_config.get("cpus_per_task", 8)}
#SBATCH --mem={slurm_config.get("mem", "16G")}
#SBATCH --time={slurm_config.get("time", "02:00:00")}
#SBATCH --job-name=vina_{name}
#SBATCH --output={output_dir}/vina_{name}_%j.out
#SBATCH --error={output_dir}/vina_{name}_%j.err

# Load modules
module load miniconda3/24.11.1
source activate $SCRATCH/conda-envs/autodock

# Run AutoDock Vina
vina \\
  --config {config_file} \\
  --out {output_dir}/{name}_out.pdbqt \\
  --log {output_dir}/{name}_log.txt

# Exit with the command's exit code
exit $?
"""
        
        # Write script to temporary file
        fd, script_path = tempfile.mkstemp(suffix=".sh")
        with os.fdopen(fd, "w") as f:
            f.write(script_content)
        
        return script_path
