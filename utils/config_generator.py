"""
Configuration file generator for F.A.D.E framework.
"""

import os
import json
import yaml
from typing import Any, Dict, List, Optional
import jinja2


class ConfigGenerator:
    """
    Generates configuration files and job scripts for downstream processes
    in the F.A.D.E framework.
    """
    
    def __init__(self, template_dir: Optional[str] = None) -> None:
        """
        Initialize the configuration generator.
        
        Args:
            template_dir: Directory containing template files. If not provided,
                          defaults to the 'templates' directory within the framework.
        """
        if template_dir is None:
            # Get the directory where this script is located
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # Go one level up to find the framework root directory
            root_dir = os.path.dirname(current_dir)
            template_dir = os.path.join(root_dir, "workflows", "job_templates")
        
        self.template_dir = template_dir
        self.jinja_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir),
            trim_blocks=True,
            lstrip_blocks=True
        )
    
    def create_alphafold_config(
        self, 
        protein_name: str,
        sequence_file: str,
        output_dir: str,
        max_template_date: str = "2025-01-01",
        use_gpu: bool = True,
        num_recycle: int = 3,
        enable_amber_relax: bool = True,
    ) -> Dict[str, Any]:
        """
        Create configuration for AlphaFold structure prediction.
        
        Args:
            protein_name: Name of the protein target.
            sequence_file: Path to the FASTA file containing the protein sequence.
            output_dir: Directory where AlphaFold output will be saved.
            max_template_date: Maximum date for template consideration.
            use_gpu: Whether to use GPU for prediction.
            num_recycle: Number of recycle iterations for AlphaFold.
            enable_amber_relax: Whether to perform AMBER relaxation on the predicted structure.
            
        Returns:
            AlphaFold configuration dictionary.
        """
        config = {
            "protein_name": protein_name,
            "sequence_file": sequence_file,
            "output_dir": output_dir,
            "max_template_date": max_template_date,
            "use_gpu": use_gpu,
            "model_preset": "monomer",
            "db_preset": "full_dbs",
            "num_recycle": num_recycle,
            "enable_amber_relax": enable_amber_relax,
        }
        
        return config
    
    def create_docking_config(
        self,
        receptor_file: str,
        ligands_file: str,
        output_dir: str,
        binding_site_center: Optional[List[float]] = None,
        binding_site_size: Optional[List[float]] = None,
        exhaustiveness: int = 8,
        num_modes: int = 9,
    ) -> Dict[str, Any]:
        """
        Create configuration for molecular docking.
        
        Args:
            receptor_file: Path to the receptor PDB file.
            ligands_file: Path to the ligands file (SDF or MOL2).
            output_dir: Directory where docking output will be saved.
            binding_site_center: Optional XYZ coordinates of the binding site center.
            binding_site_size: Optional XYZ dimensions of the binding site box.
            exhaustiveness: Search exhaustiveness parameter.
            num_modes: Number of binding modes to generate.
            
        Returns:
            Docking configuration dictionary.
        """
        config = {
            "receptor_file": receptor_file,
            "ligands_file": ligands_file,
            "output_dir": output_dir,
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
        }
        
        if binding_site_center and binding_site_size:
            config["binding_site_center"] = binding_site_center
            config["binding_site_size"] = binding_site_size
        
        return config
    
    def save_config(
        self, 
        config: Dict[str, Any], 
        output_path: str, 
        format: str = "json"
    ) -> None:
        """
        Save a configuration dictionary to a file.
        
        Args:
            config: Configuration dictionary to save.
            output_path: Path where the configuration file will be saved.
            format: Format of the configuration file ('json' or 'yaml').
        """
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        if format.lower() == "json":
            with open(output_path, "w") as f:
                json.dump(config, f, indent=2)
        elif format.lower() == "yaml":
            with open(output_path, "w") as f:
                yaml.dump(config, f, default_flow_style=False)
        else:
            raise ValueError(f"Unsupported configuration format: {format}")
    
    def generate_job_script(
        self, 
        template_name: str, 
        output_path: str, 
        **kwargs
    ) -> None:
        """
        Generate a job script from a template.
        
        Args:
            template_name: Name of the template file.
            output_path: Path where the job script will be saved.
            **kwargs: Variables to be passed to the template.
        """
        template = self.jinja_env.get_template(template_name)
        script_content = template.render(**kwargs)
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, "w") as f:
            f.write(script_content)
        
        # Make the script executable
        os.chmod(output_path, 0o755)
    
    def create_alphafold_job(
        self,
        protein_name: str,
        sequence_file: str,
        output_dir: str,
        config_file: str,
        job_script_path: str,
        partition: str = "gpu",
        memory: str = "16G",
        cpus: int = 8,
        gpus: int = 1,
        time: str = "24:00:00",
        email: Optional[str] = None,
    ) -> None:
        """
        Create AlphaFold job configuration and script.
        
        Args:
            protein_name: Name of the protein target.
            sequence_file: Path to the FASTA file containing the protein sequence.
            output_dir: Directory where AlphaFold output will be saved.
            config_file: Path where the configuration file will be saved.
            job_script_path: Path where the job script will be saved.
            partition: SLURM partition to use.
            memory: Memory allocation for the job.
            cpus: Number of CPUs to request.
            gpus: Number of GPUs to request.
            time: Time limit for the job.
            email: Email address for job notifications.
        """
        # Create AlphaFold configuration
        config = self.create_alphafold_config(
            protein_name=protein_name,
            sequence_file=sequence_file,
            output_dir=output_dir
        )
        
        # Save configuration file
        self.save_config(config, config_file)
        
        # Generate job script
        self.generate_job_script(
            template_name="alphafold_job.sh.j2",
            output_path=job_script_path,
            protein_name=protein_name,
            config_file=config_file,
            partition=partition,
            memory=memory,
            cpus=cpus,
            gpus=gpus,
            time=time,
            email=email
        )
    
    def create_docking_job(
        self,
        receptor_name: str,
        receptor_file: str,
        ligands_file: str,
        output_dir: str,
        config_file: str,
        job_script_path: str,
        binding_site_center: Optional[List[float]] = None,
        binding_site_size: Optional[List[float]] = None,
        partition: str = "cpu",
        memory: str = "8G",
        cpus: int = 16,
        time: str = "12:00:00",
        email: Optional[str] = None,
    ) -> None:
        """
        Create molecular docking job configuration and script.
        
        Args:
            receptor_name: Name of the receptor protein.
            receptor_file: Path to the receptor PDB file.
            ligands_file: Path to the ligands file (SDF or MOL2).
            output_dir: Directory where docking output will be saved.
            config_file: Path where the configuration file will be saved.
            job_script_path: Path where the job script will be saved.
            binding_site_center: Optional XYZ coordinates of the binding site center.
            binding_site_size: Optional XYZ dimensions of the binding site box.
            partition: SLURM partition to use.
            memory: Memory allocation for the job.
            cpus: Number of CPUs to request.
            time: Time limit for the job.
            email: Email address for job notifications.
        """
        # Create docking configuration
        config = self.create_docking_config(
            receptor_file=receptor_file,
            ligands_file=ligands_file,
            output_dir=output_dir,
            binding_site_center=binding_site_center,
            binding_site_size=binding_site_size
        )
        
        # Save configuration file
        self.save_config(config, config_file)
        
        # Generate job script
        self.generate_job_script(
            template_name="docking_job.sh.j2",
            output_path=job_script_path,
            receptor_name=receptor_name,
            config_file=config_file,
            partition=partition,
            memory=memory,
            cpus=cpus,
            time=time,
            email=email
        )
