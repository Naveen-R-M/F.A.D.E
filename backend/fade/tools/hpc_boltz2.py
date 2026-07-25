"""
HPC Boltz2 client for protein structure prediction and binding-affinity screening.

This module provides functions to run Boltz2 on the HPC cluster for:
  1. De novo protein structure prediction (predict_structure / wait_for_structure)
  2. Protein-ligand binding-affinity prediction (predict_affinity / wait_for_affinity)

The affinity mode is the shared screening engine used by all three branches of the
F.A.D.E pipeline (holo / apo / de-novo). Once each branch has converged on a protein
structure, a selected pocket, and a set of candidate ligand SMILES, screening is
identical regardless of how the structure was obtained.
"""

import os
import time
import json
import math
import logging
import tempfile
from typing import Dict, Any, Optional, Tuple, List
from pathlib import Path

import paramiko
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_boltz2")


class HPCBoltz2Client:
    """Client for running Boltz2 structure prediction on HPC."""
    
    def __init__(self, 
                 host: str = None,
                 username: str = None,
                 password: str = None,
                 key_path: str = None):
        """
        Initialize HPC Boltz2 client.
        
        Args:
            host: HPC host address
            username: SSH username
            password: SSH password (if not using key)
            key_path: Path to SSH private key
        """
        self.host = host or os.getenv("HPC_HOST", "explorer.northeastern.edu")
        self.username = username or os.getenv("HPC_USER", "rajagopalmohanraj.n")
        self.password = password or os.getenv("HPC_PASSWORD")
        self.key_path = key_path or os.getenv("HPC_KEY_PATH")
        
        # HPC paths
        self.remote_base = f"/scratch/{self.username}/F.A.D.E/boltz2_jobs"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles/"
        
        # Boltz2 settings
        self.boltz2_module = "boltz/2"
        self.max_wait_time = config.BOLTZ2_MAX_WAIT_TIME if hasattr(config, 'BOLTZ2_MAX_WAIT_TIME') else 1800
        self.poll_interval = config.BOLTZ2_POLL_INTERVAL if hasattr(config, 'BOLTZ2_POLL_INTERVAL') else 30
        self.partition = config.BOLTZ2_GPU_PARTITION if hasattr(config, 'BOLTZ2_GPU_PARTITION') else "gpu"
        self.slurm_account = config.BOLTZ2_SLURM_ACCOUNT if hasattr(config, 'BOLTZ2_SLURM_ACCOUNT') else None
        self.use_cpu = config.BOLTZ2_USE_CPU if hasattr(config, 'BOLTZ2_USE_CPU') else False
        self.BOLTZ_CACHE_DIR = config.BOLTZ_CACHE_DIR if hasattr(config, 'BOLTZ_CACHE_DIR') else "/scratch/rajagopalmohanraj.n/boltz_cache"

        # Boltz2 affinity-screening settings (fall back to structure settings when unset)
        self.affinity_remote_base = f"/scratch/{self.username}/F.A.D.E/boltz2_affinity_jobs"
        self.affinity_max_wait_time = config.BOLTZ2_AFFINITY_MAX_WAIT_TIME if hasattr(config, 'BOLTZ2_AFFINITY_MAX_WAIT_TIME') else self.max_wait_time
        self.affinity_poll_interval = config.BOLTZ2_AFFINITY_POLL_INTERVAL if hasattr(config, 'BOLTZ2_AFFINITY_POLL_INTERVAL') else self.poll_interval
        self.affinity_partition = config.BOLTZ2_AFFINITY_GPU_PARTITION if hasattr(config, 'BOLTZ2_AFFINITY_GPU_PARTITION') else self.partition

        self.ssh_client = None
        self.sftp_client = None
    
    def predict_structure(self, 
                         sequence: str,
                         job_name: str = None,
                         job_id: str = None,
                         use_msa: bool = True,
                         num_models: int = 1) -> Dict[str, Any]:
        """
        Predict protein structure using Boltz2.
        
        Args:
            sequence: Protein sequence (amino acids)
            job_name: Optional job name
            job_id: Optional job ID (reuse existing)
            use_msa: Whether to use MSA (Multiple Sequence Alignment)
            num_models: Number of models to generate
            
        Returns:
            Dictionary with job_id and status
        """
        # Generate job ID if not provided
        if not job_id:
            job_id = f"boltz2_{int(time.time())}"
        
        if not job_name:
            job_name = f"fade_boltz2_{job_id}"
        
        logger.info(f"Running Boltz2 prediction with job ID: {job_id}")
        
        try:
            # Connect to HPC
            self._connect()
            
            # Create remote directory
            remote_dir = f"{self.remote_base}/{job_id}"
            self._create_remote_directory(remote_dir)
            
            # Create and upload YAML file (Boltz2 likely expects YAML like AlphaFold3)
            yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {sequence}
"""
            yaml_path = f"{remote_dir}/input.yaml"
            self._upload_content(yaml_content, yaml_path)
            
            # Create and upload SLURM script
            script_content = self._create_boltz2_script(
                job_id=job_id,
                job_name=job_name,
                input_file="input.yaml",
                remote_dir=remote_dir,
                use_msa=use_msa,
                num_models=num_models
            )
            script_path = f"{remote_dir}/run_boltz2.sh"
            self._upload_content(script_content, script_path)
            
            # Submit job
            slurm_job_id = self._submit_job(script_path)
            
            logger.info(f"Boltz2 job submitted with SLURM ID: {slurm_job_id}")
            
            return {
                "job_id": job_id,
                "slurm_job_id": slurm_job_id,
                "status": "submitted",
                "remote_dir": remote_dir
            }
            
        except Exception as e:
            logger.error(f"Failed to submit Boltz2 job: {e}")
            raise
        finally:
            self._disconnect()
    
    def wait_for_structure(self, 
                          job_id: str,
                          max_wait: int = None,
                          poll_interval: int = None) -> Optional[str]:
        """
        Wait for Boltz2 job to complete and return structure.
        
        Args:
            job_id: Job ID to wait for
            max_wait: Maximum wait time in seconds
            poll_interval: Polling interval in seconds
            
        Returns:
            PDB content as string or None if failed
        """
        max_wait = max_wait or self.max_wait_time
        poll_interval = poll_interval or self.poll_interval
        
        logger.info(f"Waiting for Boltz2 job {job_id} (max {max_wait}s)")
        
        start_time = time.time()
        
        try:
            self._connect()
            remote_dir = f"{self.remote_base}/{job_id}"
            
            while time.time() - start_time < max_wait:
                # Check for different possible output file locations
                # Boltz2 creates files in output/boltz_results_input/predictions/input/
                possible_outputs = [
                    f"{remote_dir}/output/predicted_structure.cif",  # Our copied file
                    f"{remote_dir}/output/boltz_results_input/predictions/input/input_model_0.cif",  # Boltz2 default
                    f"{remote_dir}/output/protein.cif",
                    f"{remote_dir}/output/protein.pdb",
                ]
                
                for output_file in possible_outputs:
                    if self._file_exists(output_file):
                        logger.info(f"Boltz2 structure ready for job {job_id} at {output_file}")
                        
                        # Download structure file
                        structure_content = self._download_file(output_file)
                        
                        # Get confidence scores if available
                        confidence_files = [
                            f"{remote_dir}/output/confidence_scores.json",
                            f"{remote_dir}/output/boltz_results_input/predictions/input/confidence_input_model_0.json",
                        ]
                        confidence_data = None
                        for conf_file in confidence_files:
                            if self._file_exists(conf_file):
                                confidence_data = self._download_file(conf_file)
                                logger.info(f"Downloaded confidence scores from {conf_file}")
                                break
                        
                        return structure_content
                
                # Check for error
                error_file = f"{remote_dir}/boltz2.err"
                if self._file_exists(error_file):
                    error_content = self._download_file(error_file)
                    if error_content and len(error_content.strip()) > 0:
                        logger.error(f"Boltz2 error: {error_content[:500]}")
                        return None
                
                logger.debug(f"Waiting for Boltz2... ({int(time.time() - start_time)}s elapsed)")
                time.sleep(poll_interval)
            
            logger.error(f"Boltz2 job {job_id} timed out after {max_wait}s")
            return None
            
        except Exception as e:
            logger.error(f"Error waiting for Boltz2: {e}")
            return None
        finally:
            self._disconnect()
    
    def check_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check status of a Boltz2 job.
        
        Args:
            job_id: Job ID to check
            
        Returns:
            Status dictionary
        """
        try:
            self._connect()
            remote_dir = f"{self.remote_base}/{job_id}"
            
            status = {
                "job_id": job_id,
                "status": "unknown",
                "has_output": False,
                "has_error": False
            }
            
            # Check for output - Boltz2 creates in boltz_results_input directory
            possible_output_files = [
                f"{remote_dir}/output/predicted_structure.cif",
                f"{remote_dir}/output/boltz_results_input/predictions/input/input_model_0.cif",
            ]
            
            for output_file in possible_output_files:
                if self._file_exists(output_file):
                    status["status"] = "completed"
                    status["has_output"] = True
                    status["output_file"] = output_file
                    break
            
            # Check for error
            error_file = f"{remote_dir}/boltz2.err"
            if self._file_exists(error_file):
                error_content = self._download_file(error_file)
                if error_content and len(error_content.strip()) > 0:
                    status["status"] = "failed"
                    status["has_error"] = True
                    status["error_message"] = error_content[:500]
            
            # Check if still running
            log_file = f"{remote_dir}/boltz2.out"
            if self._file_exists(log_file) and status["status"] == "unknown":
                status["status"] = "running"
            
            return status
            
        except Exception as e:
            logger.error(f"Error checking job status: {e}")
            return {"job_id": job_id, "status": "error", "error": str(e)}
        finally:
            self._disconnect()
    
    # ------------------------------------------------------------------ #
    # Binding-affinity screening                                          #
    # ------------------------------------------------------------------ #
    def predict_affinity(self,
                         sequence: str,
                         ligand_smiles: str,
                         job_name: str = None,
                         job_id: str = None,
                         use_msa: bool = True,
                         ligand_id: str = "LIG",
                         protein_chain: str = "A") -> Dict[str, Any]:
        """
        Submit a Boltz2 protein-ligand binding-affinity prediction job.

        Boltz2 scores exactly ONE ligand per job, so batch screening must submit
        one job per candidate molecule. This mirrors ``predict_structure`` but the
        input YAML additionally carries a ligand block and a ``properties: affinity``
        request against that ligand.

        Args:
            sequence: Target protein sequence (amino acids).
            ligand_smiles: Candidate ligand as a SMILES string.
            job_name: Optional SLURM job name.
            job_id: Optional job ID (generated if omitted).
            use_msa: Whether to use the MSA server.
            ligand_id: Chain/entity id used for the ligand in the YAML.
            protein_chain: Chain id used for the protein in the YAML.

        Returns:
            Dictionary with job_id, slurm_job_id, status and remote_dir.
        """
        if not job_id:
            job_id = f"boltz2aff_{int(time.time())}"
        if not job_name:
            job_name = f"fade_boltz2aff_{job_id}"

        logger.info(f"Running Boltz2 affinity prediction with job ID: {job_id}")

        try:
            self._connect()

            remote_dir = f"{self.affinity_remote_base}/{job_id}"
            self._create_remote_directory(remote_dir)

            # Escape any embedded quotes in the SMILES for YAML safety.
            safe_smiles = ligand_smiles.replace('"', '\\"')

            # Boltz2 affinity YAML: protein + ligand + affinity property.
            yaml_content = f"""version: 1
sequences:
  - protein:
      id: {protein_chain}
      sequence: {sequence}
  - ligand:
      id: {ligand_id}
      smiles: "{safe_smiles}"
properties:
  - affinity:
      binder: {ligand_id}
"""
            yaml_path = f"{remote_dir}/input.yaml"
            self._upload_content(yaml_content, yaml_path)

            script_content = self._create_affinity_script(
                job_id=job_id,
                job_name=job_name,
                input_file="input.yaml",
                remote_dir=remote_dir,
                use_msa=use_msa,
            )
            script_path = f"{remote_dir}/run_boltz2_affinity.sh"
            self._upload_content(script_content, script_path)

            slurm_job_id = self._submit_job(script_path)

            logger.info(f"Boltz2 affinity job submitted with SLURM ID: {slurm_job_id}")

            return {
                "job_id": job_id,
                "slurm_job_id": slurm_job_id,
                "status": "submitted",
                "remote_dir": remote_dir,
                "ligand_smiles": ligand_smiles,
            }

        except Exception as e:
            logger.error(f"Failed to submit Boltz2 affinity job: {e}")
            raise
        finally:
            self._disconnect()

    def wait_for_affinity(self,
                          job_id: str,
                          max_wait: int = None,
                          poll_interval: int = None) -> Optional[Dict[str, Any]]:
        """
        Wait for a Boltz2 affinity job to finish and return the parsed prediction.

        Args:
            job_id: Job ID returned by ``predict_affinity``.
            max_wait: Maximum wait time in seconds.
            poll_interval: Polling interval in seconds.

        Returns:
            Dict with the parsed affinity prediction, or None on timeout/failure:
              {
                "affinity_pred_value": float,      # Boltz2 raw: log10(IC50 in uM), lower = stronger
                "affinity_probability_binary": float,  # 0-1 probability the ligand is a binder
                "binding_affinity_kcal_mol": float,    # derived, informational (see _ic50_to_kcal)
                "predicted_ic50_uM": float,            # derived from affinity_pred_value
                "raw": {...},                          # full JSON payload from Boltz2
                "pose_path": Optional[str],            # local path to the predicted complex, if any
              }
        """
        max_wait = max_wait or self.affinity_max_wait_time
        poll_interval = poll_interval or self.affinity_poll_interval

        logger.info(f"Waiting for Boltz2 affinity job {job_id} (max {max_wait}s)")

        start_time = time.time()

        try:
            self._connect()
            remote_dir = f"{self.affinity_remote_base}/{job_id}"

            # Boltz2 writes affinity_<input_stem>.json next to the predicted structure.
            possible_affinity_files = [
                f"{remote_dir}/output/affinity_scores.json",  # our copied file
                f"{remote_dir}/output/boltz_results_input/predictions/input/affinity_input.json",
            ]
            possible_structure_files = [
                f"{remote_dir}/output/predicted_complex.cif",
                f"{remote_dir}/output/boltz_results_input/predictions/input/input_model_0.cif",
            ]

            while time.time() - start_time < max_wait:
                for affinity_file in possible_affinity_files:
                    if self._file_exists(affinity_file):
                        logger.info(f"Boltz2 affinity result ready for job {job_id} at {affinity_file}")
                        raw_json = self._download_file(affinity_file)
                        parsed = self._parse_affinity_json(raw_json)

                        # Best-effort: also grab the predicted complex for pose visualisation.
                        pose_path = None
                        for struct_file in possible_structure_files:
                            if self._file_exists(struct_file):
                                try:
                                    content = self._download_file(struct_file)
                                    pose_path = self._save_affinity_pose(content, job_id)
                                except Exception as pose_err:
                                    logger.warning(f"Could not download affinity pose: {pose_err}")
                                break
                        parsed["pose_path"] = pose_path
                        return parsed

                # Check for a fatal error early.
                error_file = f"{remote_dir}/boltz2_affinity.err"
                if self._file_exists(error_file):
                    error_content = self._download_file(error_file)
                    if error_content and len(error_content.strip()) > 0:
                        logger.error(f"Boltz2 affinity error: {error_content[:500]}")
                        return None

                logger.debug(f"Waiting for Boltz2 affinity... ({int(time.time() - start_time)}s elapsed)")
                time.sleep(poll_interval)

            logger.error(f"Boltz2 affinity job {job_id} timed out after {max_wait}s")
            return None

        except Exception as e:
            logger.error(f"Error waiting for Boltz2 affinity: {e}")
            return None
        finally:
            self._disconnect()

    def _parse_affinity_json(self, raw_json: str) -> Dict[str, Any]:
        """
        Parse Boltz2's affinity JSON payload into a normalised result dict.

        Boltz2 reports:
          - affinity_pred_value: predicted binding affinity as log10(IC50) with IC50
            in micromolar; LOWER means STRONGER binding.
          - affinity_probability_binary: probability (0-1) that the ligand is a binder
            (vs a decoy), used for hit discovery / ranking.

        Some Boltz2 builds also emit per-model ensemble variants
        (affinity_pred_value1/2, affinity_probability_binary1/2); we keep the primary
        value and expose the raw payload for callers that want the ensemble.
        """
        try:
            data = json.loads(raw_json)
        except Exception as e:
            logger.error(f"Failed to parse affinity JSON: {e}")
            return {
                "affinity_pred_value": None,
                "affinity_probability_binary": None,
                "binding_affinity_kcal_mol": None,
                "predicted_ic50_uM": None,
                "raw": {},
                "parse_error": str(e),
            }

        pred_value = data.get("affinity_pred_value")
        prob_binary = data.get("affinity_probability_binary")

        predicted_ic50_uM = None
        binding_affinity_kcal_mol = None
        if pred_value is not None:
            try:
                predicted_ic50_uM = float(10 ** float(pred_value))
                binding_affinity_kcal_mol = self._ic50_to_kcal(predicted_ic50_uM)
            except Exception as conv_err:
                logger.warning(f"Affinity value conversion failed: {conv_err}")

        return {
            "affinity_pred_value": pred_value,
            "affinity_probability_binary": prob_binary,
            "binding_affinity_kcal_mol": binding_affinity_kcal_mol,
            "predicted_ic50_uM": predicted_ic50_uM,
            "raw": data,
        }

    @staticmethod
    def _ic50_to_kcal(ic50_uM: float, temperature_k: float = 298.15) -> Optional[float]:
        """
        Convert an IC50 (micromolar) to an approximate binding free energy in kcal/mol.

        Uses dG = RT * ln(Kd) with the common approximation Kd ~= IC50. This is an
        INFORMATIONAL convenience only: Boltz2's native output is log10(IC50), and
        IC50 is assay-dependent, so treat kcal/mol here as a rough guide, not a
        calibrated free energy. Ranking should use affinity_pred_value /
        affinity_probability_binary directly.
        """
        try:
            if ic50_uM is None or ic50_uM <= 0:
                return None
            R = 1.987204259e-3  # kcal/(mol*K)
            ic50_molar = ic50_uM * 1e-6
            return R * temperature_k * math.log(ic50_molar)
        except Exception:
            return None

    def _save_affinity_pose(self, content: str, job_id: str) -> str:
        """Save a predicted protein-ligand complex to the local poses directory."""
        poses_dir = config.DATA_DIR / "poses"
        poses_dir.mkdir(parents=True, exist_ok=True)
        # Preserve the CIF extension Boltz2 emits.
        filepath = poses_dir / f"boltz2_affinity_{job_id}.cif"
        with open(filepath, "w") as f:
            f.write(content)
        logger.info(f"Saved affinity pose to {filepath}")
        return str(filepath)

    def _create_affinity_script(self,
                                job_id: str,
                                job_name: str,
                                input_file: str,
                                remote_dir: str,
                                use_msa: bool) -> str:
        """
        Create the SLURM script for a Boltz2 affinity job.

        Structurally identical to the structure-prediction script but writes to the
        affinity out/err files and copies out the affinity JSON in addition to the
        predicted complex.
        """
        device_val = "cpu" if self.use_cpu else "gpu"
        msa_arg_val = "--use_msa_server" if use_msa else ""

        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={remote_dir}/boltz2_affinity.out",
            f"#SBATCH --error={remote_dir}/boltz2_affinity.err",
            f"#SBATCH --partition={self.affinity_partition}",
        ]

        if self.slurm_account:
            script_lines.append(f"#SBATCH --account={self.slurm_account}")

        if self.use_cpu:
            script_lines.extend([
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks=4",
                "#SBATCH --cpus-per-task=4",
                "#SBATCH --mem=32G",
            ])
        else:
            script_lines.extend([
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks-per-node=4",
                "#SBATCH --mem=32G",
                "#SBATCH --gres=gpu:1",
            ])

        script_lines.extend([
            "#SBATCH --time=02:00:00",
            "",
            "# Load modules",
            f"module use {self.module_path}",
            f"module load {self.boltz2_module}",
            "",
            "# Set working directory",
            f"cd {remote_dir}",
            "",
            "# Create output directory",
            "mkdir -p output",
            "",
            "echo \"Input file content:\"",
            f"head -20 {input_file}",
            "echo \"\"",
            "echo \"Starting Boltz2 affinity prediction...\"",
            f"echo \"Device: {device_val}\"",
            f"echo \"Cache: {self.BOLTZ_CACHE_DIR}\"",
            "",
            "boltz predict \\",
            f"    --cache \"{self.BOLTZ_CACHE_DIR}\" \\",
            f"    --out_dir output/ \\",
        ])

        if msa_arg_val:
            script_lines.append(f"    {msa_arg_val} \\")

        script_lines.extend([
            f"    --accelerator {device_val} \\",
            "    --seed 42 \\",
            f"    \"{input_file}\"",
            "",
            "EXIT_CODE=$?",
            "echo \"Boltz2 affinity exit code: $EXIT_CODE\"",
            "",
            "if [ $EXIT_CODE -eq 0 ]; then",
            "    echo \"Boltz2 affinity completed successfully\"",
            "    if [ -d \"output/boltz_results_input\" ]; then",
            "        echo \"Found Boltz2 results directory\"",
            "        # Copy the affinity JSON (Boltz2 names it affinity_<input_stem>.json)",
            "        AFF_FILE=$(find output/boltz_results_input -name 'affinity_*.json' | head -1)",
            "        if [ -f \"$AFF_FILE\" ]; then",
            "            cp \"$AFF_FILE\" output/affinity_scores.json",
            "            echo \"Affinity scores copied from $AFF_FILE\"",
            "        else",
            "            echo \"Warning: no affinity_*.json produced\"",
            "        fi",
            "        # Copy the predicted complex for pose inspection",
            "        CIF_FILE=$(find output/boltz_results_input -name '*_model_0.cif' | head -1)",
            "        if [ -f \"$CIF_FILE\" ]; then",
            "            cp \"$CIF_FILE\" output/predicted_complex.cif",
            "            echo \"Predicted complex copied from $CIF_FILE\"",
            "        fi",
            "    else",
            "        echo \"Warning: Expected Boltz2 output directory not found\"",
            "    fi",
            "else",
            "    echo \"Boltz2 affinity failed with exit code $EXIT_CODE\"",
            "fi",
            "",
            "echo \"Job completed at $(date)\"",
            "exit $EXIT_CODE",
        ])

        return "\n".join(script_lines)

    def _create_boltz2_script(self,
                            job_id: str,
                            job_name: str,
                            input_file: str,
                            remote_dir: str,
                            use_msa: bool,
                            num_models: int) -> str:
        """
        Create SLURM script for Boltz2.
        
        Args:
            job_id: Job ID
            job_name: SLURM job name
            input_file: Input filename (YAML format)
            remote_dir: Remote working directory
            use_msa: Whether to use MSA
            num_models: Number of models to generate
            
        Returns:
            SLURM script content with Unix line endings
        """
        # Determine DEVICE based on the flag
        device_val = "cpu" if self.use_cpu else "gpu"

        # Determine MSA argument string based on the flag
        msa_arg_val = "--use_msa_server" if use_msa else "" # If false, becomes empty string

        # --- Start building the script ---
        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={remote_dir}/boltz2.out",
            f"#SBATCH --error={remote_dir}/boltz2.err",
            f"#SBATCH --partition={self.partition}",
        ]

        # Add account if specified
        if self.slurm_account:
            script_lines.append(f"#SBATCH --account={self.slurm_account}")

        # Add resource requirements based on CPU/GPU mode
        if self.use_cpu:
            script_lines.extend([
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks=4",
                "#SBATCH --cpus-per-task=4",
                "#SBATCH --mem=32G", # Adjust memory as needed for CPU runs
            ])
        else:
            script_lines.extend([
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks-per-node=4", # Adjust if needed
                "#SBATCH --mem=32G", # Adjust memory based on model size/GPU needs
                "#SBATCH --gres=gpu:1", # Request 1 GPU
            ])

        script_lines.extend([
            "#SBATCH --time=02:00:00", # Adjust time as needed
            "",
            "# Load modules",
            f"module use {self.module_path}",
            f"module load {self.boltz2_module}",
            "",
            "# Set working directory",
            f"cd {remote_dir}",
            "",
            "# Create output directory",
            "mkdir -p output",
            "",
            "# Debug: Check input file content",
            "echo \"Input file content:\"",
            f"head -10 {input_file}",
            "echo \"\"",
            "# Run Boltz2",
            "echo \"Starting Boltz2 prediction...\"",
            f"echo \"Input: {input_file}\"",
            "echo \"Output directory: output/\"",
            f"echo \"Number of models: {num_models}\"",
            f"echo \"Use MSA: {use_msa}\"",
            f"echo \"Device: {device_val}\"",
            f"echo \"Cache: {self.BOLTZ_CACHE_DIR}\"",
            "",
            "# Construct the boltz command carefully",
            "boltz predict \\", # Start command
            f"    --cache \"{self.BOLTZ_CACHE_DIR}\" \\", # Cache dir
            f"    --out_dir output/ \\", # Output relative dir
            f"    --diffusion_samples {num_models} \\", # Number of samples
        ])

        # Only add the MSA flag if msa_arg_val is not empty
        if msa_arg_val:
            script_lines.append(f"    {msa_arg_val} \\")

        script_lines.extend([
            f"    --accelerator {device_val} \\", # Device type
            f"    --seed 42 \\", # Seed
            f"    \"{input_file}\"", # Input YAML file (LAST positional argument)
            "",
            "# Capture exit code",
            "EXIT_CODE=$?",
            "echo \"Boltz2 exit code: $EXIT_CODE\"",
            "",
            "# Check for output files regardless of exit code",
            "# The error might be non-fatal since we see 100% progress",
            "echo \"Checking for all generated files...\"",
            "find . -type f -newer $0 2>/dev/null | head -20",
            "echo \"\"",
            "",
            "# Check if successful",
            "if [ $EXIT_CODE -eq 0 ]; then",
            "    echo \"Boltz2 completed successfully\"",
            "    ",
            "    # Debug: List output directory contents",
            "    echo \"Contents of output directory:\"",
            "    ls -la output/ 2>/dev/null || echo \"No output directory\"",
            "    echo \"\"",
            "    ",
            "    # Debug: Find any generated files",
            "    echo \"Looking for any structure files:\"",
            "    find . -name '*.cif' -o -name '*.pdb' 2>/dev/null | head -10",
            "    echo \"\"",
            "    ",
            "    # Check for Boltz2 output structure",
            "    # Boltz2 creates output in boltz_results_<input_name> directory",
            "    if [ -d \"output/boltz_results_input\" ]; then",
            "        echo \"Found Boltz2 results directory\"",
            "        # Look for the CIF file in the predictions subdirectory",
            "        CIF_FILE=\"output/boltz_results_input/predictions/input/input_model_0.cif\"",
            "        if [ -f \"$CIF_FILE\" ]; then",
            "            cp \"$CIF_FILE\" output/predicted_structure.cif",
            "            echo \"Structure copied from $CIF_FILE to output/predicted_structure.cif\"",
            "        else",
            "            # Try to find any CIF file in the boltz results",
            "            STRUCTURE_FILE=$(find output/boltz_results_input -name '*.cif' | head -1)",
            "            if [ -f \"$STRUCTURE_FILE\" ]; then",
            "                cp \"$STRUCTURE_FILE\" output/predicted_structure.cif",
            "                echo \"Structure copied from $STRUCTURE_FILE\"",
            "            fi",
            "        fi",
            "        ",
            "        # Also copy confidence data if available",
            "        CONF_FILE=\"output/boltz_results_input/predictions/input/confidence_input_model_0.json\"",
            "        if [ -f \"$CONF_FILE\" ]; then",
            "            cp \"$CONF_FILE\" output/confidence_scores.json",
            "            echo \"Confidence scores copied to output/confidence_scores.json\"",
            "        fi",
            "    else",
            "        echo \"Warning: Expected Boltz2 output directory not found\"",
            "    fi",
            "else",
            "    echo \"Boltz2 failed with exit code $EXIT_CODE\"",
            "fi",
            "",
            "echo \"Job completed at $(date)\"",
            "",
            "# Ensure Slurm gets the correct exit code",
            "exit $EXIT_CODE",
        ])

        # Join lines with Unix line endings
        script = "\n".join(script_lines)
        return script
    
    def _connect(self):
        """Establish SSH connection to HPC."""
        if self.ssh_client and self.ssh_client.get_transport() and self.ssh_client.get_transport().is_active():
            return
        
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        connect_kwargs = {
            'hostname': self.host,
            'username': self.username,
        }
        
        if self.key_path and os.path.exists(self.key_path):
            connect_kwargs['key_filename'] = self.key_path
            logger.info(f"Connecting to {self.host} with key file")
        elif self.password:
            connect_kwargs['password'] = self.password
            logger.info(f"Connecting to {self.host} with password")
        else:
            raise ValueError("No authentication method available (no key or password)")
        
        self.ssh_client.connect(**connect_kwargs)
        self.sftp_client = self.ssh_client.open_sftp()
        logger.info(f"Successfully connected to {self.host}")
    
    def _disconnect(self):
        """Close SSH connection."""
        if self.sftp_client:
            self.sftp_client.close()
            self.sftp_client = None
        if self.ssh_client:
            self.ssh_client.close()
            self.ssh_client = None
        logger.debug("SSH connection closed")
    
    def _create_remote_directory(self, path: str):
        """Create directory on remote server."""
        stdin, stdout, stderr = self.ssh_client.exec_command(f"mkdir -p {path}")
        stdout.read()
        logger.debug(f"Created remote directory: {path}")
    
    def _upload_content(self, content: str, remote_path: str):
        """Upload content as file to remote server."""
        # Convert Windows line endings to Unix line endings
        content = content.replace('\r\n', '\n')
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt', newline='\n') as f:
            f.write(content)
            temp_path = f.name
        
        try:
            self.sftp_client.put(temp_path, remote_path)
            logger.debug(f"Uploaded content to {remote_path}")
        finally:
            os.unlink(temp_path)
    
    def _download_file(self, remote_path: str) -> str:
        """Download file from remote server."""
        with tempfile.NamedTemporaryFile(mode='r', delete=False) as f:
            temp_path = f.name
        
        try:
            self.sftp_client.get(remote_path, temp_path)
            with open(temp_path, 'r') as f:
                content = f.read()
            return content
        finally:
            os.unlink(temp_path)
    
    def _file_exists(self, remote_path: str) -> bool:
        """Check if file exists on remote server."""
        try:
            self.sftp_client.stat(remote_path)
            return True
        except FileNotFoundError:
            return False
    
    def _submit_job(self, script_path: str) -> str:
        """Submit SLURM job and return job ID."""
        # First, convert the script file to Unix format on the remote server
        convert_cmd = f"dos2unix {script_path} 2>/dev/null || sed -i 's/\r$//' {script_path}"
        stdin, stdout, stderr = self.ssh_client.exec_command(convert_cmd)
        stdout.read()  # Wait for completion
        
        # Now submit the job
        submit_cmd = f"cd {os.path.dirname(script_path)} && sbatch {os.path.basename(script_path)}"
        stdin, stdout, stderr = self.ssh_client.exec_command(submit_cmd)
        output = stdout.read().decode()
        error = stderr.read().decode()
        
        # Check for line ending errors specifically
        if "DOS line breaks" in error:
            logger.error("Script has DOS line breaks. Attempting to fix...")
            # Try to fix and resubmit
            fix_cmd = f"tr -d '\r' < {script_path} > {script_path}.fixed && mv {script_path}.fixed {script_path}"
            stdin, stdout, stderr = self.ssh_client.exec_command(fix_cmd)
            stdout.read()
            
            # Retry submission
            stdin, stdout, stderr = self.ssh_client.exec_command(submit_cmd)
            output = stdout.read().decode()
            error = stderr.read().decode()
        
        if error and "Submitted batch job" not in error:
            logger.warning(f"SLURM submission warning: {error}")
        
        # Parse SLURM job ID from output or error (sometimes it's in stderr)
        import re
        combined_output = output + error
        match = re.search(r'Submitted batch job (\d+)', combined_output)
        if match:
            logger.info(f"SLURM job submitted successfully: {match.group(1)}")
            return match.group(1)
        else:
            raise RuntimeError(f"Failed to submit job. Output: {output} Error: {error}")


# Singleton instance
_hpc_boltz2_client = None

def get_hpc_boltz2_client() -> HPCBoltz2Client:
    """Get or create HPC Boltz2 client singleton."""
    global _hpc_boltz2_client
    if _hpc_boltz2_client is None:
        _hpc_boltz2_client = HPCBoltz2Client()
    return _hpc_boltz2_client
