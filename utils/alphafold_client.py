"""
AlphaFold Client for F.A.D.E

This module provides an interface to interact with AlphaFold3 for
protein structure prediction.
"""

import os
import json
import glob
from typing import Any, Dict, List, Optional, Tuple, Union

from utils.logging import get_logger


class AlphaFoldClient:
    """
    Client for interacting with AlphaFold.
    """
    
    def __init__(self) -> None:
        """Initialize the AlphaFold client."""
        self.logger = get_logger("fade.utils.alphafold_client")
    
    def get_best_model_path(self, output_dir: str, protein_name: str) -> Optional[str]:
        """
        Get the path to the best model in the AlphaFold output directory.
        
        Args:
            output_dir: Base path to the AlphaFold output directory (e.g., .../structures/KRAS).
            protein_name: The name of the protein (e.g., "KRAS").
            
        Returns:
            Path to the best model PDB file, or None if not found.
        """
        # AlphaFold3 often creates a subdirectory named after the input JSON basename,
        # and then another subdirectory named after the protein itself.
        # Example: .../structures/KRAS/KRAS_af3_input/kras/
        
        # First, try to find the subdirectory created by AlphaFold3 based on the input JSON name
        # This is typically the basename of the input JSON file (e.g., KRAS_af3_input)
        # We need to glob for this as the exact name might vary slightly or be unknown.
        potential_json_output_dirs = glob.glob(os.path.join(output_dir, f"{protein_name}_af3_input*"))
        
        effective_output_dir = None
        if potential_json_output_dirs:
            # Assuming the first match is the correct one, then look for the protein_name subdirectory
            temp_dir = potential_json_output_dirs[0]
            if os.path.isdir(os.path.join(temp_dir, protein_name.lower())):
                effective_output_dir = os.path.join(temp_dir, protein_name.lower())
            elif os.path.isdir(os.path.join(temp_dir, protein_name)):
                effective_output_dir = os.path.join(temp_dir, protein_name)
            else:
                # Fallback to just the first potential_json_output_dir if protein_name subdir not found
                effective_output_dir = temp_dir
        else:
            # If no specific JSON output dir found, assume output_dir itself might contain the protein_name subdir
            if os.path.isdir(os.path.join(output_dir, protein_name.lower())):
                effective_output_dir = os.path.join(output_dir, protein_name.lower())
            elif os.path.isdir(os.path.join(output_dir, protein_name)):
                effective_output_dir = os.path.join(output_dir, protein_name)
            else:
                effective_output_dir = output_dir # Fallback to original output_dir

        if not effective_output_dir or not os.path.exists(effective_output_dir):
            self.logger.warning(f"Effective AlphaFold output directory not found for {protein_name}: {effective_output_dir}")
            return None
            
        self.logger.info(f"Searching for best model in effective directory: {effective_output_dir}")

        # Try to find the best model based on ranking JSON
        ranking_path = os.path.join(effective_output_dir, "ranking_debug.json")
        
        if os.path.exists(ranking_path):
            try:
                with open(ranking_path, "r") as f:
                    ranking_data = json.load(f)
                    
                # Get the model name with the highest pLDDT score
                if "plddts" in ranking_data:
                    best_model = max(ranking_data["plddts"].items(), key=lambda x: x[1])[0]
                    
                    # Convert to PDB file path
                    best_model_path = os.path.join(effective_output_dir, f"unrelaxed_{best_model}.pdb")
                    
                    if os.path.exists(best_model_path):
                        self.logger.info(f"Found best model: {best_model_path}")
                        return best_model_path
            except Exception as e:
                self.logger.error(f"Error parsing ranking file {ranking_path}: {e}")
        
        # Fallback: Check for PDB files and return the one with the highest pLDDT
        # (based on model_* naming convention)
        pdb_files = glob.glob(os.path.join(effective_output_dir, "unrelaxed_model_*.pdb"))
        
        if not pdb_files:
            # Try looking for relaxed models
            pdb_files = glob.glob(os.path.join(effective_output_dir, "relaxed_model_*.pdb"))
            
        if not pdb_files:
            # Look for any PDB files
            pdb_files = glob.glob(os.path.join(effective_output_dir, "*.pdb"))
            
        if pdb_files:
            # Sort by modification time (newest first) and return the first one as a heuristic
            pdb_files.sort(key=os.path.getmtime, reverse=True)
            best_model_path = pdb_files[0]
            self.logger.info(f"Found model (fallback): {best_model_path}")
            return best_model_path
            
        # No PDB files found
        self.logger.warning(f"No PDB files found in effective directory: {effective_output_dir}")
        return None
    
    def get_confidence_scores(self, output_dir: str, protein_name: str) -> Optional[Dict[str, Any]]:
        """
        Get confidence scores from the AlphaFold output directory.
        
        Args:
            output_dir: Base path to the AlphaFold output directory.
            protein_name: The name of the protein.
            
        Returns:
            Dictionary containing confidence scores, or None if not found.
        """
        # Determine the effective output directory, similar to get_best_model_path
        potential_json_output_dirs = glob.glob(os.path.join(output_dir, f"{protein_name}_af3_input*"))
        
        effective_output_dir = None
        if potential_json_output_dirs:
            temp_dir = potential_json_output_dirs[0]
            if os.path.isdir(os.path.join(temp_dir, protein_name.lower())):
                effective_output_dir = os.path.join(temp_dir, protein_name.lower())
            elif os.path.isdir(os.path.join(temp_dir, protein_name)):
                effective_output_dir = os.path.join(temp_dir, protein_name)
            else:
                effective_output_dir = temp_dir
        else:
            if os.path.isdir(os.path.join(output_dir, protein_name.lower())):
                effective_output_dir = os.path.join(output_dir, protein_name.lower())
            elif os.path.isdir(os.path.join(output_dir, protein_name)):
                effective_output_dir = os.path.join(output_dir, protein_name)
            else:
                effective_output_dir = output_dir

        if not effective_output_dir or not os.path.exists(effective_output_dir):
            self.logger.warning(f"Effective AlphaFold output directory not found for {protein_name}: {effective_output_dir}")
            return None
            
        self.logger.info(f"Searching for confidence scores in effective directory: {effective_output_dir}")

        # Try to find the confidence scores from ranking JSON
        ranking_path = os.path.join(effective_output_dir, "ranking_debug.json")
        
        if os.path.exists(ranking_path):
            try:
                with open(ranking_path, "r") as f:
                    ranking_data = json.load(f)
                    
                # Extract confidence scores
                confidence_scores = {}
                
                if "plddts" in ranking_data:
                    confidence_scores["plddts"] = ranking_data["plddts"]
                    
                    # Calculate average pLDDT
                    if ranking_data["plddts"]:
                        confidence_scores["avg_plddt"] = sum(ranking_data["plddts"].values()) / len(ranking_data["plddts"])
                
                if "ptms" in ranking_data:
                    confidence_scores["ptms"] = ranking_data["ptms"]
                    
                    # Calculate average pTM
                    if ranking_data["ptms"]:
                        confidence_scores["avg_ptm"] = sum(ranking_data["ptms"].values()) / len(ranking_data["ptms"])
                
                if confidence_scores:
                    return confidence_scores
            except Exception as e:
                self.logger.error(f"Error parsing ranking file {ranking_path}: {e}")
        
        # Fallback: Return None if no confidence scores found
        self.logger.warning(f"No confidence scores found in effective directory: {effective_output_dir}")
        return None
    
    def get_multimer_info(self, output_dir: str, protein_name: str) -> Optional[Dict[str, Any]]:
        """
        Get multimer information from the AlphaFold output directory.
        
        Args:
            output_dir: Base path to the AlphaFold output directory.
            protein_name: The name of the protein.
            
        Returns:
            Dictionary containing multimer information, or None if not found.
        """
        # Determine the effective output directory, similar to get_best_model_path
        potential_json_output_dirs = glob.glob(os.path.join(output_dir, f"{protein_name}_af3_input*"))
        
        effective_output_dir = None
        if potential_json_output_dirs:
            temp_dir = potential_json_output_dirs[0]
            if os.path.isdir(os.path.join(temp_dir, protein_name.lower())):
                effective_output_dir = os.path.join(temp_dir, protein_name.lower())
            elif os.path.isdir(os.path.join(temp_dir, protein_name)):
                effective_output_dir = os.path.join(temp_dir, protein_name)
            else:
                effective_output_dir = temp_dir
        else:
            if os.path.isdir(os.path.join(output_dir, protein_name.lower())):
                effective_output_dir = os.path.join(output_dir, protein_name.lower())
            elif os.path.isdir(os.path.join(output_dir, protein_name)):
                effective_output_dir = os.path.join(output_dir, protein_name)
            else:
                effective_output_dir = output_dir

        if not effective_output_dir or not os.path.exists(effective_output_dir):
            self.logger.warning(f"Effective AlphaFold output directory not found for {protein_name}: {effective_output_dir}")
            return None
            
        self.logger.info(f"Searching for multimer info in effective directory: {effective_output_dir}")

        # Try to find the multimer info from ranking JSON
        ranking_path = os.path.join(effective_output_dir, "ranking_debug.json")
        
        if os.path.exists(ranking_path):
            try:
                with open(ranking_path, "r") as f:
                    ranking_data = json.load(f)
                    
                # Extract multimer information
                multimer_info = {}
                
                if "iptm+ptm" in ranking_data:
                    multimer_info["iptm_ptm"] = ranking_data["iptm+ptm"]
                
                if "interface" in ranking_data:
                    multimer_info["interface"] = ranking_data["interface"]
                
                if multimer_info:
                    return multimer_info
            except Exception as e:
                self.logger.error(f"Error parsing ranking file {ranking_path}: {e}")
        
        # Fallback: Return None if no multimer info found
        self.logger.warning(f"No multimer information found in effective directory: {effective_output_dir}")
        return None
