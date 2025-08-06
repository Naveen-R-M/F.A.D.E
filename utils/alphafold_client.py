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
        
        self.logger.debug(f"get_best_model_path: Received output_dir='{output_dir}', protein_name='{protein_name}'")

        # First, try to find the subdirectory created by AlphaFold3 based on the input JSON name
        # This is typically the basename of the input JSON file (e.g., KRAS_af3_input)
        # We need to glob for this as the exact name might vary slightly or be unknown.
        potential_json_output_dirs = glob.glob(os.path.join(output_dir, f"{protein_name}_af3_input*"))
        self.logger.debug(f"get_best_model_path: Found potential_json_output_dirs: {potential_json_output_dirs}")
        
        effective_output_dir = None
        if potential_json_output_dirs:
            # AlphaFold3 creates a timestamped subdirectory within the protein_name_af3_input directory
            # We need to find the latest one.
            base_af3_input_dir = potential_json_output_dirs[0]
            self.logger.debug(f"get_best_model_path: Selected base_af3_input_dir: {base_af3_input_dir}")
            
            # Find all timestamped directories within this base directory, or the direct protein_name.lower() directory
            timestamped_or_direct_dirs = glob.glob(os.path.join(base_af3_input_dir, f"{protein_name.lower()}_*"))
            # Also check for the direct protein_name.lower() directory
            direct_protein_dir = os.path.join(base_af3_input_dir, protein_name.lower())
            if os.path.isdir(direct_protein_dir):
                timestamped_or_direct_dirs.append(direct_protein_dir)

            self.logger.debug(f"get_best_model_path: Found timestamped_or_direct_dirs: {timestamped_or_direct_dirs}")
            
            if timestamped_or_direct_dirs:
                # Sort by modification time (newest first) to get the latest run
                timestamped_or_direct_dirs.sort(key=os.path.getmtime, reverse=True)
                effective_output_dir = timestamped_or_direct_dirs[0]
                self.logger.debug(f"get_best_model_path: Selected latest effective_output_dir: {effective_output_dir}")
            else:
                # Fallback if no timestamped or direct dirs found, maybe it's directly in base_af3_input_dir
                effective_output_dir = base_af3_input_dir
                self.logger.debug(f"get_best_model_path: No specific subdirs, falling back to base_af3_input_dir: {effective_output_dir}")
        else:
            # If no specific JSON output dir found, assume output_dir itself might contain the protein_name subdir
            # This part remains largely the same as a fallback, though less likely for AF3
            self.logger.debug(f"get_best_model_path: No potential_json_output_dirs found, checking direct subdirs.")
            if os.path.isdir(os.path.join(output_dir, protein_name.lower())):
                effective_output_dir = os.path.join(output_dir, protein_name.lower())
                self.logger.debug(f"get_best_model_path: Found protein_name.lower() subdir: {effective_output_dir}")
            elif os.path.isdir(os.path.join(output_dir, protein_name)):
                effective_output_dir = os.path.join(output_dir, protein_name)
                self.logger.debug(f"get_best_model_path: Found protein_name subdir: {effective_output_dir}")
            else:
                effective_output_dir = output_dir # Fallback to original output_dir
                self.logger.debug(f"get_best_model_path: Falling back to original output_dir: {effective_output_dir}")

        if not effective_output_dir or not os.path.exists(effective_output_dir):
            self.logger.warning(f"Effective AlphaFold output directory not found for {protein_name}: {effective_output_dir}")
            return None
            
        self.logger.info(f"Searching for best model in effective directory: {effective_output_dir}")

        # Try to find the best model based on ranking JSON (AlphaFold2 format)
        ranking_path = os.path.join(effective_output_dir, "ranking_debug.json")
        
        if os.path.exists(ranking_path):
            try:
                with open(ranking_path, "r") as f:
                    ranking_data = json.load(f)
                    
                if "plddts" in ranking_data:
                    best_model_name = max(ranking_data["plddts"].items(), key=lambda x: x[1])[0]
                    best_model_path = os.path.join(effective_output_dir, f"unrelaxed_{best_model_name}.pdb")
                    if os.path.exists(best_model_path):
                        self.logger.info(f"Found best model (AF2 ranking): {best_model_path}")
                        return best_model_path
            except Exception as e:
                self.logger.error(f"Error parsing AF2 ranking file {ranking_path}: {e}")
        
        # Try to find the best model based on summary confidences (AlphaFold3 format)
        summary_confidences_path = os.path.join(effective_output_dir, f"{protein_name.lower()}_summary_confidences.json")
        
        if os.path.exists(summary_confidences_path):
            try:
                with open(summary_confidences_path, "r") as f:
                    summary_data = json.load(f)
                
                # AlphaFold3's kras_summary_confidences.json has a "ranking_score"
                if "ranking_score" in summary_data:
                    # For AF3, the main model is usually kras_model.cif
                    model_cif_path = os.path.join(effective_output_dir, f"{protein_name.lower()}_model.cif")
                    if os.path.exists(model_cif_path):
                        self.logger.info(f"Found best model (AF3 summary): {model_cif_path}")
                        return model_cif_path
            except Exception as e:
                self.logger.error(f"Error parsing AF3 summary confidences file {summary_confidences_path}: {e}")

        # Fallback: Check for any .cif or .pdb files and return the newest one
        cif_files = glob.glob(os.path.join(effective_output_dir, "*.cif"))
        pdb_files = glob.glob(os.path.join(effective_output_dir, "*.pdb"))
        
        all_model_files = cif_files + pdb_files
        
        if all_model_files:
            all_model_files.sort(key=os.path.getmtime, reverse=True)
            best_model_path = all_model_files[0]
            self.logger.info(f"Found model (fallback, newest .cif/.pdb): {best_model_path}")
            return best_model_path
            
        self.logger.warning(f"No model files (.cif or .pdb) found in effective directory: {effective_output_dir}")
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
        self.logger.debug(f"get_confidence_scores: Received output_dir='{output_dir}', protein_name='{protein_name}'")
        potential_json_output_dirs = glob.glob(os.path.join(output_dir, f"{protein_name}_af3_input*"))
        self.logger.debug(f"get_confidence_scores: Found potential_json_output_dirs: {potential_json_output_dirs}")
        
        effective_output_dir = None
        if potential_json_output_dirs:
            base_af3_input_dir = potential_json_output_dirs[0]
            self.logger.debug(f"get_confidence_scores: Selected base_af3_input_dir: {base_af3_input_dir}")
            timestamped_or_direct_dirs = glob.glob(os.path.join(base_af3_input_dir, f"{protein_name.lower()}_*"))
            direct_protein_dir = os.path.join(base_af3_input_dir, protein_name.lower())
            if os.path.isdir(direct_protein_dir):
                timestamped_or_direct_dirs.append(direct_protein_dir)
            self.logger.debug(f"get_confidence_scores: Found timestamped_or_direct_dirs: {timestamped_or_direct_dirs}")
            if timestamped_or_direct_dirs:
                timestamped_or_direct_dirs.sort(key=os.path.getmtime, reverse=True)
                effective_output_dir = timestamped_or_direct_dirs[0]
                self.logger.debug(f"get_confidence_scores: Selected latest effective_output_dir: {effective_output_dir}")
            else:
                effective_output_dir = base_af3_input_dir
                self.logger.debug(f"get_confidence_scores: No specific subdirs, falling back to base_af3_input_dir: {effective_output_dir}")
        else:
            self.logger.debug(f"get_confidence_scores: No potential_json_output_dirs found, checking direct subdirs.")
            if os.path.isdir(os.path.join(output_dir, protein_name.lower())):
                effective_output_dir = os.path.join(output_dir, protein_name.lower())
                self.logger.debug(f"get_confidence_scores: Found protein_name.lower() subdir: {effective_output_dir}")
            elif os.path.isdir(os.path.join(output_dir, protein_name)):
                effective_output_dir = os.path.join(output_dir, protein_name)
                self.logger.debug(f"get_confidence_scores: Found protein_name subdir: {effective_output_dir}")
            else:
                effective_output_dir = output_dir
                self.logger.debug(f"get_confidence_scores: Falling back to original output_dir: {effective_output_dir}")

        if not effective_output_dir or not os.path.exists(effective_output_dir):
            self.logger.warning(f"Effective AlphaFold output directory not found for {protein_name}: {effective_output_dir}")
            return None
            
        self.logger.info(f"Searching for confidence scores in effective directory: {effective_output_dir}")

        # Try to find the confidence scores from summary confidences (AlphaFold3 format)
        summary_confidences_path = os.path.join(effective_output_dir, f"{protein_name.lower()}_summary_confidences.json")
        confidences_path = os.path.join(effective_output_dir, f"{protein_name.lower()}_confidences.json")
        
        confidence_scores = {}

        if os.path.exists(summary_confidences_path):
            try:
                with open(summary_confidences_path, "r") as f:
                    summary_data = json.load(f)
                
                if "ranking_score" in summary_data:
                    confidence_scores["overall"] = summary_data["ranking_score"]
                if "ptm" in summary_data:
                    confidence_scores["ptm"] = summary_data["ptm"]
            except Exception as e:
                self.logger.error(f"Error parsing AF3 summary confidences file {summary_confidences_path}: {e}")
        
        if os.path.exists(confidences_path):
            try:
                with open(confidences_path, "r") as f:
                    confidences_data = json.load(f)
                
                if "atom_plddts" in confidences_data and confidences_data["atom_plddts"]:
                    confidence_scores["avg_plddt"] = sum(confidences_data["atom_plddts"]) / len(confidences_data["atom_plddts"])
            except Exception as e:
                self.logger.error(f"Error parsing AF3 confidences file {confidences_path}: {e}")

        if confidence_scores:
            return confidence_scores
        
        # Fallback to AlphaFold2 ranking_debug.json if AF3 files not found or parsed
        ranking_path = os.path.join(effective_output_dir, "ranking_debug.json")
        if os.path.exists(ranking_path):
            try:
                with open(ranking_path, "r") as f:
                    ranking_data = json.load(f)
                    
                if "plddts" in ranking_data:
                    confidence_scores["plddts"] = ranking_data["plddts"]
                    if ranking_data["plddts"]:
                        confidence_scores["avg_plddt"] = sum(ranking_data["plddts"].values()) / len(ranking_data["plddts"])
                
                if "ptms" in ranking_data:
                    confidence_scores["ptms"] = ranking_data["ptms"]
                    if ranking_data["ptms"]:
                        confidence_scores["avg_ptm"] = sum(ranking_data["ptms"].values()) / len(ranking_data["ptms"])
                
                if confidence_scores:
                    return confidence_scores
            except Exception as e:
                self.logger.error(f"Error parsing AF2 ranking file {ranking_path}: {e}")
        
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
        self.logger.debug(f"get_confidence_scores: Received output_dir='{output_dir}', protein_name='{protein_name}'")
        potential_json_output_dirs = glob.glob(os.path.join(output_dir, f"{protein_name}_af3_input*"))
        self.logger.debug(f"get_confidence_scores: Found potential_json_output_dirs: {potential_json_output_dirs}")
        
        effective_output_dir = None
        if potential_json_output_dirs:
            base_af3_input_dir = potential_json_output_dirs[0]
            self.logger.debug(f"get_confidence_scores: Selected base_af3_input_dir: {base_af3_input_dir}")
            timestamped_or_direct_dirs = glob.glob(os.path.join(base_af3_input_dir, f"{protein_name.lower()}_*"))
            direct_protein_dir = os.path.join(base_af3_input_dir, protein_name.lower())
            if os.path.isdir(direct_protein_dir):
                timestamped_or_direct_dirs.append(direct_protein_dir)
            self.logger.debug(f"get_confidence_scores: Found timestamped_or_direct_dirs: {timestamped_or_direct_dirs}")
            if timestamped_or_direct_dirs:
                timestamped_or_direct_dirs.sort(key=os.path.getmtime, reverse=True)
                effective_output_dir = timestamped_or_direct_dirs[0]
                self.logger.debug(f"get_confidence_scores: Selected latest effective_output_dir: {effective_output_dir}")
            else:
                effective_output_dir = base_af3_input_dir
                self.logger.debug(f"get_confidence_scores: No specific subdirs, falling back to base_af3_input_dir: {effective_output_dir}")
        else:
            self.logger.debug(f"get_confidence_scores: No potential_json_output_dirs found, checking direct subdirs.")
            if os.path.isdir(os.path.join(output_dir, protein_name.lower())):
                effective_output_dir = os.path.join(output_dir, protein_name.lower())
                self.logger.debug(f"get_confidence_scores: Found protein_name.lower() subdir: {effective_output_dir}")
            elif os.path.isdir(os.path.join(output_dir, protein_name)):
                effective_output_dir = os.path.join(output_dir, protein_name)
                self.logger.debug(f"get_confidence_scores: Found protein_name subdir: {effective_output_dir}")
            else:
                effective_output_dir = output_dir
                self.logger.debug(f"get_confidence_scores: Falling back to original output_dir: {effective_output_dir}")

        if not effective_output_dir or not os.path.exists(effective_output_dir):
            self.logger.warning(f"Effective AlphaFold output directory not found for {protein_name}: {effective_output_dir}")
            return None
            
        self.logger.info(f"Searching for multimer info in effective directory: {effective_output_dir}")

        # Try to find the multimer info from summary confidences (AlphaFold3 format)
        summary_confidences_path = os.path.join(effective_output_dir, f"{protein_name.lower()}_summary_confidences.json")
        
        multimer_info = {}

        if os.path.exists(summary_confidences_path):
            try:
                with open(summary_confidences_path, "r") as f:
                    summary_data = json.load(f)
                
                if "iptm" in summary_data and summary_data["iptm"] is not None:
                    multimer_info["iptm"] = summary_data["iptm"]
                if "chain_pair_iptm" in summary_data:
                    multimer_info["chain_pair_iptm"] = summary_data["chain_pair_iptm"]
                
                if multimer_info:
                    return multimer_info
            except Exception as e:
                self.logger.error(f"Error parsing AF3 summary confidences file {summary_confidences_path}: {e}")

        # Fallback to AlphaFold2 ranking_debug.json if AF3 files not found or parsed
        ranking_path = os.path.join(effective_output_dir, "ranking_debug.json")
        if os.path.exists(ranking_path):
            try:
                with open(ranking_path, "r") as f:
                    ranking_data = json.load(f)
                    
                if "iptm+ptm" in ranking_data:
                    multimer_info["iptm_ptm"] = ranking_data["iptm+ptm"]
                
                if "interface" in ranking_data:
                    multimer_info["interface"] = ranking_data["interface"]
                
                if multimer_info:
                    return multimer_info
            except Exception as e:
                self.logger.error(f"Error parsing ranking file {ranking_path}: {e}")
        
        self.logger.warning(f"No multimer information found in effective directory: {effective_output_dir}")
        return None
