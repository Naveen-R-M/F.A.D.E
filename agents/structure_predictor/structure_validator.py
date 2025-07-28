"""
Structure Validator for F.A.D.E

This module validates protein structures, checking for issues like
missing residues, poor geometry, and low confidence regions.
"""

import os
from typing import Any, Dict, List, Optional, Tuple, Union
import json
import numpy as np

from Bio.PDB import PDBParser, Structure
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

# Suppress Bio.PDB warnings about missing atoms
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


class StructureValidator:
    """
    Class for validating protein structures.
    """
    
    def __init__(self, llm_client=None) -> None:
        """
        Initialize the structure validator.
        
        Args:
            llm_client: Optional client for LLM integration.
        """
        self.parser = PDBParser(QUIET=True)
        self.llm_client = llm_client
    
    def validate(
        self,
        pdb_file: str,
        sequence: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Validate a protein structure.
        
        Args:
            pdb_file: Path to the PDB file.
            sequence: Optional original sequence for comparison.
            
        Returns:
            Dictionary containing validation results.
        """
        # Parse the PDB file
        structure = self.parser.get_structure("structure", pdb_file)
        
        # Run validation checks
        missing_residues = self._check_missing_residues(structure, sequence)
        geometry_issues = self._check_geometry(structure)
        confidence_scores = self._extract_confidence_scores(structure)
        
        # Calculate overall score
        validation_scores = {
            "missing_residues_score": 1.0 - (len(missing_residues) / len(sequence) if sequence else 0),
            "geometry_score": 1.0 - (len(geometry_issues) / 100),  # Normalize to 0-1
            "plddt": confidence_scores.get("confidence", 0.0)
        }
        
        # PTM (predicted TM-score) is often available in the AlphaFold output
        # For now, we'll estimate it based on pLDDT
        ptm_estimate = min(1.0, confidence_scores.get("confidence", 0.0) * 1.2)
        validation_scores["ptm"] = ptm_estimate
        
        # Overall score is weighted average
        overall_score = (
            validation_scores["missing_residues_score"] * 0.3 +
            validation_scores["geometry_score"] * 0.3 +
            validation_scores["plddt"] * 0.4
        )
        validation_scores["overall_score"] = overall_score
        
        # Add interpretations
        validation_interpretations = self._get_interpretations(validation_scores)
        
        # Create result dictionary
        result = {
            "pdb_file": pdb_file,
            "missing_residues": missing_residues,
            "geometry_issues": geometry_issues,
            "confidence_scores": confidence_scores,
            "validation_scores": validation_scores,
            "interpretations": validation_interpretations,
            "is_valid": overall_score >= 0.7,  # Consider valid if score is at least 0.7
            "overall_score": overall_score
        }
        
        # Optionally use LLM for deeper analysis
        if self.llm_client and overall_score < 0.9:
            llm_analysis = self._analyze_with_llm(result)
            result["llm_analysis"] = llm_analysis
        
        return result
    
    def _check_missing_residues(
        self,
        structure: Structure,
        sequence: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Check for missing residues in the structure.
        
        Args:
            structure: BioPython Structure object.
            sequence: Optional original sequence for comparison.
            
        Returns:
            List of dictionaries describing missing residues.
        """
        missing_residues = []
        
        # If no sequence is provided, we can only check for gaps in numbering
        if not sequence:
            # Check for gaps in residue numbering
            for chain in structure.get_chains():
                residue_ids = [residue.get_id()[1] for residue in chain if residue.get_id()[0] == " "]
                
                if not residue_ids:
                    continue
                    
                # Check for gaps
                min_id = min(residue_ids)
                max_id = max(residue_ids)
                
                for i in range(min_id, max_id + 1):
                    if i not in residue_ids:
                        missing_residues.append({
                            "chain_id": chain.get_id(),
                            "residue_id": i,
                            "residue_name": "UNK",  # Unknown residue type
                            "description": f"Gap in residue numbering at position {i}"
                        })
        else:
            # More sophisticated check using provided sequence
            # For each chain in the structure
            for chain in structure.get_chains():
                # Extract residue sequence from structure
                struct_sequence = ""
                residue_mapping = {}
                
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        res_id = residue.get_id()[1]
                        res_name = residue.get_resname()
                        res_code = self._three_to_one(res_name)
                        struct_sequence += res_code
                        residue_mapping[len(struct_sequence) - 1] = res_id
                
                # Simple sequence alignment to find missing residues
                # Note: For production, use a proper sequence alignment algorithm
                if len(struct_sequence) < len(sequence):
                    # Structure is missing residues
                    if struct_sequence in sequence:
                        # Simple case: structure is a continuous substring
                        start_idx = sequence.find(struct_sequence)
                        
                        # Report missing residues at the beginning
                        for i in range(start_idx):
                            missing_residues.append({
                                "chain_id": chain.get_id(),
                                "residue_id": i + 1,  # 1-indexed
                                "residue_name": sequence[i],
                                "description": f"Missing residue at beginning of chain"
                            })
                            
                        # Report missing residues at the end
                        for i in range(start_idx + len(struct_sequence), len(sequence)):
                            missing_residues.append({
                                "chain_id": chain.get_id(),
                                "residue_id": i + 1,  # 1-indexed
                                "residue_name": sequence[i],
                                "description": f"Missing residue at end of chain"
                            })
                    else:
                        # More complex case: fragmented alignment
                        # For this example, just report the difference in length
                        missing_residues.append({
                            "chain_id": chain.get_id(),
                            "residue_id": 0,
                            "residue_name": "UNK",
                            "description": f"Structure is missing {len(sequence) - len(struct_sequence)} residues"
                        })
        
        return missing_residues
    
    def _check_geometry(self, structure: Structure) -> List[Dict[str, Any]]:
        """
        Check for geometry issues in the structure.
        
        Args:
            structure: BioPython Structure object.
            
        Returns:
            List of dictionaries describing geometry issues.
        """
        geometry_issues = []
        
        # Check bond lengths and angles
        # This is a placeholder - in a real implementation, you would use
        # more sophisticated geometry validation
        
        return geometry_issues
    
    def _extract_confidence_scores(self, structure: Structure) -> Dict[str, float]:
        """
        Extract confidence scores from B-factors (for AlphaFold models).
        
        Args:
            structure: BioPython Structure object.
            
        Returns:
            Dictionary containing confidence scores.
        """
        # Get all B-factors
        b_factors = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":  # Skip hetero residues
                        for atom in residue:
                            b_factors.append(atom.get_bfactor())
        
        # Calculate statistics
        if b_factors:
            mean_b = np.mean(b_factors)
            max_b = np.max(b_factors)
            min_b = np.min(b_factors)
            
            # For AlphaFold, B-factors represent pLDDT scores (0-100)
            # Higher is better
            confidence = mean_b / 100.0
            
            return {
                "mean_b_factor": float(mean_b),
                "max_b_factor": float(max_b),
                "min_b_factor": float(min_b),
                "confidence": float(confidence)
            }
        
        return {
            "mean_b_factor": 0.0,
            "max_b_factor": 0.0,
            "min_b_factor": 0.0,
            "confidence": 0.0
        }
    
    def _get_interpretations(self, scores: Dict[str, float]) -> Dict[str, str]:
        """
        Get human-readable interpretations of validation scores.
        
        Args:
            scores: Dictionary of validation scores.
            
        Returns:
            Dictionary mapping score names to interpretations.
        """
        interpretations = {}
        
        # Interpret missing residues score
        missing_residues_score = scores.get("missing_residues_score", 0.0)
        if missing_residues_score >= 0.95:
            interpretations["missing_residues"] = "Excellent: The structure contains almost all residues."
        elif missing_residues_score >= 0.9:
            interpretations["missing_residues"] = "Good: The structure is missing a few residues."
        elif missing_residues_score >= 0.8:
            interpretations["missing_residues"] = "Fair: The structure is missing several residues."
        else:
            interpretations["missing_residues"] = "Poor: The structure is missing many residues."
        
        # Interpret geometry score
        geometry_score = scores.get("geometry_score", 0.0)
        if geometry_score >= 0.95:
            interpretations["geometry"] = "Excellent: The structure has very good geometry."
        elif geometry_score >= 0.9:
            interpretations["geometry"] = "Good: The structure has good geometry with minor issues."
        elif geometry_score >= 0.8:
            interpretations["geometry"] = "Fair: The structure has several geometry issues."
        else:
            interpretations["geometry"] = "Poor: The structure has many geometry issues."
        
        # Interpret pLDDT score
        plddt = scores.get("plddt", 0.0)
        if plddt >= 0.9:
            interpretations["plddt"] = "Very high confidence (pLDDT > 90)."
        elif plddt >= 0.7:
            interpretations["plddt"] = "High confidence (pLDDT > 70)."
        elif plddt >= 0.5:
            interpretations["plddt"] = "Medium confidence (pLDDT > 50)."
        else:
            interpretations["plddt"] = "Low confidence (pLDDT < 50)."
        
        # Interpret overall score
        overall_score = scores.get("overall_score", 0.0)
        if overall_score >= 0.9:
            interpretations["overall"] = "Excellent: This structure is very reliable for further analysis."
        elif overall_score >= 0.8:
            interpretations["overall"] = "Good: This structure is reliable for most analyses."
        elif overall_score >= 0.7:
            interpretations["overall"] = "Fair: This structure may be suitable for analysis, but caution is advised."
        elif overall_score >= 0.5:
            interpretations["overall"] = "Poor: This structure has significant issues and should be used with caution."
        else:
            interpretations["overall"] = "Very poor: This structure is not reliable for further analysis."
        
        return interpretations
    
    def _analyze_with_llm(self, validation_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Use LLM to analyze validation results in more depth.
        
        Args:
            validation_results: Dictionary of validation results.
            
        Returns:
            Dictionary containing LLM analysis.
        """
        if not self.llm_client:
            return {"error": "LLM client not available"}
            
        # Prepare prompt
        prompt = f"""
        Analyze the following protein structure validation results:
        
        PDB File: {validation_results.get("pdb_file", "Unknown")}
        Overall Score: {validation_results.get("overall_score", 0.0):.2f}
        pLDDT Score: {validation_results.get("validation_scores", {}).get("plddt", 0.0):.2f}
        Missing Residues Score: {validation_results.get("validation_scores", {}).get("missing_residues_score", 0.0):.2f}
        Geometry Score: {validation_results.get("validation_scores", {}).get("geometry_score", 0.0):.2f}
        
        Missing Residues: {len(validation_results.get("missing_residues", []))}
        Geometry Issues: {len(validation_results.get("geometry_issues", []))}
        
        Please analyze these validation results and provide:
        1. An assessment of the overall quality of the structure
        2. Recommendations for how to improve the structure (if needed)
        3. Whether the structure is suitable for drug discovery applications
        4. Key areas of concern (if any)
        
        Format your response as a JSON with the following structure:
        {{
            "quality_assessment": "Your assessment here",
            "recommendations": ["Recommendation 1", "Recommendation 2", ...],
            "suitability": "Your assessment of suitability for drug discovery",
            "concerns": ["Concern 1", "Concern 2", ...],
            "confidence": 0.8  // Your confidence in this analysis (0.0 to 1.0)
        }}
        """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.2)
            
            # Extract JSON from response
            try:
                # Find JSON in the response
                json_start = response.find('{')
                json_end = response.rfind('}') + 1
                if json_start >= 0 and json_end > json_start:
                    json_str = response[json_start:json_end]
                    analysis = json.loads(json_str)
                    return analysis
                else:
                    return {"error": "Could not find JSON in LLM response", "raw_response": response}
            except json.JSONDecodeError:
                return {"error": "Could not parse JSON from LLM response", "raw_response": response}
        except Exception as e:
            return {"error": f"Error calling LLM: {e}"}
    
    def _three_to_one(self, three_letter_code: str) -> str:
        """
        Convert three-letter amino acid code to one-letter code.
        
        Args:
            three_letter_code: Three-letter amino acid code.
            
        Returns:
            One-letter amino acid code.
        """
        # Conversion dictionary
        three_to_one = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "SEC": "U", "PYL": "O", "MSE": "M"  # Selenocysteine, Pyrrolysine, Selenomethionine
        }
        
        return three_to_one.get(three_letter_code, "X")
