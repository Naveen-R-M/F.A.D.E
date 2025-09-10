"""
Sequence validator component for the Target Selector agent.

This component validates protein sequences for scientific accuracy.
"""

from typing import Any, Dict, List, Optional, Tuple, Set
import json
import os
import re
from utils.gemini_client import GeminiClient

class SequenceValidator:
    """
    Validates protein sequences for scientific accuracy.
    """
    
    def __init__(
        self, 
        llm_client: Optional[GeminiClient] = None,
        reference_dir: Optional[str] = None
    ) -> None:
        """
        Initialize the sequence validator.
        
        Args:
            llm_client: Client for the LLM (e.g., GeminiClient)
            reference_dir: Directory containing reference data
        """
        self.llm_client = llm_client
        
        # Set up reference directory
        if reference_dir is None:
            # Determine project root directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            agents_dir = os.path.dirname(current_dir)
            project_root = os.path.dirname(agents_dir)
            reference_dir = os.path.join(project_root, "data", "references", "protein_sequences")
        
        self.reference_dir = reference_dir
        
        # Load reference data if available
        self.reference_data = self._load_reference_data()
        
        # Define standard validation checks
        self.standard_checks = {
            'valid_amino_acids': self._check_valid_amino_acids,
            'sequence_length': self._check_sequence_length,
            'key_residues': self._check_key_residues,
            'motifs': self._check_motifs
        }
    
    def validate(
        self, 
        protein_name: str, 
        sequence: str,
        metadata: Optional[Dict[str, Any]] = None
    ) -> Tuple[bool, str, Dict[str, Any]]:
        """
        Validate if a sequence matches expected properties of the named protein.
        
        Args:
            protein_name: Name of the protein
            sequence: Protein sequence to validate
            metadata: Additional metadata about the protein
            
        Returns:
            Tuple of (is_valid, reason, details)
        """
        if metadata is None:
            metadata = {}
        
        # Initialize validation results
        validation_results = {
            'protein_name': protein_name,
            'sequence_length': len(sequence),
            'checks': {},
            'overall_score': 0.0,
            'warnings': [],
            'errors': []
        }
        
        # Normalize protein name for lookup
        normalized_name = self._normalize_protein_name(protein_name)
        
        # Get reference data for this protein if available
        reference = self.reference_data.get(normalized_name)
        
        # Run standard validation checks
        for check_name, check_func in self.standard_checks.items():
            result = check_func(sequence, reference, protein_name, metadata)
            validation_results['checks'][check_name] = result
            
            if not result['pass'] and result['severity'] == 'error':
                validation_results['errors'].append(f"{check_name}: {result['message']}")
            elif not result['pass'] and result['severity'] == 'warning':
                validation_results['warnings'].append(f"{check_name}: {result['message']}")
        
        # Calculate overall score (0.0 to 1.0)
        if validation_results['checks']:
            passed_checks = sum(1 for r in validation_results['checks'].values() if r['pass'])
            total_checks = len(validation_results['checks'])
            validation_results['overall_score'] = passed_checks / total_checks
        
        # Use LLM for advanced validation if available and reference data is missing
        if self.llm_client and not reference:
            llm_validation = self._validate_with_llm(protein_name, sequence, metadata)
            validation_results['llm_validation'] = llm_validation
            
            # Incorporate LLM findings
            if not llm_validation['is_valid']:
                for issue in llm_validation['issues']:
                    if issue['severity'] == 'error':
                        validation_results['errors'].append(f"LLM: {issue['description']}")
                    else:
                        validation_results['warnings'].append(f"LLM: {issue['description']}")
            
            # Adjust overall score
            if 'confidence_score' in llm_validation:
                # Weight LLM score at 40% if we have other checks
                if validation_results['checks']:
                    validation_results['overall_score'] = (
                        0.6 * validation_results['overall_score'] + 
                        0.4 * llm_validation['confidence_score']
                    )
                else:
                    validation_results['overall_score'] = llm_validation['confidence_score']
        
        # Determine final validity
        is_valid = (
            len(validation_results['errors']) == 0 and 
            validation_results['overall_score'] >= 0.7
        )
        
        # Create reason message
        if is_valid:
            if validation_results['warnings']:
                reason = f"Valid with warnings: {validation_results['warnings'][0]}"
            else:
                reason = "Sequence passed all validation checks"
        else:
            if validation_results['errors']:
                reason = validation_results['errors'][0]
            elif validation_results['overall_score'] < 0.7:
                reason = f"Low confidence score: {validation_results['overall_score']:.2f}"
            else:
                reason = "Failed validation for unknown reasons"
        
        return is_valid, reason, validation_results
    
    def _load_reference_data(self) -> Dict[str, Dict[str, Any]]:
        """
        Load reference data for protein sequences.
        
        Returns:
            Dictionary mapping protein names to reference data
        """
        reference_data = {}
        
        # Check if reference directory exists
        if not os.path.exists(self.reference_dir):
            os.makedirs(self.reference_dir, exist_ok=True)
            return reference_data
        
        # Load JSON files from reference directory
        for filename in os.listdir(self.reference_dir):
            if filename.endswith('.json'):
                try:
                    with open(os.path.join(self.reference_dir, filename), 'r') as f:
                        protein_data = json.load(f)
                        
                        # Add to reference data using normalized names
                        if 'name' in protein_data:
                            normalized_name = self._normalize_protein_name(protein_data['name'])
                            reference_data[normalized_name] = protein_data
                            
                            # Add alternative names if available
                            for alt_name in protein_data.get('alternative_names', []):
                                normalized_alt = self._normalize_protein_name(alt_name)
                                reference_data[normalized_alt] = protein_data
                except Exception as e:
                    # Skip invalid files
                    continue
        
        return reference_data
    
    def _normalize_protein_name(self, name: str) -> str:
        """
        Normalize protein name for consistent lookup.
        
        Args:
            name: Original protein name
            
        Returns:
            Normalized name
        """
        # Convert to lowercase
        normalized = name.lower()
        
        # Remove common prefixes
        prefixes = ['human ', 'mouse ', 'rat ', 'h-', 'm-', 'r-']
        for prefix in prefixes:
            if normalized.startswith(prefix):
                normalized = normalized[len(prefix):]
                break
        
        # Remove whitespace and special characters
        normalized = re.sub(r'[^a-z0-9]', '', normalized)
        
        return normalized
    
    def _check_valid_amino_acids(
        self, 
        sequence: str, 
        reference: Optional[Dict[str, Any]],
        protein_name: str,
        metadata: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Check if sequence contains only valid amino acid codes.
        
        Args:
            sequence: Protein sequence
            reference: Reference data for the protein
            protein_name: Name of the protein
            metadata: Additional metadata
            
        Returns:
            Validation result
        """
        # Valid amino acid codes (including ambiguous ones)
        valid_codes = set('ACDEFGHIKLMNPQRSTVWY')
        extended_codes = set('ACDEFGHIKLMNPQRSTVWYBZX*-')
        
        # Check if all characters in sequence are valid amino acids
        invalid_chars = [c for c in sequence.upper() if c not in extended_codes]
        
        if invalid_chars:
            return {
                'pass': False,
                'message': f"Sequence contains invalid amino acid codes: {''.join(set(invalid_chars))}",
                'severity': 'error',
                'details': {
                    'invalid_chars': list(set(invalid_chars)),
                    'positions': [i for i, c in enumerate(sequence.upper()) if c not in extended_codes]
                }
            }
        
        # Warn about ambiguous or extended codes
        ambiguous_chars = [c for c in sequence.upper() if c not in valid_codes]
        if ambiguous_chars:
            return {
                'pass': True,
                'message': f"Sequence contains ambiguous amino acid codes: {''.join(set(ambiguous_chars))}",
                'severity': 'warning',
                'details': {
                    'ambiguous_chars': list(set(ambiguous_chars)),
                    'positions': [i for i, c in enumerate(sequence.upper()) if c not in valid_codes]
                }
            }
        
        return {
            'pass': True,
            'message': "Sequence contains only valid amino acid codes",
            'severity': 'none',
            'details': {}
        }
    
    def _check_sequence_length(
        self, 
        sequence: str, 
        reference: Optional[Dict[str, Any]],
        protein_name: str,
        metadata: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Check if sequence length is within expected range.
        
        Args:
            sequence: Protein sequence
            reference: Reference data for the protein
            protein_name: Name of the protein
            metadata: Additional metadata
            
        Returns:
            Validation result
        """
        # If we have reference data with length information, use it
        if reference and 'length' in reference:
            expected_length = reference['length']
            tolerance = reference.get('length_tolerance', 5)
            
            min_length = expected_length - tolerance
            max_length = expected_length + tolerance
            
            if min_length <= len(sequence) <= max_length:
                return {
                    'pass': True,
                    'message': f"Sequence length ({len(sequence)}) is within expected range ({min_length}-{max_length})",
                    'severity': 'none',
                    'details': {
                        'length': len(sequence),
                        'expected_length': expected_length,
                        'tolerance': tolerance
                    }
                }
            else:
                return {
                    'pass': False,
                    'message': f"Sequence length ({len(sequence)}) is outside expected range ({min_length}-{max_length})",
                    'severity': 'warning',
                    'details': {
                        'length': len(sequence),
                        'expected_length': expected_length,
                        'tolerance': tolerance,
                        'difference': abs(len(sequence) - expected_length)
                    }
                }
        
        # If no reference data, use heuristics based on protein name
        # Known proteins with typical lengths
        protein_length_heuristics = {
            'kras': (188, 190),  # KRAS is typically around 189 amino acids
            'braf': (765, 775),  # BRAF is typically around 766 amino acids
            'egfr': (1200, 1220),  # EGFR is typically around 1210 amino acids
            'her2': (1250, 1260),  # HER2 is typically around 1255 amino acids
            'p53': (390, 400),  # p53 is typically around 393 amino acids
            'pten': (400, 410)  # PTEN is typically around 403 amino acids
        }
        
        # Check if protein name matches any of our known proteins
        for known_protein, (min_len, max_len) in protein_length_heuristics.items():
            if known_protein in self._normalize_protein_name(protein_name):
                if min_len <= len(sequence) <= max_len:
                    return {
                        'pass': True,
                        'message': f"Sequence length ({len(sequence)}) is within expected range for {known_protein.upper()} ({min_len}-{max_len})",
                        'severity': 'none',
                        'details': {
                            'length': len(sequence),
                            'expected_range': (min_len, max_len),
                            'protein_match': known_protein.upper()
                        }
                    }
                else:
                    return {
                        'pass': False,
                        'message': f"Sequence length ({len(sequence)}) is outside expected range for {known_protein.upper()} ({min_len}-{max_len})",
                        'severity': 'warning',
                        'details': {
                            'length': len(sequence),
                            'expected_range': (min_len, max_len),
                            'protein_match': known_protein.upper(),
                            'difference': min(abs(len(sequence) - min_len), abs(len(sequence) - max_len))
                        }
                    }
        
        # General validation for unknown proteins
        if len(sequence) < 30:
            return {
                'pass': False,
                'message': f"Sequence is unusually short ({len(sequence)} aa), may be incomplete",
                'severity': 'warning',
                'details': {'length': len(sequence)}
            }
        elif len(sequence) > 3000:
            return {
                'pass': False,
                'message': f"Sequence is unusually long ({len(sequence)} aa), may include non-protein regions",
                'severity': 'warning',
                'details': {'length': len(sequence)}
            }
        
        return {
            'pass': True,
            'message': f"Sequence length ({len(sequence)} aa) is within reasonable range for a protein",
            'severity': 'none',
            'details': {'length': len(sequence)}
        }
    
    def _check_key_residues(
        self, 
        sequence: str, 
        reference: Optional[Dict[str, Any]],
        protein_name: str,
        metadata: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Check if key residues are present at expected positions.
        
        Args:
            sequence: Protein sequence
            reference: Reference data for the protein
            protein_name: Name of the protein
            metadata: Additional metadata
            
        Returns:
            Validation result
        """
        # If reference data with key residues is available, use it
        if reference and 'key_residues' in reference:
            key_residues = reference['key_residues']
            
            # Check each key residue
            missing_residues = []
            for residue_info in key_residues:
                position = residue_info.get('position')
                expected = residue_info.get('residue')
                
                # Skip if position or expected residue is not specified
                if position is None or not expected:
                    continue
                
                # Adjust for 0-based indexing
                position_idx = position - 1
                
                # Check if position is within sequence
                if position_idx < 0 or position_idx >= len(sequence):
                    missing_residues.append({
                        'position': position,
                        'expected': expected,
                        'issue': 'Position out of range'
                    })
                    continue
                
                # Check if residue matches
                actual = sequence[position_idx].upper()
                if actual != expected.upper():
                    missing_residues.append({
                        'position': position,
                        'expected': expected,
                        'actual': actual,
                        'issue': 'Residue mismatch'
                    })
            
            if missing_residues:
                return {
                    'pass': False,
                    'message': f"Sequence is missing {len(missing_residues)} key residues",
                    'severity': 'error' if len(missing_residues) > len(key_residues) // 2 else 'warning',
                    'details': {
                        'missing_residues': missing_residues,
                        'total_key_residues': len(key_residues)
                    }
                }
            
            return {
                'pass': True,
                'message': f"Sequence contains all {len(key_residues)} key residues",
                'severity': 'none',
                'details': {
                    'total_key_residues': len(key_residues)
                }
            }
        
        # Known key residues for common proteins
        protein_key_residues = {
            'kras': [
                {'position': 12, 'residue': 'G', 'note': 'Common mutation site G12D in cancer'},
                {'position': 13, 'residue': 'G', 'note': 'Adjacent to G12 mutation site'},
                {'position': 61, 'residue': 'Q', 'note': 'Common mutation site Q61 in cancer'}
            ],
            'braf': [
                {'position': 600, 'residue': 'V', 'note': 'Common mutation site V600E in cancer'}
            ],
            'egfr': [
                {'position': 790, 'residue': 'T', 'note': 'Common mutation site T790M in cancer'}
            ]
        }
        
        # Check if protein name matches any of our known proteins with key residues
        for known_protein, residues in protein_key_residues.items():
            if known_protein in self._normalize_protein_name(protein_name):
                # Check each key residue
                missing_residues = []
                for residue_info in residues:
                    position = residue_info.get('position')
                    expected = residue_info.get('residue')
                    
                    # Adjust for 0-based indexing
                    position_idx = position - 1
                    
                    # Check if position is within sequence
                    if position_idx < 0 or position_idx >= len(sequence):
                        missing_residues.append({
                            'position': position,
                            'expected': expected,
                            'issue': 'Position out of range'
                        })
                        continue
                    
                    # Check if residue matches
                    actual = sequence[position_idx].upper()
                    if actual != expected.upper():
                        missing_residues.append({
                            'position': position,
                            'expected': expected,
                            'actual': actual,
                            'issue': 'Residue mismatch'
                        })
                
                if missing_residues:
                    return {
                        'pass': False,
                        'message': f"Sequence is missing {len(missing_residues)} key residues for {known_protein.upper()}",
                        'severity': 'error' if len(missing_residues) > len(residues) // 2 else 'warning',
                        'details': {
                            'missing_residues': missing_residues,
                            'total_key_residues': len(residues),
                            'protein_match': known_protein.upper()
                        }
                    }
                
                return {
                    'pass': True,
                    'message': f"Sequence contains all {len(residues)} key residues for {known_protein.upper()}",
                    'severity': 'none',
                    'details': {
                        'total_key_residues': len(residues),
                        'protein_match': known_protein.upper()
                    }
                }
        
        # Check for mutations mentioned in metadata
        if 'mutations' in metadata:
            mutations = metadata.get('mutations', [])
            
            # Verify each mutation
            invalid_mutations = []
            for mutation in mutations:
                original = mutation.get('original_residue')
                position = mutation.get('position')
                mutated = mutation.get('mutated_residue')
                
                # Skip if any required information is missing
                if not original or not position or not mutated:
                    continue
                
                # Adjust for 0-based indexing
                position_idx = position - 1
                
                # Check if position is within sequence
                if position_idx < 0 or position_idx >= len(sequence):
                    invalid_mutations.append({
                        'mutation': f"{original}{position}{mutated}",
                        'issue': 'Position out of range'
                    })
                    continue
                
                # Check if the current residue matches the expected mutated residue
                actual = sequence[position_idx].upper()
                if actual != mutated.upper():
                    invalid_mutations.append({
                        'mutation': f"{original}{position}{mutated}",
                        'actual': actual,
                        'issue': 'Mutated residue not found'
                    })
            
            if invalid_mutations:
                return {
                    'pass': False,
                    'message': f"Sequence does not contain {len(invalid_mutations)} expected mutations",
                    'severity': 'warning',
                    'details': {
                        'invalid_mutations': invalid_mutations,
                        'total_mutations': len(mutations)
                    }
                }
        
        # Default to passing for unknown proteins with no specific checks
        return {
            'pass': True,
            'message': "No key residues defined for this protein",
            'severity': 'none',
            'details': {}
        }
    
    def _check_motifs(
        self, 
        sequence: str, 
        reference: Optional[Dict[str, Any]],
        protein_name: str,
        metadata: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Check if expected motifs are present in the sequence.
        
        Args:
            sequence: Protein sequence
            reference: Reference data for the protein
            protein_name: Name of the protein
            metadata: Additional metadata
            
        Returns:
            Validation result
        """
        # If reference data with motifs is available, use it
        if reference and 'motifs' in reference:
            motifs = reference['motifs']
            
            # Check each motif
            missing_motifs = []
            for motif_info in motifs:
                pattern = motif_info.get('pattern')
                description = motif_info.get('description', 'Unknown motif')
                
                # Skip if pattern is not specified
                if not pattern:
                    continue
                
                # Check if motif is present
                if not re.search(pattern, sequence):
                    missing_motifs.append({
                        'pattern': pattern,
                        'description': description
                    })
            
            if missing_motifs:
                return {
                    'pass': False,
                    'message': f"Sequence is missing {len(missing_motifs)} expected motifs",
                    'severity': 'warning',
                    'details': {
                        'missing_motifs': missing_motifs,
                        'total_motifs': len(motifs)
                    }
                }
            
            return {
                'pass': True,
                'message': f"Sequence contains all {len(motifs)} expected motifs",
                'severity': 'none',
                'details': {
                    'total_motifs': len(motifs)
                }
            }
        
        # Known motifs for common proteins
        protein_motifs = {
            'kras': [
                {'pattern': 'GAGGVGKS', 'description': 'P-loop (G1) motif'},
                {'pattern': 'DTAGQ', 'description': 'G3 motif'}
            ],
            'kinase': [
                {'pattern': 'G.G..G.V', 'description': 'ATP-binding glycine-rich loop'},
                {'pattern': 'HRD.K', 'description': 'Catalytic loop'}
            ]
        }
        
        # Check if protein name matches any of our known proteins with motifs
        for known_protein, motifs in protein_motifs.items():
            if known_protein in self._normalize_protein_name(protein_name):
                # Check each motif
                missing_motifs = []
                for motif_info in motifs:
                    pattern = motif_info.get('pattern')
                    description = motif_info.get('description', 'Unknown motif')
                    
                    # Check if motif is present
                    if not re.search(pattern, sequence):
                        missing_motifs.append({
                            'pattern': pattern,
                            'description': description
                        })
                
                if missing_motifs:
                    return {
                        'pass': False,
                        'message': f"Sequence is missing {len(missing_motifs)} expected motifs for {known_protein.upper()}",
                        'severity': 'warning',
                        'details': {
                            'missing_motifs': missing_motifs,
                            'total_motifs': len(motifs),
                            'protein_match': known_protein.upper()
                        }
                    }
                
                return {
                    'pass': True,
                    'message': f"Sequence contains all {len(motifs)} expected motifs for {known_protein.upper()}",
                    'severity': 'none',
                    'details': {
                        'total_motifs': len(motifs),
                        'protein_match': known_protein.upper()
                    }
                }
        
        # Common protein family motifs
        if 'kinase' in self._normalize_protein_name(protein_name):
            # Check for common kinase motifs
            kinase_motifs = [
                {'pattern': 'G.G..G.V', 'description': 'ATP-binding glycine-rich loop'},
                {'pattern': 'HRD.K', 'description': 'Catalytic loop'}
            ]
            
            missing_motifs = []
            for motif_info in kinase_motifs:
                pattern = motif_info.get('pattern')
                description = motif_info.get('description')
                
                if not re.search(pattern, sequence):
                    missing_motifs.append({
                        'pattern': pattern,
                        'description': description
                    })
            
            if missing_motifs:
                return {
                    'pass': False,
                    'message': f"Sequence is missing {len(missing_motifs)} common kinase motifs",
                    'severity': 'warning',
                    'details': {
                        'missing_motifs': missing_motifs,
                        'total_motifs': len(kinase_motifs)
                    }
                }
            
            return {
                'pass': True,
                'message': "Sequence contains common kinase motifs",
                'severity': 'none',
                'details': {
                    'total_motifs': len(kinase_motifs)
                }
            }
        
        # Default to passing for unknown proteins with no specific motifs
        return {
            'pass': True,
            'message': "No specific motifs defined for this protein",
            'severity': 'none',
            'details': {}
        }
    
    def _validate_with_llm(
        self, 
        protein_name: str, 
        sequence: str,
        metadata: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Use LLM to validate the protein sequence.
        
        Args:
            protein_name: Name of the protein
            sequence: Protein sequence
            metadata: Additional metadata
            
        Returns:
            LLM validation results
        """
        # Truncate sequence for prompt if too long
        display_sequence = sequence
        if len(sequence) > 500:
            display_sequence = sequence[:200] + "..." + sequence[-200:]
        
        # Format metadata
        metadata_str = json.dumps(metadata, indent=2)
        
        # Construct the prompt
        prompt = f"""
        As a protein sequence validator, please analyze this protein sequence for {protein_name}.
        
        Protein: {protein_name}
        
        Sequence (first few and last few amino acids shown if long):
        {display_sequence}
        
        Sequence length: {len(sequence)} amino acids
        
        Additional metadata:
        {metadata_str}
        
        Please analyze whether this sequence appears to be a valid sequence for {protein_name}.
        Consider:
        1. Typical length of {protein_name}
        2. Key residues or motifs expected for {protein_name}
        3. Any inconsistencies with known properties of {protein_name}
        4. Whether any mutations mentioned in the metadata are correctly reflected
        
        Provide your response as a JSON object with the following structure:
        {{
            "is_valid": true/false,
            "confidence_score": 0.8, // 0.0 to 1.0 indicating confidence
            "issues": [
                {{
                    "description": "Issue description",
                    "severity": "error/warning"
                }}
            ],
            "key_observations": [
                "Observation 1",
                "Observation 2"
            ],
            "suggestions": [
                "Suggestion 1",
                "Suggestion 2"
            ]
        }}
        
        IMPORTANT: Ensure the JSON is valid and complete. If you don't have enough information, indicate lower confidence.
        """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.3)
            
            # Extract JSON from response
            json_start = response.find('{')
            json_end = response.rfind('}') + 1
            
            if json_start >= 0 and json_end > json_start:
                json_str = response[json_start:json_end]
                validation_result = json.loads(json_str)
                
                # Ensure required fields are present
                if 'is_valid' not in validation_result:
                    validation_result['is_valid'] = True
                if 'confidence_score' not in validation_result:
                    validation_result['confidence_score'] = 0.5
                if 'issues' not in validation_result:
                    validation_result['issues'] = []
                if 'key_observations' not in validation_result:
                    validation_result['key_observations'] = []
                if 'suggestions' not in validation_result:
                    validation_result['suggestions'] = []
                
                return validation_result
            else:
                # Default validation result if JSON extraction fails
                return {
                    'is_valid': True,
                    'confidence_score': 0.5,
                    'issues': [],
                    'key_observations': ['Failed to parse LLM response'],
                    'suggestions': []
                }
        
        except Exception as e:
            # Default validation result if LLM fails
            return {
                'is_valid': True,
                'confidence_score': 0.4,
                'issues': [],
                'key_observations': [f'Error using LLM for validation: {str(e)}'],
                'suggestions': []
            }
    
    def add_reference_protein(
        self, 
        name: str, 
        sequence: str, 
        length: Optional[int] = None,
        key_residues: Optional[List[Dict[str, Any]]] = None,
        motifs: Optional[List[Dict[str, Any]]] = None,
        alternative_names: Optional[List[str]] = None,
        organism: Optional[str] = None
    ) -> bool:
        """
        Add a reference protein to the database.
        
        Args:
            name: Protein name
            sequence: Protein sequence
            length: Expected length (defaults to len(sequence))
            key_residues: List of key residues
            motifs: List of motifs
            alternative_names: Alternative names for the protein
            organism: Organism the protein is from
            
        Returns:
            True if protein was added successfully, False otherwise
        """
        try:
            # Create reference directory if it doesn't exist
            if not os.path.exists(self.reference_dir):
                os.makedirs(self.reference_dir, exist_ok=True)
            
            # Create reference data
            reference_data = {
                'name': name,
                'sequence': sequence,
                'length': length or len(sequence),
                'length_tolerance': 5,
                'key_residues': key_residues or [],
                'motifs': motifs or [],
                'alternative_names': alternative_names or [],
                'organism': organism
            }
            
            # Save to file
            normalized_name = self._normalize_protein_name(name)
            file_path = os.path.join(self.reference_dir, f"{normalized_name}.json")
            
            with open(file_path, 'w') as f:
                json.dump(reference_data, f, indent=2)
            
            # Add to in-memory reference data
            self.reference_data[normalized_name] = reference_data
            
            # Add alternative names
            if alternative_names:
                for alt_name in alternative_names:
                    normalized_alt = self._normalize_protein_name(alt_name)
                    self.reference_data[normalized_alt] = reference_data
            
            return True
            
        except Exception as e:
            # Failed to add reference protein
            return False
