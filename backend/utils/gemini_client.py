"""
Gemini API client for F.A.D.E
"""

import os
import json
import logging
from typing import Any, Dict, List, Optional, Union
import re

import google.generativeai as genai
from dotenv import load_dotenv


class GeminiClient:
    """
    Client for interacting with the Google Gemini API.
    """
    
    def __init__(self, api_key: Optional[str] = None, model: Optional[str] = None) -> None:
        """
        Initialize the Gemini client.
        
        Args:
            api_key: Google Gemini API key. If not provided, it will be loaded from
                environment variables.
            model: Gemini model to use. If not provided, it will be loaded from
                environment variables.
        """
        # Load environment variables from .env file
        load_dotenv()
        
        # Set up logger
        self.logger = logging.getLogger("fade.gemini")
        
        # Use provided API key or load from environment variables
        self.api_key = api_key or os.getenv("GEMINI_API_KEY")
        
        if not self.api_key:
            raise ValueError(
                "Gemini API key not provided. Set GEMINI_API_KEY environment variable or "
                "provide it directly."
            )
        
        # Configure the Gemini API
        genai.configure(api_key=self.api_key)
        
        # Get available models
        self.models = genai.list_models()
        
        # Use provided model or load from environment variables, with fallback
        self.default_model = model or os.getenv("GEMINI_MODEL") or "models/gemini-1.5-flash"
        
    def get_available_models(self) -> List[str]:
        """
        Get a list of available Gemini models.
        
        Returns:
            List of model names.
        """
        return [model.name for model in self.models]
    
    def generate_text(
        self, 
        prompt: str, 
        model: str = None,
        temperature: float = 0.7,
        max_tokens: int = 1024,
        top_p: float = 0.95,
        top_k: int = 40,
    ) -> str:
        """
        Generate text using the Gemini API.
        
        Args:
            prompt: The input text prompt.
            model: The model to use. Defaults to the default model.
            temperature: Controls randomness of output. Higher values (e.g., 0.8) make
                output more random, lower values (e.g., 0.2) make it more deterministic.
            max_tokens: Maximum number of tokens to generate.
            top_p: Nucleus sampling parameter. Keep the most likely tokens whose 
                cumulative probability exceeds this value.
            top_k: Only sample from the top K options for each next token.
            
        Returns:
            The generated text response.
        """
        model_name = model or self.default_model
        generation_config = {
            "temperature": temperature,
            "top_p": top_p,
            "top_k": top_k,
            "max_output_tokens": max_tokens,
        }
        
        # Add retry logic for API requests
        max_retries = 3
        attempt = 0
        
        while attempt < max_retries:
            try:
                model = genai.GenerativeModel(model_name)
                response = model.generate_content(
                    prompt,
                    generation_config=generation_config
                )
                return response.text
            except Exception as e:
                attempt += 1
                self.logger.warning(f"API request failed (attempt {attempt}/{max_retries}): {e}")
                if attempt >= max_retries:
                    raise
        
        # This should not be reached due to the raise in the loop
        return ""
    
    def extract_structured_data(
        self, 
        prompt: str, 
        schema: Dict[str, Any],
        model: str = None,
        temperature: float = 0.2,  # Lower temperature for more deterministic results
        max_attempts: int = 3,  # Maximum number of attempts to get valid JSON
    ) -> Dict[str, Any]:
        """
        Extract structured data from text using the Gemini API.
        
        Args:
            prompt: The input text prompt.
            schema: JSON schema for the structured data to extract.
            model: The model to use. Defaults to the default model.
            temperature: Controls randomness of output. Lower value for more deterministic
                results.
            max_attempts: Maximum number of attempts to get valid JSON.
                
        Returns:
            Structured data according to the provided schema.
        """
        # Create a stronger prompt that emphasizes the need for complete, valid JSON
        extraction_prompt = f"""
        Extract structured information from the following text according to the schema below.
        
        SCHEMA:
        {json.dumps(schema, indent=2)}
        
        TEXT TO ANALYZE:
        {prompt}
        
        CRITICAL INSTRUCTIONS:
        - Extract all relevant information that matches the schema
        - Provide output as a complete, valid JSON object following the schema exactly
        - If information is not present, use null values
        - DO NOT include any explanations or text outside the JSON object
        - Ensure the output is a COMPLETE and VALID JSON object with all closing brackets
        - Check that all strings are properly quoted and terminated
        - VERIFY that your response is a complete JSON object before returning
        
        JSON OUTPUT:
        """
        
        model_name = model or self.default_model
        
        # Try multiple attempts with decreasing temperature
        for attempt in range(max_attempts):
            current_temp = max(0.01, temperature - (attempt * 0.05))
            self.logger.info(f"Attempt {attempt+1}/{max_attempts} with temperature {current_temp}")
            
            try:
                # Generate the structured data response
                response_text = self.generate_text(
                    extraction_prompt,
                    model=model_name,
                    temperature=current_temp,
                    max_tokens=2048  # Increase max tokens to ensure we get complete response
                )
                
                # Clean up the response text to ensure it's valid JSON
                cleaned_response = response_text.strip()
                
                # Remove markdown code block syntax if present
                if cleaned_response.startswith("```json"):
                    cleaned_response = cleaned_response[7:]
                elif cleaned_response.startswith("```"):
                    cleaned_response = cleaned_response[3:]
                
                if cleaned_response.endswith("```"):
                    cleaned_response = cleaned_response[:-3]
                
                # Attempt to fix common JSON issues
                cleaned_response = self._fix_json(cleaned_response)
                
                # Try to parse the JSON
                structured_data = json.loads(cleaned_response.strip())
                return structured_data
                
            except json.JSONDecodeError as e:
                self.logger.warning(f"JSON parsing failed (attempt {attempt+1}/{max_attempts}): {e}")
                
                # If this is the last attempt, try to recreate a valid JSON
                if attempt == max_attempts - 1:
                    try:
                        # Attempt to recreate a valid JSON structure from partial response
                        fixed_json = self._reconstruct_json(cleaned_response, schema)
                        if fixed_json:
                            self.logger.info("Successfully reconstructed JSON from partial response")
                            return fixed_json
                    except Exception as fix_error:
                        self.logger.error(f"Failed to reconstruct JSON: {fix_error}")
                        
                    # If all attempts fail, raise a detailed error
                    raise ValueError(
                        f"Failed to parse JSON response after {max_attempts} attempts: {e}\n"
                        f"Response: {response_text}"
                    )
            
            except Exception as e:
                self.logger.error(f"Unexpected error: {e}")
                if attempt == max_attempts - 1:
                    raise
    
    def _fix_json(self, json_str: str) -> str:
        """
        Attempt to fix common JSON formatting issues.
        
        Args:
            json_str: JSON string to fix
            
        Returns:
            Fixed JSON string
        """
        # Remove any trailing commas before closing brackets or braces
        json_str = re.sub(r',\s*}', '}', json_str)
        json_str = re.sub(r',\s*]', ']', json_str)
        
        # Check for unclosed JSON object
        if json_str.count('{') > json_str.count('}'):
            json_str += '}' * (json_str.count('{') - json_str.count('}'))
        
        # Check for unclosed JSON array
        if json_str.count('[') > json_str.count(']'):
            json_str += ']' * (json_str.count('[') - json_str.count(']'))
        
        return json_str
    
    def _reconstruct_json(self, partial_json: str, schema: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Attempt to reconstruct a valid JSON object from a partial response.
        
        Args:
            partial_json: Partial JSON string
            schema: JSON schema to use for reconstruction
            
        Returns:
            Reconstructed JSON object or None if reconstruction failed
        """
        # Simple case: if we can fix it with closing brackets, do that
        try:
            fixed_json = self._fix_json(partial_json)
            return json.loads(fixed_json)
        except json.JSONDecodeError:
            pass
        
        # More complex case: try to extract the valid parts
        # Extract any valid key-value pairs we can find
        try:
            # Extract protein_targets if present
            protein_targets = []
            protein_match = re.search(r'"protein_targets"\s*:\s*(\[.*?\])', partial_json, re.DOTALL)
            if protein_match:
                try:
                    protein_json = protein_match.group(1)
                    # Fix unterminated arrays
                    if protein_json.count('[') > protein_json.count(']'):
                        protein_json += ']' * (protein_json.count('[') - protein_json.count(']'))
                    protein_targets = json.loads(protein_json)
                except:
                    # If we can't parse the entire array, try to get at least the first item
                    try:
                        first_item_match = re.search(r'\[\s*({.*?})', protein_json, re.DOTALL)
                        if first_item_match:
                            first_item = first_item_match.group(1)
                            # Ensure it's properly terminated
                            if first_item.count('{') > first_item.count('}'):
                                first_item += '}' * (first_item.count('{') - first_item.count('}'))
                            protein_targets = [json.loads(first_item)]
                    except:
                        pass
            
            # Extract molecule_properties if present
            molecule_properties = {}
            properties_match = re.search(r'"molecule_properties"\s*:\s*({.*?})', partial_json, re.DOTALL)
            if properties_match:
                try:
                    properties_json = properties_match.group(1)
                    # Fix unterminated objects
                    if properties_json.count('{') > properties_json.count('}'):
                        properties_json += '}' * (properties_json.count('{') - properties_json.count('}'))
                    molecule_properties = json.loads(properties_json)
                except:
                    pass
            
            # Extract confidence if present
            confidence = 0.5  # Default
            confidence_match = re.search(r'"confidence"\s*:\s*([\d.]+)', partial_json)
            if confidence_match:
                try:
                    confidence = float(confidence_match.group(1))
                except:
                    pass
            
            # Reconstruct a minimal valid JSON
            reconstructed = {
                "protein_targets": protein_targets,
                "molecule_properties": molecule_properties,
                "confidence": confidence,
                "clarification_questions": []
            }
            
            return reconstructed
        except Exception as e:
            self.logger.error(f"Failed to reconstruct JSON: {e}")
            return None
    
    def extract_protein_info(self, query: str) -> Dict[str, Any]:
        """
        Extract protein target information from a natural language query.
        
        Args:
            query: The natural language query about protein targets.
            
        Returns:
            Structured information about protein targets and requirements.
        """
        schema = {
            "protein_targets": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "name": {"type": "string"},
                        "full_name": {"type": "string"},
                        "organism": {"type": "string"},
                        "mutations": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "original_residue": {"type": "string"},
                                    "position": {"type": "integer"},
                                    "mutated_residue": {"type": "string"}
                                }
                            }
                        },
                        "binding_sites": {
                            "type": "array",
                            "items": {"type": "string"}
                        }
                    },
                    "required": ["name"]
                }
            },
            "molecule_properties": {
                "type": "object",
                "properties": {
                    "blood_brain_barrier_permeability": {"type": "boolean"},
                    "lipinski_rule_of_five": {"type": "boolean"},
                    "toxicity_requirements": {"type": "string"},
                    "solubility_requirements": {"type": "string"},
                    "other_requirements": {"type": "string"}
                }
            },
            "confidence": {
                "type": "number",
                "description": "Confidence in the extraction (0.0 to 1.0)"
            },
            "clarification_questions": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Questions to ask for clarification if information is ambiguous"
            }
        }
        
        return self.extract_structured_data(query, schema, temperature=0.1)
