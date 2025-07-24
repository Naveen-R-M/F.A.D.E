"""
Gemini API client for F.A.D.E
"""

import os
import json
from typing import Any, Dict, List, Optional, Union

import google.generativeai as genai
from dotenv import load_dotenv


class GeminiClient:
    """
    Client for interacting with the Google Gemini API.
    """
    
    def __init__(self, api_key: Optional[str] = None) -> None:
        """
        Initialize the Gemini client.
        
        Args:
            api_key: Google Gemini API key. If not provided, it will be loaded from
                environment variables.
        """
        # Load environment variables from .env file
        load_dotenv()
        
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
        # Using a Flash model which has higher quotas
        self.default_model = "models/gemini-1.5-flash"
        
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
        
        model = genai.GenerativeModel(model_name)
        response = model.generate_content(
            prompt,
            generation_config=generation_config
        )
        
        return response.text
    
    def extract_structured_data(
        self, 
        prompt: str, 
        schema: Dict[str, Any],
        model: str = None,
        temperature: float = 0.2,  # Lower temperature for more deterministic results
    ) -> Dict[str, Any]:
        """
        Extract structured data from text using the Gemini API.
        
        Args:
            prompt: The input text prompt.
            schema: JSON schema for the structured data to extract.
            model: The model to use. Defaults to the default model.
            temperature: Controls randomness of output. Lower value for more deterministic
                results.
                
        Returns:
            Structured data according to the provided schema.
        """
        extraction_prompt = f"""
        Extract structured information from the following text according to the schema below.
        
        SCHEMA:
        {json.dumps(schema, indent=2)}
        
        TEXT TO ANALYZE:
        {prompt}
        
        INSTRUCTIONS:
        - Extract all relevant information that matches the schema
        - Provide output as a valid JSON object following the schema
        - If information is not present, use null values
        - DO NOT include any explanations or text outside the JSON
        - Ensure the output is a valid JSON object
        
        JSON OUTPUT:
        """
        
        model_name = model or self.default_model
        
        # Generate the structured data response
        response_text = self.generate_text(
            extraction_prompt,
            model=model_name,
            temperature=temperature
        )
        
        # Parse the JSON response
        try:
            # Clean up the response text to ensure it's valid JSON
            cleaned_response = response_text.strip()
            # Remove markdown code block syntax if present
            if cleaned_response.startswith("```json"):
                cleaned_response = cleaned_response[7:]
            if cleaned_response.endswith("```"):
                cleaned_response = cleaned_response[:-3]
            
            structured_data = json.loads(cleaned_response.strip())
            return structured_data
        except json.JSONDecodeError as e:
            raise ValueError(f"Failed to parse JSON response: {e}\nResponse: {response_text}")
    
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
