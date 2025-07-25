"""
Error analyzer component for the Target Selector agent.

This component analyzes API errors using LLM to determine the cause and suggest solutions.
"""

from typing import Any, Dict, List, Optional
import json
import re
from utils.gemini_client import GeminiClient

class ErrorAnalyzer:
    """
    Analyzes API errors using LLM to determine cause and solution.
    """
    
    def __init__(self, llm_client: Optional[GeminiClient] = None) -> None:
        """
        Initialize the error analyzer.
        
        Args:
            llm_client: Client for the LLM (e.g., GeminiClient)
        """
        self.llm_client = llm_client
    
    def analyze(self, error_message: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Analyze error message and return diagnosis and suggested fixes.
        
        Args:
            error_message: The error message from the API
            context: Additional context about the operation that caused the error
            
        Returns:
            Dictionary containing diagnosis and suggested fixes
        """
        if context is None:
            context = {}
        
        # Determine the error source based on message and context
        error_source = self._determine_error_source(error_message, context)
        
        # Extract relevant parts of the error message
        error_details = self._extract_error_details(error_message, error_source)
        
        # If LLM client is available, use it for advanced analysis
        if self.llm_client:
            analysis = self._analyze_with_llm(error_message, error_details, error_source, context)
        else:
            # Default analysis based on pattern matching
            analysis = self._analyze_with_patterns(error_message, error_details, error_source, context)
        
        return analysis
    
    def _determine_error_source(self, error_message: str, context: Dict[str, Any]) -> str:
        """
        Determine the source of the error based on the message and context.
        
        Args:
            error_message: The error message
            context: Additional context
            
        Returns:
            Source of the error (e.g., 'uniprot', 'network', 'parsing')
        """
        # Check for common sources in the error message
        if 'UniProt' in error_message or 'uniprot.org' in error_message:
            return 'uniprot'
        elif 'JSON' in error_message or 'json' in error_message or 'parsing' in error_message:
            return 'parsing'
        elif 'timed out' in error_message or 'connection' in error_message or 'network' in error_message:
            return 'network'
        elif 'response code' in error_message or 'status code' in error_message:
            return 'http'
        elif 'sequence' in error_message or 'validation' in error_message:
            return 'validation'
        elif 'API key' in error_message or 'authentication' in error_message:
            return 'authentication'
        
        # Check context for operation type
        operation = context.get('operation', '').lower()
        if 'uniprot' in operation:
            return 'uniprot'
        elif 'sequence' in operation:
            return 'sequence'
        elif 'validation' in operation:
            return 'validation'
        
        # Default to generic
        return 'generic'
    
    def _extract_error_details(self, error_message: str, error_source: str) -> Dict[str, Any]:
        """
        Extract relevant details from the error message.
        
        Args:
            error_message: The error message
            error_source: Source of the error
            
        Returns:
            Dictionary with extracted details
        """
        details = {
            'source': error_source,
            'full_message': error_message
        }
        
        # Extract HTTP status code if present
        status_match = re.search(r'status(?:\s+code)?(?:\:|\s+)?\s*(\d+)', error_message, re.IGNORECASE)
        if status_match:
            details['status_code'] = int(status_match.group(1))
        
        # Extract request URL if present
        url_match = re.search(r'(https?://[^\s\'"]+)', error_message)
        if url_match:
            details['url'] = url_match.group(1)
        
        # Extract protein/gene identifiers if present
        id_match = re.search(r'(?:protein|gene|accession|identifier)(?:\s+id)?(?:\:|\s+)?\s*[\'"]?([A-Za-z0-9_\-\.]+)[\'"]?', 
                            error_message, re.IGNORECASE)
        if id_match:
            details['identifier'] = id_match.group(1)
        
        return details
    
    def _analyze_with_llm(
        self, 
        error_message: str, 
        error_details: Dict[str, Any], 
        error_source: str, 
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Analyze error using LLM.
        
        Args:
            error_message: The error message
            error_details: Extracted error details
            error_source: Source of the error
            context: Additional context
            
        Returns:
            Analysis dictionary
        """
        # Format context information
        context_info = json.dumps(context, indent=2)
        error_details_info = json.dumps(error_details, indent=2)
        
        # Construct prompt based on error source
        if error_source == 'uniprot':
            prompt = f"""
            Analyze the following UniProt API error:
            
            Error Message: {error_message}
            
            Error Details: {error_details_info}
            
            Context: {context_info}
            
            Please analyze this error and suggest solutions. Common issues with UniProt API include:
            - Invalid protein/gene identifiers
            - Malformed query parameters
            - Rate limiting
            - API endpoint changes
            - Network issues
            
            Provide your response as a JSON object with the following structure:
            {{
                "error_type": "The specific type of error",
                "probable_cause": "Your analysis of the most likely cause",
                "suggested_actions": [
                    {{
                        "action": "action_name", // e.g., retry, modify_query, use_alternative_endpoint
                        "parameters": {{ 
                            // Parameters specific to the action
                        }},
                        "reason": "Why this action might help"
                    }}
                ],
                "confidence": 0.7 // Your confidence in this analysis (0.0 to 1.0)
            }}
            
            IMPORTANT: Ensure the JSON is valid and complete. Focus on actionable solutions.
            """
        elif error_source == 'parsing':
            prompt = f"""
            Analyze the following data parsing error:
            
            Error Message: {error_message}
            
            Error Details: {error_details_info}
            
            Context: {context_info}
            
            Please analyze this error and suggest solutions. Common issues with parsing include:
            - Malformed JSON/XML responses
            - Unexpected data structure
            - Missing fields
            - Type mismatches
            
            Provide your response as a JSON object with the following structure:
            {{
                "error_type": "The specific type of error",
                "probable_cause": "Your analysis of the most likely cause",
                "suggested_actions": [
                    {{
                        "action": "action_name", // e.g., retry, modify_parsing, handle_missing_fields
                        "parameters": {{ 
                            // Parameters specific to the action
                        }},
                        "reason": "Why this action might help"
                    }}
                ],
                "confidence": 0.7 // Your confidence in this analysis (0.0 to 1.0)
            }}
            
            IMPORTANT: Ensure the JSON is valid and complete. Focus on actionable solutions.
            """
        else:
            # Generic prompt for other error types
            prompt = f"""
            Analyze the following API error:
            
            Error Message: {error_message}
            
            Error Details: {error_details_info}
            
            Context: {context_info}
            
            Please analyze this error and suggest solutions. Consider common issues like:
            - Invalid parameters or inputs
            - Authentication problems
            - Rate limiting
            - Network issues
            - API changes
            
            Provide your response as a JSON object with the following structure:
            {{
                "error_type": "The specific type of error",
                "probable_cause": "Your analysis of the most likely cause",
                "suggested_actions": [
                    {{
                        "action": "action_name", // e.g., retry, modify_parameters, use_alternative
                        "parameters": {{ 
                            // Parameters specific to the action
                        }},
                        "reason": "Why this action might help"
                    }}
                ],
                "confidence": 0.7 // Your confidence in this analysis (0.0 to 1.0)
            }}
            
            IMPORTANT: Ensure the JSON is valid and complete. Focus on actionable solutions.
            """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.2)
            
            # Extract JSON from response
            json_start = response.find('{')
            json_end = response.rfind('}') + 1
            
            if json_start >= 0 and json_end > json_start:
                json_str = response[json_start:json_end]
                analysis = json.loads(json_str)
                
                # Ensure required fields are present
                if 'error_type' not in analysis:
                    analysis['error_type'] = error_source
                if 'probable_cause' not in analysis:
                    analysis['probable_cause'] = 'Unknown'
                if 'suggested_actions' not in analysis:
                    analysis['suggested_actions'] = [{'action': 'retry', 'parameters': {}, 'reason': 'Simple retry'}]
                if 'confidence' not in analysis:
                    analysis['confidence'] = 0.5
                
                return analysis
            else:
                # Fallback to pattern matching if JSON extraction fails
                return self._analyze_with_patterns(error_message, error_details, error_source, context)
        
        except Exception as e:
            # Fallback to pattern matching if LLM fails
            return self._analyze_with_patterns(error_message, error_details, error_source, context)
    
    def _analyze_with_patterns(
        self, 
        error_message: str, 
        error_details: Dict[str, Any], 
        error_source: str, 
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Analyze error using pattern matching.
        
        Args:
            error_message: The error message
            error_details: Extracted error details
            error_source: Source of the error
            context: Additional context
            
        Returns:
            Analysis dictionary
        """
        # Default analysis structure
        analysis = {
            'error_type': error_source,
            'probable_cause': 'Unknown',
            'suggested_actions': [
                {
                    'action': 'retry',
                    'parameters': {},
                    'reason': 'Simple retry might resolve transient issues'
                }
            ],
            'confidence': 0.3
        }
        
        # Pattern matching based on error source
        if error_source == 'uniprot':
            # Check for common UniProt errors
            if 'not found' in error_message.lower():
                analysis['error_type'] = 'resource_not_found'
                analysis['probable_cause'] = 'The requested protein or gene was not found in UniProt'
                analysis['suggested_actions'] = [
                    {
                        'action': 'modify_query',
                        'parameters': {'use_search': True},
                        'reason': 'Direct lookup failed, try searching instead'
                    },
                    {
                        'action': 'use_alternative_source',
                        'parameters': {'source': 'ncbi'},
                        'reason': 'Try an alternative database'
                    }
                ]
                analysis['confidence'] = 0.7
            elif 'rate limit' in error_message.lower() or '429' in error_message:
                analysis['error_type'] = 'rate_limit'
                analysis['probable_cause'] = 'Rate limit exceeded for UniProt API'
                analysis['suggested_actions'] = [
                    {
                        'action': 'retry',
                        'parameters': {'backoff': True, 'wait_time': 30},
                        'reason': 'Wait longer before retrying'
                    }
                ]
                analysis['confidence'] = 0.8
            elif 'timeout' in error_message.lower():
                analysis['error_type'] = 'timeout'
                analysis['probable_cause'] = 'Request timed out'
                analysis['suggested_actions'] = [
                    {
                        'action': 'retry',
                        'parameters': {'timeout': 60},
                        'reason': 'Retry with longer timeout'
                    }
                ]
                analysis['confidence'] = 0.6
        elif error_source == 'network':
            analysis['error_type'] = 'network_error'
            analysis['probable_cause'] = 'Network connectivity issue'
            analysis['suggested_actions'] = [
                {
                    'action': 'retry',
                    'parameters': {'backoff': True, 'max_retries': 5},
                    'reason': 'Network issues are often transient'
                }
            ]
            analysis['confidence'] = 0.7
        elif error_source == 'parsing':
            analysis['error_type'] = 'parsing_error'
            analysis['probable_cause'] = 'Failed to parse API response'
            analysis['suggested_actions'] = [
                {
                    'action': 'modify_parsing',
                    'parameters': {'format': 'text'},
                    'reason': 'Try parsing as plain text first'
                },
                {
                    'action': 'retry',
                    'parameters': {},
                    'reason': 'API response might be temporarily malformed'
                }
            ]
            analysis['confidence'] = 0.5
        elif error_source == 'validation':
            analysis['error_type'] = 'validation_error'
            analysis['probable_cause'] = 'Sequence validation failed'
            analysis['suggested_actions'] = [
                {
                    'action': 'skip_validation',
                    'parameters': {'warning': True},
                    'reason': 'Continue with a warning about validation failure'
                },
                {
                    'action': 'use_alternative_source',
                    'parameters': {'source': 'ncbi'},
                    'reason': 'Try retrieving the sequence from an alternative source'
                }
            ]
            analysis['confidence'] = 0.6
        
        # Check for HTTP status codes in details
        if 'status_code' in error_details:
            status_code = error_details['status_code']
            if status_code == 404:
                analysis['error_type'] = 'resource_not_found'
                analysis['probable_cause'] = 'The requested resource was not found'
                analysis['suggested_actions'] = [
                    {
                        'action': 'modify_query',
                        'parameters': {'use_search': True},
                        'reason': 'Direct lookup failed, try searching instead'
                    }
                ]
                analysis['confidence'] = 0.8
            elif status_code == 429:
                analysis['error_type'] = 'rate_limit'
                analysis['probable_cause'] = 'Rate limit exceeded'
                analysis['suggested_actions'] = [
                    {
                        'action': 'retry',
                        'parameters': {'backoff': True, 'wait_time': 30},
                        'reason': 'Wait longer before retrying'
                    }
                ]
                analysis['confidence'] = 0.8
            elif status_code >= 500:
                analysis['error_type'] = 'server_error'
                analysis['probable_cause'] = 'Server-side error'
                analysis['suggested_actions'] = [
                    {
                        'action': 'retry',
                        'parameters': {'backoff': True, 'max_retries': 3},
                        'reason': 'Server errors are often temporary'
                    }
                ]
                analysis['confidence'] = 0.7
        
        return analysis
