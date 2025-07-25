"""
Agentic mixin for F.A.D.E agents

This mixin provides agentic capabilities to base agents, enabling them
to make autonomous decisions, recover from errors, and learn from interactions.
"""

from typing import Any, Dict, List, Optional, Tuple, Union
import logging
import json
import os
import time
from collections import defaultdict

class AgenticMixin:
    """
    Mixin class that adds agentic capabilities to any agent.
    
    This mixin provides common functionality for:
    - Error recovery
    - LLM-based decision making
    - Learning from experience
    """
    
    def initialize_agentic_components(
        self, 
        llm_client=None, 
        memory_file: Optional[str] = None,
        max_attempts: int = 3,
        backoff_factor: float = 1.5
    ) -> None:
        """
        Initialize components needed for agentic behavior.
        
        Args:
            llm_client: Client for interacting with LLM (e.g., GeminiClient)
            memory_file: Path to file for persisting learned interactions
            max_attempts: Maximum number of retry attempts
            backoff_factor: Factor to increase wait time between retries
        """
        self.llm_client = llm_client
        self.memory_file = memory_file
        self.max_attempts = max_attempts
        self.backoff_factor = backoff_factor
        
        # Initialize interaction memory
        self.interaction_memory = self._load_memory() if memory_file else defaultdict(list)
        
        # Make sure logger is available
        if not hasattr(self, 'logger'):
            self.logger = logging.getLogger(f"fade.agent.{self.__class__.__name__}")
        
        self.logger.info("Initialized agentic components")
    
    def handle_error(
        self, 
        error: Exception, 
        context: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Handle errors using LLM-based analysis and recovery.
        
        Args:
            error: The error that occurred
            context: Additional context about the error
            
        Returns:
            Dictionary with error analysis and recovery suggestions
        """
        if context is None:
            context = {}
            
        error_message = str(error)
        error_type = type(error).__name__
        
        self.logger.info(f"Handling error: {error_type} - {error_message}")
        
        # Check if we've seen this error before
        if 'operation' in context:
            operation_key = context['operation']
            for memory in self.interaction_memory.get(operation_key, []):
                if memory['error_type'] == error_type and memory['error_message'] == error_message:
                    self.logger.info(f"Found similar error in memory: {error_type}")
                    if memory['resolution_successful']:
                        return memory['resolution_strategy']
        
        # If we haven't seen this error before or previous resolutions failed,
        # use LLM to analyze the error
        if self.llm_client:
            analysis = self._analyze_error_with_llm(error_message, error_type, context)
        else:
            # Default analysis if no LLM is available
            analysis = {
                'error_type': error_type,
                'error_message': error_message,
                'probable_cause': 'Unknown',
                'suggested_actions': [
                    {'action': 'retry', 'parameters': {}}
                ],
                'confidence': 0.3
            }
        
        self.logger.info(f"Error analysis: {analysis}")
        return analysis
    
    def make_decision(
        self, 
        options: List[Dict[str, Any]], 
        context: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Use LLM to make decisions between multiple options.
        
        Args:
            options: List of option dictionaries
            context: Context information for the decision
            
        Returns:
            Selected option with additional reasoning
        """
        if not options:
            raise ValueError("No options provided for decision making")
            
        if context is None:
            context = {}
            
        self.logger.info(f"Making decision between {len(options)} options")
        
        # Check if we've made a similar decision before
        if 'decision_type' in context:
            decision_key = context['decision_type']
            context_key = json.dumps(context.get('key_factors', {}))
            
            for memory in self.interaction_memory.get(f"decision:{decision_key}", []):
                if memory['context_key'] == context_key and memory['success']:
                    matching_option = next((opt for opt in options if self._option_matches(opt, memory['selected_option'])), None)
                    if matching_option:
                        self.logger.info(f"Found successful similar decision in memory")
                        return {
                            'selected_option': matching_option,
                            'reasoning': memory['reasoning'],
                            'confidence': memory['confidence'],
                            'from_memory': True
                        }
        
        # If no memory match or LLM not available, use simple heuristic
        if not self.llm_client:
            # Default to first option if no LLM is available
            selected = options[0]
            return {
                'selected_option': selected,
                'reasoning': 'Selected based on order (no LLM available)',
                'confidence': 0.5,
                'from_memory': False
            }
        
        # Use LLM to make the decision
        decision = self._make_decision_with_llm(options, context)
        self.logger.info(f"Decision made with confidence {decision['confidence']}")
        return decision
    
    def learn_from_interaction(
        self, 
        interaction_data: Dict[str, Any]
    ) -> None:
        """
        Update internal knowledge based on interaction outcomes.
        
        Args:
            interaction_data: Data about the interaction and its outcome
        """
        if not interaction_data:
            return
            
        self.logger.info(f"Learning from interaction: {interaction_data.get('interaction_type', 'unknown')}")
        
        # Extract key elements
        interaction_type = interaction_data.get('interaction_type', 'unknown')
        success = interaction_data.get('success', False)
        
        # Store in memory based on interaction type
        if interaction_type in ['error_recovery', 'decision', 'search_strategy', 'validation']:
            key = interaction_data.get('key', interaction_type)
            self.interaction_memory[key].append({
                'timestamp': time.time(),
                'success': success,
                **interaction_data
            })
            
            # Limit memory size
            if len(self.interaction_memory[key]) > 100:
                self.interaction_memory[key] = self.interaction_memory[key][-100:]
        
        # Persist to file if configured
        if self.memory_file:
            self._save_memory()
    
    def _analyze_error_with_llm(
        self, 
        error_message: str, 
        error_type: str, 
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Use LLM to analyze an error and suggest recovery actions.
        
        Args:
            error_message: The error message
            error_type: The type of error
            context: Additional context about the error
            
        Returns:
            Dictionary with error analysis and recovery suggestions
        """
        # Construct a prompt for the LLM
        operation = context.get('operation', 'unknown operation')
        additional_context = json.dumps(context.get('additional_context', {}))
        
        prompt = f"""
        Analyze the following error that occurred during {operation}:
        
        Error Type: {error_type}
        Error Message: {error_message}
        
        Additional Context:
        {additional_context}
        
        Please analyze the error and suggest potential recovery actions. Provide your response as a JSON object with the following structure:
        {{
            "error_type": "The type of error",
            "error_message": "The error message",
            "probable_cause": "Your analysis of the most likely cause",
            "suggested_actions": [
                {{
                    "action": "Name of action (e.g., retry, modify_query, use_alternative_endpoint)",
                    "parameters": {{ 
                        // Parameters specific to the action
                    }}
                }}
            ],
            "confidence": 0.7 // Your confidence in this analysis (0.0 to 1.0)
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
                    self.logger.warning(f"Could not find JSON in LLM response")
                    # Return a basic analysis
                    return {
                        'error_type': error_type,
                        'error_message': error_message,
                        'probable_cause': 'Unknown (could not parse LLM response)',
                        'suggested_actions': [
                            {'action': 'retry', 'parameters': {}}
                        ],
                        'confidence': 0.3
                    }
            except json.JSONDecodeError:
                self.logger.warning(f"Could not parse JSON from LLM response")
                # Return a basic analysis
                return {
                    'error_type': error_type,
                    'error_message': error_message,
                    'probable_cause': 'Unknown (could not parse LLM response)',
                    'suggested_actions': [
                        {'action': 'retry', 'parameters': {}}
                    ],
                    'confidence': 0.3
                }
        except Exception as e:
            self.logger.error(f"Error calling LLM: {e}")
            # Return a basic analysis
            return {
                'error_type': error_type,
                'error_message': error_message,
                'probable_cause': 'Unknown (error calling LLM)',
                'suggested_actions': [
                    {'action': 'retry', 'parameters': {}}
                ],
                'confidence': 0.2
            }
    
    def _make_decision_with_llm(
        self, 
        options: List[Dict[str, Any]], 
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Use LLM to make a decision between multiple options.
        
        Args:
            options: List of option dictionaries
            context: Context information for the decision
            
        Returns:
            Selected option with additional reasoning
        """
        # Construct a prompt for the LLM
        decision_type = context.get('decision_type', 'unknown decision')
        criteria = json.dumps(context.get('criteria', {}))
        options_json = json.dumps(options)
        
        prompt = f"""
        You need to make a decision about {decision_type} based on the following criteria:
        
        Criteria:
        {criteria}
        
        Available options:
        {options_json}
        
        Please select the best option and explain your reasoning. Provide your response as a JSON object with the following structure:
        {{
            "selected_option_index": 0, // Index of the selected option in the provided list
            "reasoning": "Detailed explanation of why this option was selected",
            "confidence": 0.8 // Your confidence in this decision (0.0 to 1.0)
        }}
        """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.3)
            
            # Extract JSON from response
            try:
                # Find JSON in the response
                json_start = response.find('{')
                json_end = response.rfind('}') + 1
                if json_start >= 0 and json_end > json_start:
                    json_str = response[json_start:json_end]
                    decision_data = json.loads(json_str)
                    
                    # Get the selected option
                    selected_index = decision_data.get('selected_option_index', 0)
                    if selected_index < 0 or selected_index >= len(options):
                        selected_index = 0
                    
                    return {
                        'selected_option': options[selected_index],
                        'reasoning': decision_data.get('reasoning', 'No reasoning provided'),
                        'confidence': decision_data.get('confidence', 0.5),
                        'from_memory': False
                    }
                else:
                    self.logger.warning(f"Could not find JSON in LLM response")
                    # Default to first option
                    return {
                        'selected_option': options[0],
                        'reasoning': 'Default selection (could not parse LLM response)',
                        'confidence': 0.3,
                        'from_memory': False
                    }
            except json.JSONDecodeError:
                self.logger.warning(f"Could not parse JSON from LLM response")
                # Default to first option
                return {
                    'selected_option': options[0],
                    'reasoning': 'Default selection (could not parse LLM response)',
                    'confidence': 0.3,
                    'from_memory': False
                }
        except Exception as e:
            self.logger.error(f"Error calling LLM: {e}")
            # Default to first option
            return {
                'selected_option': options[0],
                'reasoning': f'Default selection (error calling LLM: {e})',
                'confidence': 0.2,
                'from_memory': False
            }
    
    def _option_matches(self, option1: Dict[str, Any], option2: Dict[str, Any]) -> bool:
        """
        Check if two options are similar enough to be considered the same.
        
        Args:
            option1: First option
            option2: Second option
            
        Returns:
            True if options match, False otherwise
        """
        # Simple matching based on key fields
        # Can be extended with more sophisticated matching
        if 'id' in option1 and 'id' in option2:
            return option1['id'] == option2['id']
        
        # Check key fields for similarity
        similarity_count = 0
        total_fields = 0
        
        for key in option1:
            if key in option2:
                total_fields += 1
                if option1[key] == option2[key]:
                    similarity_count += 1
        
        # Consider match if 80% of fields match
        return total_fields > 0 and (similarity_count / total_fields) >= 0.8
    
    def _load_memory(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Load interaction memory from file.
        
        Returns:
            Dictionary of interaction memories
        """
        if not self.memory_file or not os.path.exists(self.memory_file):
            return defaultdict(list)
            
        try:
            with open(self.memory_file, 'r') as f:
                data = json.load(f)
                # Convert to defaultdict
                memory = defaultdict(list)
                for key, value in data.items():
                    memory[key] = value
                return memory
        except Exception as e:
            self.logger.error(f"Error loading memory file: {e}")
            return defaultdict(list)
    
    def _save_memory(self) -> None:
        """Save interaction memory to file."""
        if not self.memory_file:
            return
            
        try:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(self.memory_file), exist_ok=True)
            
            # Convert defaultdict to regular dict for serialization
            memory_dict = dict(self.interaction_memory)
            
            with open(self.memory_file, 'w') as f:
                json.dump(memory_dict, f, indent=2)
        except Exception as e:
            self.logger.error(f"Error saving memory file: {e}")
    
    def execute_with_retry(
        self, 
        func, 
        *args, 
        operation_name: str = "operation", 
        max_attempts: Optional[int] = None, 
        **kwargs
    ) -> Any:
        """
        Execute a function with automatic retry and error handling.
        
        Args:
            func: Function to execute
            *args: Arguments to pass to the function
            operation_name: Name of the operation for logging
            max_attempts: Maximum number of attempts (overrides default)
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Result of the function execution
        """
        attempts = 0
        max_retry = max_attempts if max_attempts is not None else self.max_attempts
        wait_time = 1.0  # Initial wait time in seconds
        
        while attempts < max_retry:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                attempts += 1
                self.logger.warning(f"{operation_name} failed (attempt {attempts}/{max_retry}): {e}")
                
                if attempts >= max_retry:
                    self.logger.error(f"{operation_name} failed after {max_retry} attempts")
                    raise
                
                # Analyze error and decide how to retry
                error_analysis = self.handle_error(e, {'operation': operation_name})
                
                # Extract action from analysis
                suggested_actions = error_analysis.get('suggested_actions', [])
                if not suggested_actions:
                    # Default to simple retry
                    action = {'action': 'retry', 'parameters': {}}
                else:
                    action = suggested_actions[0]
                
                action_type = action.get('action', 'retry')
                parameters = action.get('parameters', {})
                
                # Apply action
                if action_type == 'retry':
                    # Simple retry with backoff
                    wait_time *= self.backoff_factor
                    self.logger.info(f"Retrying {operation_name} after {wait_time:.2f}s")
                    time.sleep(wait_time)
                elif action_type == 'modify_parameters':
                    # Update kwargs with new parameters
                    for param_key, param_value in parameters.items():
                        kwargs[param_key] = param_value
                    self.logger.info(f"Retrying {operation_name} with modified parameters")
                elif action_type == 'abort':
                    # Abort retry
                    self.logger.info(f"Aborting {operation_name} based on error analysis")
                    raise
                else:
                    # Unknown action, just retry
                    wait_time *= self.backoff_factor
                    self.logger.info(f"Retrying {operation_name} with unknown action: {action_type}")
                    time.sleep(wait_time)
                
                # Record the interaction for learning
                self.learn_from_interaction({
                    'interaction_type': 'error_recovery',
                    'key': operation_name,
                    'error_type': type(e).__name__,
                    'error_message': str(e),
                    'action_taken': action,
                    'success': False  # We'll only know success after retry
                })
