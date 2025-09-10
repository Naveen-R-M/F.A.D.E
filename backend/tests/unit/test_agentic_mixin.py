"""
Unit tests for the AgenticMixin class.
"""

import os
import unittest
from unittest.mock import MagicMock, patch
import json
import tempfile

from agents.base.agentic_mixin import AgenticMixin


class TestAgent(AgenticMixin):
    """Test agent that uses the AgenticMixin."""
    
    def __init__(self, llm_client=None, memory_file=None):
        self.logger = MagicMock()
        self.initialize_agentic_components(llm_client, memory_file)


class TestAgenticMixin(unittest.TestCase):
    """Tests for the AgenticMixin class."""
    
    def setUp(self):
        # Create a mock LLM client
        self.mock_llm_client = MagicMock()
        self.mock_llm_client.generate_text.return_value = '{"selected_option_index": 0, "reasoning": "This is the best option", "confidence": 0.8}'
        
        # Create a temporary file for memory
        self.temp_file = tempfile.NamedTemporaryFile(delete=False)
        self.temp_file.close()
        
        # Create the test agent
        self.agent = TestAgent(self.mock_llm_client, self.temp_file.name)
    
    def tearDown(self):
        # Clean up the temporary file
        os.unlink(self.temp_file.name)
    
    def test_initialization(self):
        """Test that the agent is initialized correctly."""
        self.assertEqual(self.agent.llm_client, self.mock_llm_client)
        self.assertEqual(self.agent.memory_file, self.temp_file.name)
        self.assertEqual(self.agent.max_attempts, 3)  # Default value
        
    def test_handle_error(self):
        """Test the error handling functionality."""
        # Create a test error
        test_error = ValueError("Test error")
        
        # Handle the error
        result = self.agent.handle_error(test_error, {"operation": "test_operation"})
        
        # Check that the LLM client was called
        self.mock_llm_client.generate_text.assert_called()
        
        # Check that the result has the expected structure
        self.assertIn("error_type", result)
        self.assertIn("error_message", result)
        self.assertIn("suggested_actions", result)
        
    def test_make_decision(self):
        """Test the decision making functionality."""
        # Create test options
        options = [
            {"id": 1, "name": "Option 1"},
            {"id": 2, "name": "Option 2"}
        ]
        
        # Make a decision
        result = self.agent.make_decision(options, {"decision_type": "test_decision"})
        
        # Check that the LLM client was called
        self.mock_llm_client.generate_text.assert_called()
        
        # Check that the result has the expected structure
        self.assertIn("selected_option", result)
        self.assertIn("reasoning", result)
        self.assertIn("confidence", result)
        
        # Check that the selected option is one of the provided options
        self.assertIn(result["selected_option"], options)
        
    def test_learn_from_interaction(self):
        """Test the learning functionality."""
        # Create test interaction data
        interaction_data = {
            "interaction_type": "test_interaction",
            "key": "test_key",
            "success": True,
            "details": {"test": "data"}
        }
        
        # Learn from the interaction
        self.agent.learn_from_interaction(interaction_data)
        
        # Check that the interaction was stored in memory
        self.assertIn("test_key", self.agent.interaction_memory)
        self.assertEqual(len(self.agent.interaction_memory["test_key"]), 1)
        
        # Check that the memory was saved to file
        with open(self.temp_file.name, "r") as f:
            memory_data = json.load(f)
            self.assertIn("test_key", memory_data)
            self.assertEqual(len(memory_data["test_key"]), 1)
            
    def test_execute_with_retry_success(self):
        """Test the execute_with_retry function with a successful execution."""
        # Create a mock function that succeeds
        mock_func = MagicMock(return_value="success")
        
        # Execute the function with retry
        result = self.agent.execute_with_retry(mock_func, "arg1", "arg2", operation_name="test_operation")
        
        # Check that the function was called once
        mock_func.assert_called_once_with("arg1", "arg2")
        
        # Check that the result is correct
        self.assertEqual(result, "success")
        
    def test_execute_with_retry_failure(self):
        """Test the execute_with_retry function with failures that eventually succeed."""
        # Create a mock function that fails twice then succeeds
        mock_func = MagicMock(side_effect=[ValueError("Error 1"), ValueError("Error 2"), "success"])
        
        # Execute the function with retry
        result = self.agent.execute_with_retry(mock_func, "arg1", operation_name="test_operation")
        
        # Check that the function was called three times
        self.assertEqual(mock_func.call_count, 3)
        
        # Check that the result is correct
        self.assertEqual(result, "success")
        
    def test_execute_with_retry_max_failures(self):
        """Test the execute_with_retry function with too many failures."""
        # Create a mock function that always fails
        mock_func = MagicMock(side_effect=ValueError("Error"))
        
        # Execute the function with retry (should raise an exception)
        with self.assertRaises(ValueError):
            self.agent.execute_with_retry(mock_func, "arg1", max_attempts=2, operation_name="test_operation")
        
        # Check that the function was called twice
        self.assertEqual(mock_func.call_count, 2)


if __name__ == "__main__":
    unittest.main()
