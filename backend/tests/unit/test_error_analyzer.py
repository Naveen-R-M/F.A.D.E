"""
Unit tests for the ErrorAnalyzer class.
"""

import unittest
from unittest.mock import MagicMock, patch

from agents.target_selector.error_analyzer import ErrorAnalyzer


class TestErrorAnalyzer(unittest.TestCase):
    """Tests for the ErrorAnalyzer class."""
    
    def setUp(self):
        # Create a mock LLM client
        self.mock_llm_client = MagicMock()
        self.mock_llm_client.generate_text.return_value = '''
        {
            "error_type": "network_error",
            "probable_cause": "Connection timeout",
            "suggested_actions": [
                {
                    "action": "retry",
                    "parameters": {"backoff": true, "wait_time": 30},
                    "reason": "Network issues are often transient"
                }
            ],
            "confidence": 0.8
        }
        '''
        
        # Create the error analyzer
        self.error_analyzer = ErrorAnalyzer(self.mock_llm_client)
    
    def test_determine_error_source(self):
        """Test error source determination."""
        # Test UniProt error source
        error_message = "Error accessing UniProt API: 404 Not Found"
        context = {"operation": "get_protein"}
        source = self.error_analyzer._determine_error_source(error_message, context)
        self.assertEqual(source, "uniprot")
        
        # Test parsing error source
        error_message = "Error parsing JSON response"
        context = {"operation": "process_response"}
        source = self.error_analyzer._determine_error_source(error_message, context)
        self.assertEqual(source, "parsing")
        
        # Test network error source
        error_message = "Connection timed out"
        context = {"operation": "api_call"}
        source = self.error_analyzer._determine_error_source(error_message, context)
        self.assertEqual(source, "network")
        
        # Test validation error source
        error_message = "Sequence validation failed"
        context = {"operation": "validate"}
        source = self.error_analyzer._determine_error_source(error_message, context)
        self.assertEqual(source, "validation")
        
        # Test generic error source
        error_message = "Unknown error"
        context = {"operation": "unknown"}
        source = self.error_analyzer._determine_error_source(error_message, context)
        self.assertEqual(source, "generic")
    
    def test_extract_error_details(self):
        """Test error details extraction."""
        # Test HTTP status code extraction
        error_message = "API request failed with status code 404"
        source = "http"
        details = self.error_analyzer._extract_error_details(error_message, source)
        self.assertEqual(details["status_code"], 404)
        
        # Test URL extraction
        error_message = "Failed to access https://www.uniprot.org/api/proteins"
        source = "network"
        details = self.error_analyzer._extract_error_details(error_message, source)
        self.assertEqual(details["url"], "https://www.uniprot.org/api/proteins")
        
        # Test identifier extraction
        error_message = "Protein identifier 'P12345' not found"
        source = "uniprot"
        details = self.error_analyzer._extract_error_details(error_message, source)
        self.assertEqual(details["identifier"], "P12345")
    
    def test_analyze_with_llm(self):
        """Test error analysis using LLM."""
        # Test successful LLM analysis
        error_message = "Connection timed out accessing UniProt API"
        error_details = {"source": "network", "full_message": error_message}
        error_source = "network"
        context = {"operation": "fetch_protein"}
        
        analysis = self.error_analyzer._analyze_with_llm(
            error_message, error_details, error_source, context
        )
        
        # Check that the LLM client was called
        self.mock_llm_client.generate_text.assert_called()
        
        # Check that the analysis has the expected structure
        self.assertEqual(analysis["error_type"], "network_error")
        self.assertEqual(analysis["probable_cause"], "Connection timeout")
        self.assertEqual(len(analysis["suggested_actions"]), 1)
        self.assertEqual(analysis["suggested_actions"][0]["action"], "retry")
        self.assertEqual(analysis["confidence"], 0.8)
    
    def test_analyze_with_patterns(self):
        """Test error analysis using pattern matching."""
        # Test UniProt "not found" error
        error_message = "Protein not found: KRAS"
        error_details = {"source": "uniprot", "full_message": error_message}
        error_source = "uniprot"
        context = {"operation": "fetch_protein"}
        
        analysis = self.error_analyzer._analyze_with_patterns(
            error_message, error_details, error_source, context
        )
        
        # Check that the analysis has the expected structure
        self.assertEqual(analysis["error_type"], "resource_not_found")
        self.assertTrue(any(action["action"] == "modify_query" for action in analysis["suggested_actions"]))
        
        # Test rate limit error
        error_message = "Rate limit exceeded: 429 Too Many Requests"
        error_details = {"source": "uniprot", "full_message": error_message, "status_code": 429}
        
        analysis = self.error_analyzer._analyze_with_patterns(
            error_message, error_details, error_source, context
        )
        
        # Check that the analysis has the expected structure
        self.assertEqual(analysis["error_type"], "rate_limit")
        self.assertTrue(any(action["action"] == "retry" for action in analysis["suggested_actions"]))
        
        # Test server error
        error_message = "Server error: 500 Internal Server Error"
        error_details = {"source": "http", "full_message": error_message, "status_code": 500}
        error_source = "http"
        
        analysis = self.error_analyzer._analyze_with_patterns(
            error_message, error_details, error_source, context
        )
        
        # Check that the analysis has the expected structure
        self.assertEqual(analysis["error_type"], "server_error")
        self.assertTrue(any(action["action"] == "retry" for action in analysis["suggested_actions"]))
    
    def test_analyze_full(self):
        """Test the full analyze method."""
        # Test with a typical error
        error_message = "Connection timed out accessing UniProt API"
        context = {"operation": "fetch_protein"}
        
        analysis = self.error_analyzer.analyze(error_message, context)
        
        # Check that the analysis has the expected structure
        self.assertIn("error_type", analysis)
        self.assertIn("probable_cause", analysis)
        self.assertIn("suggested_actions", analysis)
        self.assertTrue(len(analysis["suggested_actions"]) > 0)
        
        # Test with LLM failure (should fall back to pattern matching)
        self.mock_llm_client.generate_text.side_effect = Exception("LLM API error")
        
        analysis = self.error_analyzer.analyze(error_message, context)
        
        # Check that the analysis still has the expected structure
        self.assertIn("error_type", analysis)
        self.assertIn("probable_cause", analysis)
        self.assertIn("suggested_actions", analysis)
        self.assertTrue(len(analysis["suggested_actions"]) > 0)


if __name__ == "__main__":
    unittest.main()
