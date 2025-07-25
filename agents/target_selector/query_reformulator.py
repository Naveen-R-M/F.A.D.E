"""
Query reformulator component for the Target Selector agent.

This component reformulates failed queries into potentially valid formats.
"""

from typing import Any, Dict, List, Optional, Tuple, Union
import json
import re
from utils.gemini_client import GeminiClient

class QueryReformulator:
    """
    Reformulates failed queries into potentially valid formats.
    """
    
    def __init__(self, llm_client: Optional[GeminiClient] = None) -> None:
        """
        Initialize the query reformulator.
        
        Args:
            llm_client: Client for the LLM (e.g., GeminiClient)
        """
        self.llm_client = llm_client
        self.query_patterns = self._initialize_patterns()
    
    def reformulate(
        self, 
        original_query: Dict[str, Any], 
        error_analysis: Dict[str, Any],
        target_info: Optional[Dict[str, Any]] = None
    ) -> List[Dict[str, Any]]:
        """
        Generate alternative query formats based on error analysis.
        
        Args:
            original_query: The original failed query
            error_analysis: Analysis from ErrorAnalyzer
            target_info: Original target information (optional)
            
        Returns:
            List of alternative query formats to try
        """
        # Determine query type and operation
        query_type = original_query.get('method', 'unknown')
        operation = original_query.get('operation', 'unknown')
        
        # If LLM client is available, use it for advanced reformulation
        if self.llm_client and target_info:
            alternatives = self._reformulate_with_llm(original_query, error_analysis, target_info)
            if alternatives:
                return alternatives
        
        # Fall back to rule-based reformulation
        return self._reformulate_with_rules(original_query, error_analysis, query_type, operation)
    
    def _initialize_patterns(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Initialize common query patterns for different operations.
        
        Returns:
            Dictionary of query patterns
        """
        return {
            'gene_name': [
                # Standard gene name query
                {
                    'method': 'gene_name',
                    'use_organism': True,
                    'case_sensitive': False,
                    'description': 'Standard gene name with organism'
                },
                # Case-sensitive gene name
                {
                    'method': 'gene_name',
                    'use_organism': True,
                    'case_sensitive': True,
                    'description': 'Case-sensitive gene name with organism'
                },
                # Gene name without organism
                {
                    'method': 'gene_name',
                    'use_organism': False,
                    'case_sensitive': False,
                    'description': 'Gene name without organism constraint'
                },
                # Alternative search with gene name in query
                {
                    'method': 'search',
                    'query_template': 'gene:{gene_name}',
                    'limit': 3,
                    'description': 'Search with gene name field'
                }
            ],
            'search': [
                # Basic keyword search
                {
                    'method': 'search',
                    'query_template': '{query}',
                    'limit': 5,
                    'description': 'Basic keyword search'
                },
                # Gene-specific search
                {
                    'method': 'search',
                    'query_template': 'gene:{query}',
                    'limit': 3,
                    'description': 'Gene-specific search'
                },
                # Protein name search
                {
                    'method': 'search',
                    'query_template': 'name:"{query}"',
                    'limit': 3,
                    'description': 'Protein name exact match'
                },
                # Combined gene and organism search
                {
                    'method': 'search',
                    'query_template': 'gene:{query} AND organism:"{organism}"',
                    'limit': 3,
                    'description': 'Gene and organism combined search'
                },
                # Accession-like query
                {
                    'method': 'search',
                    'query_template': 'accession:{query}',
                    'limit': 1,
                    'description': 'Accession lookup'
                }
            ],
            'accession': [
                # Standard accession lookup
                {
                    'method': 'accession',
                    'format_accession': False,
                    'description': 'Standard accession lookup'
                },
                # Accession lookup with formatting
                {
                    'method': 'accession',
                    'format_accession': True,
                    'description': 'Accession lookup with formatting'
                },
                # Fallback to search by accession
                {
                    'method': 'search',
                    'query_template': 'accession:{accession}',
                    'limit': 1,
                    'description': 'Search by accession'
                }
            ]
        }
    
    def _reformulate_with_llm(
        self, 
        original_query: Dict[str, Any], 
        error_analysis: Dict[str, Any],
        target_info: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Use LLM to generate alternative query formats.
        
        Args:
            original_query: The original failed query
            error_analysis: Analysis from ErrorAnalyzer
            target_info: Original target information
            
        Returns:
            List of alternative query formats to try
        """
        # Format the error analysis and query information for the LLM
        query_info = json.dumps(original_query, indent=2)
        error_info = json.dumps(error_analysis, indent=2)
        target_info_str = json.dumps(target_info, indent=2)
        
        # Construct the prompt
        prompt = f"""
        I need to reformulate a query that failed when retrieving protein information.
        
        Original target information:
        {target_info_str}
        
        Failed query:
        {query_info}
        
        Error analysis:
        {error_info}
        
        Please suggest 3 alternative query formulations that might succeed.
        Consider different search strategies, parameter formats, and alternative identifiers.
        Common approaches include:
        
        1. Using gene name with/without organism constraints
        2. Using different search fields (name, gene, keyword, etc.)
        3. Trying accession-based lookups if identifiers are available
        4. Using more general or more specific search terms
        5. Using alternative formatting for the same information
        
        Provide your response as a JSON array of query objects with the following structure:
        [
            {{
                "method": "search|gene_name|accession", // Method to use
                "query": "search query string", // For search method
                "gene_name": "gene name", // For gene_name method
                "organism": "organism name", // For methods that use organism
                "accession": "accession id", // For accession method
                "limit": 5, // Number of results to return (for search)
                "description": "Description of this reformulation strategy"
            }},
            // Additional reformulated queries...
        ]
        
        IMPORTANT: Ensure the JSON is valid and complete. Each object should have the appropriate fields for its method.
        """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.4)
            
            # Extract JSON from response
            json_start = response.find('[')
            json_end = response.rfind(']') + 1
            
            if json_start >= 0 and json_end > json_start:
                json_str = response[json_start:json_end]
                alternatives = json.loads(json_str)
                
                # Validate and clean up alternatives
                valid_alternatives = []
                for alt in alternatives:
                    if 'method' in alt:
                        # Ensure required fields are present based on method
                        if alt['method'] == 'search' and 'query' not in alt:
                            alt['query'] = target_info.get('name', '')
                        elif alt['method'] == 'gene_name' and 'gene_name' not in alt:
                            alt['gene_name'] = target_info.get('name', '')
                        elif alt['method'] == 'accession' and 'accession' not in alt:
                            # Skip if no accession is available
                            continue
                        
                        # Add description if missing
                        if 'description' not in alt:
                            alt['description'] = f"Reformulated {alt['method']} query"
                        
                        valid_alternatives.append(alt)
                
                return valid_alternatives
            
            # Fall back to rule-based if JSON extraction fails
            return []
            
        except Exception as e:
            # Fall back to rule-based on any error
            return []
    
    def _reformulate_with_rules(
        self, 
        original_query: Dict[str, Any], 
        error_analysis: Dict[str, Any],
        query_type: str,
        operation: str
    ) -> List[Dict[str, Any]]:
        """
        Use rule-based approach to generate alternative query formats.
        
        Args:
            original_query: The original failed query
            error_analysis: Analysis from ErrorAnalyzer
            query_type: Type of the original query
            operation: Operation being performed
            
        Returns:
            List of alternative query formats to try
        """
        alternatives = []
        
        # Get error type and probable cause
        error_type = error_analysis.get('error_type', 'unknown')
        probable_cause = error_analysis.get('probable_cause', 'unknown')
        
        # Extract relevant parameters from original query
        gene_name = original_query.get('gene_name', '')
        organism = original_query.get('organism', '')
        query = original_query.get('query', '')
        accession = original_query.get('accession', '')
        
        # If query type is in our patterns, use those
        if query_type in self.query_patterns:
            for pattern in self.query_patterns[query_type]:
                # Create a new query based on the pattern
                new_query = pattern.copy()
                
                # Fill in specific fields based on query type
                if query_type == 'gene_name':
                    new_query['gene_name'] = gene_name
                    if new_query.get('use_organism') and organism:
                        new_query['organism'] = organism
                elif query_type == 'search':
                    # Handle query templates
                    template = new_query.pop('query_template', '{query}')
                    if '{gene_name}' in template and gene_name:
                        query_str = template.replace('{gene_name}', gene_name)
                    elif '{query}' in template:
                        query_str = template.replace('{query}', query or gene_name)
                    elif '{accession}' in template and accession:
                        query_str = template.replace('{accession}', accession)
                    else:
                        query_str = gene_name or query
                    
                    # Handle organism in template
                    if '{organism}' in query_str and organism:
                        query_str = query_str.replace('{organism}', organism)
                    
                    new_query['query'] = query_str
                elif query_type == 'accession':
                    new_query['accession'] = accession
                    # Format accession if needed
                    if new_query.get('format_accession'):
                        # Ensure proper accession format (e.g., add UniProt prefix if missing)
                        if accession and not re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', accession):
                            # Try to convert to standard format
                            if re.match(r'^[A-Z0-9_]+$', accession):
                                new_query['accession'] = accession.upper()
                
                alternatives.append(new_query)
        
        # Add specialized alternatives based on error type
        if error_type == 'resource_not_found':
            # For not found errors, try broader search approaches
            alternatives.extend([
                {
                    'method': 'search',
                    'query': gene_name.split('_')[0] if '_' in gene_name else gene_name,
                    'limit': 5,
                    'description': 'Simplified gene name search'
                },
                {
                    'method': 'search',
                    'query': f"{gene_name} {organism.split()[0] if organism else ''}".strip(),
                    'limit': 5,
                    'description': 'Combined gene and partial organism search'
                }
            ])
        elif error_type == 'parsing_error':
            # For parsing errors, try simple formats
            alternatives.extend([
                {
                    'method': 'search',
                    'query': gene_name,
                    'format': 'txt',
                    'limit': 3,
                    'description': 'Simple text format search'
                }
            ])
        elif error_type == 'validation_error':
            # For validation errors, try different sources
            alternatives.extend([
                {
                    'method': 'search',
                    'query': f"reviewed:{gene_name}",
                    'limit': 3,
                    'description': 'Reviewed entries only'
                }
            ])
        
        # Ensure unique alternatives by creating a key for each
        unique_alternatives = {}
        for alt in alternatives:
            key = json.dumps({k: v for k, v in alt.items() if k != 'description'})
            unique_alternatives[key] = alt
        
        return list(unique_alternatives.values())
