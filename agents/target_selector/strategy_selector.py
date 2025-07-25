"""
Strategy selector component for the Target Selector agent.

This component selects appropriate search strategies based on context and history.
"""

from typing import Any, Dict, List, Optional, Tuple
import json
from utils.gemini_client import GeminiClient

class SearchStrategySelector:
    """
    Selects appropriate search strategy based on context and history.
    """
    
    def __init__(self, llm_client: Optional[GeminiClient] = None) -> None:
        """
        Initialize the search strategy selector.
        
        Args:
            llm_client: Client for the LLM (e.g., GeminiClient)
        """
        self.llm_client = llm_client
        self.default_strategies = self._initialize_default_strategies()
    
    def select_strategy(
        self, 
        protein_target: Dict[str, Any], 
        previous_attempts: List[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Select search strategy based on target and previous attempts.
        
        Args:
            protein_target: Target protein information
            previous_attempts: Record of previous search attempts
            
        Returns:
            Search strategy configuration
        """
        if previous_attempts is None:
            previous_attempts = []
        
        # If LLM client is available, use it for advanced strategy selection
        if self.llm_client and previous_attempts:
            llm_strategy = self._select_with_llm(protein_target, previous_attempts)
            if llm_strategy:
                return llm_strategy
        
        # Fall back to rule-based strategy selection
        return self._select_with_rules(protein_target, previous_attempts)
    
    def _initialize_default_strategies(self) -> List[Dict[str, Any]]:
        """
        Initialize default search strategies in order of specificity.
        
        Returns:
            List of default strategies
        """
        return [
            # Strategy 1: Direct gene name lookup with organism constraint
            {
                'method': 'gene_name',
                'description': 'Direct gene name lookup with organism constraint',
                'requires': ['name', 'organism'],
                'success_rate': 0.9,
                'specificity': 'high'
            },
            # Strategy 2: Direct gene name lookup without organism constraint
            {
                'method': 'gene_name',
                'description': 'Direct gene name lookup without organism constraint',
                'requires': ['name'],
                'success_rate': 0.8,
                'specificity': 'medium'
            },
            # Strategy 3: Accession lookup
            {
                'method': 'accession',
                'description': 'Direct accession lookup',
                'requires': ['accession'],
                'success_rate': 0.95,
                'specificity': 'very high'
            },
            # Strategy 4: Search by gene name and organism
            {
                'method': 'search',
                'description': 'Search by gene name and organism',
                'query_template': 'gene:{gene_name} AND organism:"{organism}"',
                'requires': ['name', 'organism'],
                'success_rate': 0.85,
                'specificity': 'high',
                'limit': 3
            },
            # Strategy 5: Search by gene name
            {
                'method': 'search',
                'description': 'Search by gene name',
                'query_template': 'gene:{gene_name}',
                'requires': ['name'],
                'success_rate': 0.7,
                'specificity': 'medium',
                'limit': 5
            },
            # Strategy 6: Search by protein name
            {
                'method': 'search',
                'description': 'Search by protein name',
                'query_template': 'name:"{name}"',
                'requires': ['name'],
                'success_rate': 0.6,
                'specificity': 'medium',
                'limit': 5
            },
            # Strategy 7: Search by keyword
            {
                'method': 'search',
                'description': 'General keyword search',
                'query_template': '{name}',
                'requires': ['name'],
                'success_rate': 0.5,
                'specificity': 'low',
                'limit': 10
            }
        ]
    
    def _select_with_llm(
        self, 
        protein_target: Dict[str, Any], 
        previous_attempts: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """
        Use LLM to select an appropriate search strategy.
        
        Args:
            protein_target: Target protein information
            previous_attempts: Record of previous search attempts
            
        Returns:
            Selected search strategy or None if LLM fails
        """
        # Format the target and attempts information for the LLM
        target_info = json.dumps(protein_target, indent=2)
        attempts_info = json.dumps(previous_attempts, indent=2)
        strategies_info = json.dumps(self.default_strategies, indent=2)
        
        # Construct the prompt
        prompt = f"""
        I need to select an appropriate search strategy for retrieving protein information.
        
        Target protein information:
        {target_info}
        
        Previous search attempts:
        {attempts_info}
        
        Available search strategies:
        {strategies_info}
        
        Please select the most appropriate search strategy based on the target information and previous attempts.
        Consider:
        1. What strategies have already been tried and failed
        2. What information is available in the target data
        3. The specificity and success rate of each strategy
        4. How to adapt the strategy based on previous errors
        
        You can use one of the provided strategies or create a custom one.
        Provide your response as a JSON object with the following structure:
        {{
            "method": "search|gene_name|accession", // Method to use
            "description": "Description of this strategy",
            "reasoning": "Explanation of why this strategy was selected",
            // For search method:
            "query": "search query string",
            "limit": 5,
            // For gene_name method:
            "gene_name": "gene name",
            "organism": "organism name", // Optional
            // For accession method:
            "accession": "accession id"
        }}
        
        IMPORTANT: Ensure the JSON is valid and complete. The strategy should include all required fields for the selected method.
        """
        
        try:
            # Get response from LLM
            response = self.llm_client.generate_text(prompt, temperature=0.3)
            
            # Extract JSON from response
            json_start = response.find('{')
            json_end = response.rfind('}') + 1
            
            if json_start >= 0 and json_end > json_start:
                json_str = response[json_start:json_end]
                strategy = json.loads(json_str)
                
                # Validate strategy
                if 'method' in strategy:
                    # Ensure required fields are present based on method
                    if strategy['method'] == 'search' and 'query' not in strategy:
                        if 'query_template' in strategy:
                            # Apply template
                            query = strategy['query_template']
                            for key, value in protein_target.items():
                                if isinstance(value, str):
                                    query = query.replace(f"{{{key}}}", value)
                            strategy['query'] = query
                        else:
                            strategy['query'] = protein_target.get('name', '')
                    elif strategy['method'] == 'gene_name' and 'gene_name' not in strategy:
                        strategy['gene_name'] = protein_target.get('name', '')
                    elif strategy['method'] == 'accession' and 'accession' not in strategy:
                        if 'accession' in protein_target:
                            strategy['accession'] = protein_target['accession']
                        else:
                            # Fall back to gene name if no accession is available
                            strategy['method'] = 'gene_name'
                            strategy['gene_name'] = protein_target.get('name', '')
                            strategy['description'] = 'Fallback to gene name (no accession available)'
                    
                    # Add description if missing
                    if 'description' not in strategy:
                        strategy['description'] = f"LLM-selected {strategy['method']} strategy"
                    
                    return strategy
            
            # Fall back to rule-based if JSON extraction fails
            return None
            
        except Exception as e:
            # Fall back to rule-based on any error
            return None
    
    def _select_with_rules(
        self, 
        protein_target: Dict[str, Any], 
        previous_attempts: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Use rule-based approach to select an appropriate search strategy.
        
        Args:
            protein_target: Target protein information
            previous_attempts: Record of previous search attempts
            
        Returns:
            Selected search strategy
        """
        # Track attempted strategies
        attempted_methods = set()
        attempted_queries = set()
        
        for attempt in previous_attempts:
            strategy = attempt.get('strategy', {})
            method = strategy.get('method')
            
            if method:
                attempted_methods.add(method)
            
            if method == 'search':
                query = strategy.get('query')
                if query:
                    attempted_queries.add(query)
            elif method == 'gene_name':
                gene_name = strategy.get('gene_name')
                organism = strategy.get('organism')
                if gene_name and organism:
                    attempted_queries.add(f"gene_name:{gene_name}:organism:{organism}")
                elif gene_name:
                    attempted_queries.add(f"gene_name:{gene_name}")
            elif method == 'accession':
                accession = strategy.get('accession')
                if accession:
                    attempted_queries.add(f"accession:{accession}")
        
        # Check what information is available in the target
        has_name = 'name' in protein_target and protein_target['name']
        has_organism = 'organism' in protein_target and protein_target['organism']
        has_accession = 'accession' in protein_target and protein_target['accession']
        
        # Go through strategies in order of specificity
        for strategy in self.default_strategies:
            # Check if all required fields are available
            has_all_required = all(protein_target.get(req) for req in strategy.get('requires', []))
            
            if not has_all_required:
                continue
            
            # Check if this exact strategy has been attempted
            method = strategy['method']
            
            if method == 'gene_name':
                gene_name = protein_target.get('name')
                organism = protein_target.get('organism') if 'organism' in strategy.get('requires', []) else None
                
                if organism:
                    strategy_key = f"gene_name:{gene_name}:organism:{organism}"
                else:
                    strategy_key = f"gene_name:{gene_name}"
                
                if strategy_key in attempted_queries:
                    continue
                
                # Create gene name strategy
                return {
                    'method': 'gene_name',
                    'gene_name': gene_name,
                    'organism': organism,
                    'description': strategy['description']
                }
                
            elif method == 'accession' and has_accession:
                accession = protein_target.get('accession')
                strategy_key = f"accession:{accession}"
                
                if strategy_key in attempted_queries:
                    continue
                
                # Create accession strategy
                return {
                    'method': 'accession',
                    'accession': accession,
                    'description': strategy['description']
                }
                
            elif method == 'search':
                # Apply query template
                template = strategy.get('query_template', '{name}')
                query = template
                
                for key, value in protein_target.items():
                    if isinstance(value, str):
                        query = query.replace(f"{{{key}}}", value)
                
                if query in attempted_queries:
                    continue
                
                # Create search strategy
                return {
                    'method': 'search',
                    'query': query,
                    'limit': strategy.get('limit', 5),
                    'description': strategy['description']
                }
        
        # If all standard strategies have been attempted, create a more creative search
        if has_name:
            name = protein_target.get('name')
            
            # Try different variations of the name
            if '_' in name:
                # Try with underscores removed
                clean_name = name.replace('_', ' ')
                return {
                    'method': 'search',
                    'query': clean_name,
                    'limit': 10,
                    'description': 'Search with cleaned name (underscores removed)'
                }
            
            if has_organism:
                organism = protein_target.get('organism')
                organism_parts = organism.split()
                
                if len(organism_parts) > 1:
                    # Try with only genus
                    genus = organism_parts[0]
                    return {
                        'method': 'search',
                        'query': f"{name} {genus}",
                        'limit': 10,
                        'description': 'Search with name and genus only'
                    }
            
            # Try a general search with just the name
            return {
                'method': 'search',
                'query': name,
                'limit': 20,
                'description': 'Broad search with name only'
            }
        
        # Last resort: very broad search
        return {
            'method': 'search',
            'query': protein_target.get('name', 'unknown'),
            'limit': 30,
            'description': 'Last resort broad search'
        }
