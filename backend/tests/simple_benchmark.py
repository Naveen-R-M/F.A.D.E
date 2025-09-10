#!/usr/bin/env python3
"""
Simple benchmark script to verify both agentic and non-agentic Target Selector work.
"""

import os
import sys
import time

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from agents.target_selector import TargetSelector

def main():
    # Test query
    query = "Find molecules targeting KRAS G12D with good BBB permeability"
    
    print("Testing standard (non-agentic) Target Selector...")
    # Create standard agent
    standard_agent = TargetSelector(name="standard_target_selector")
    
    # Disable agentic components
    standard_agent.execute_with_retry = lambda func, *args, operation_name=None, **kwargs: func(*args, **kwargs)
    standard_agent.error_analyzer = None
    standard_agent.query_reformulator = None
    standard_agent.sequence_validator = None
    standard_agent.strategy_selector = None
    
    # Process query
    start_time = time.time()
    result = standard_agent.process(query)
    end_time = time.time()
    
    # Print results
    print(f"Standard agent execution time: {end_time - start_time:.2f}s")
    print(f"Number of proteins parsed: {len(result['parsed_data'].get('protein_targets', []))}")
    print(f"Number of sequences retrieved: {len(result['sequences'])}")
    print(f"Number of config files generated: {len(result['config_files'])}")
    print()
    
    print("Testing agentic Target Selector...")
    # Create agentic agent
    agentic_agent = TargetSelector(name="agentic_target_selector")
    
    # Process query
    start_time = time.time()
    result = agentic_agent.process(query)
    end_time = time.time()
    
    # Print results
    print(f"Agentic agent execution time: {end_time - start_time:.2f}s")
    print(f"Number of proteins parsed: {len(result['parsed_data'].get('protein_targets', []))}")
    print(f"Number of sequences retrieved: {len(result['sequences'])}")
    print(f"Number of config files generated: {len(result['config_files'])}")
    
    print("\nBenchmark completed successfully!")

if __name__ == "__main__":
    main()
