#!/usr/bin/env python3
"""
Benchmark script to compare agentic vs. non-agentic Target Selector performance.
"""

import os
import sys
import json
import time
import argparse
from typing import Dict, List, Any
import pandas as pd
import matplotlib.pyplot as plt

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from agents.target_selector import TargetSelector
from utils.gemini_client import GeminiClient
from utils.uniprot_client import UniProtClient


def load_test_queries(file_path: str) -> List[str]:
    """Load test queries from a file."""
    with open(file_path, 'r') as f:
        queries = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return queries


def create_agent(use_agentic: bool) -> TargetSelector:
    """Create a Target Selector agent."""
    agent = TargetSelector(name=f"target_selector_{'agentic' if use_agentic else 'standard'}")
    
    if not use_agentic:
        # Disable agentic components
        agent.execute_with_retry = lambda func, *args, operation_name=None, **kwargs: func(*args, **kwargs)
        agent.error_analyzer = None
        agent.query_reformulator = None
        agent.sequence_validator = None
        agent.strategy_selector = None
    
    return agent


def run_benchmark(queries: List[str], use_agentic: bool) -> List[Dict[str, Any]]:
    """Run benchmark on a list of queries."""
    agent = create_agent(use_agentic)
    results = []
    
    for i, query in enumerate(queries):
        print(f"Processing query {i+1}/{len(queries)}: {query[:50]}...")
        
        try:
            # Measure execution time
            start_time = time.time()
            result = agent.process(query)
            end_time = time.time()
            
            # Calculate metrics
            execution_time = end_time - start_time
            num_proteins = len(result['parsed_data'].get('protein_targets', []))
            num_sequences = len(result['sequences'])
            num_configs = len(result['config_files'])
            
            # Record result
            results.append({
                'query': query,
                'execution_time': execution_time,
                'success': num_sequences > 0,
                'num_proteins': num_proteins,
                'num_sequences': num_sequences,
                'num_configs': num_configs,
                'agentic': use_agentic
            })
            
            print(f"  Success: {num_sequences}/{num_proteins} proteins retrieved in {execution_time:.2f}s")
            
        except Exception as e:
            print(f"  Error: {e}")
            
            # Record failure
            results.append({
                'query': query,
                'execution_time': 0,
                'success': False,
                'num_proteins': 0,
                'num_sequences': 0,
                'num_configs': 0,
                'error': str(e),
                'agentic': use_agentic
            })
    
    return results


def generate_report(standard_results: List[Dict[str, Any]], agentic_results: List[Dict[str, Any]], output_dir: str) -> None:
    """Generate benchmark report."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert to DataFrame
    df_standard = pd.DataFrame(standard_results)
    df_agentic = pd.DataFrame(agentic_results)
    
    # Calculate aggregate metrics
    standard_metrics = {
        'total_queries': len(df_standard),
        'success_rate': df_standard['success'].mean() * 100,
        'avg_execution_time': df_standard['execution_time'].mean(),
        'total_proteins': df_standard['num_proteins'].sum(),
        'total_sequences': df_standard['num_sequences'].sum(),
        'sequence_retrieval_rate': df_standard['num_sequences'].sum() / df_standard['num_proteins'].sum() * 100 if df_standard['num_proteins'].sum() > 0 else 0
    }
    
    agentic_metrics = {
        'total_queries': len(df_agentic),
        'success_rate': df_agentic['success'].mean() * 100,
        'avg_execution_time': df_agentic['execution_time'].mean(),
        'total_proteins': df_agentic['num_proteins'].sum(),
        'total_sequences': df_agentic['num_sequences'].sum(),
        'sequence_retrieval_rate': df_agentic['num_sequences'].sum() / df_agentic['num_proteins'].sum() * 100 if df_agentic['num_proteins'].sum() > 0 else 0
    }
    
    # Generate plots
    
    # Success rate comparison
    plt.figure(figsize=(10, 6))
    plt.bar(['Standard', 'Agentic'], [standard_metrics['success_rate'], agentic_metrics['success_rate']])
    plt.title('Success Rate Comparison')
    plt.ylabel('Success Rate (%)')
    plt.ylim(0, 100)
    plt.savefig(os.path.join(output_dir, 'success_rate.png'))
    
    # Execution time comparison
    plt.figure(figsize=(10, 6))
    plt.bar(['Standard', 'Agentic'], [standard_metrics['avg_execution_time'], agentic_metrics['avg_execution_time']])
    plt.title('Average Execution Time Comparison')
    plt.ylabel('Execution Time (s)')
    plt.savefig(os.path.join(output_dir, 'execution_time.png'))
    
    # Sequence retrieval rate comparison
    plt.figure(figsize=(10, 6))
    plt.bar(['Standard', 'Agentic'], [standard_metrics['sequence_retrieval_rate'], agentic_metrics['sequence_retrieval_rate']])
    plt.title('Sequence Retrieval Rate Comparison')
    plt.ylabel('Sequence Retrieval Rate (%)')
    plt.ylim(0, 100)
    plt.savefig(os.path.join(output_dir, 'retrieval_rate.png'))
    
    # Generate CSV report
    df_standard.to_csv(os.path.join(output_dir, 'standard_results.csv'), index=False)
    df_agentic.to_csv(os.path.join(output_dir, 'agentic_results.csv'), index=False)
    
    # Generate metrics report
    metrics_df = pd.DataFrame({
        'Metric': list(standard_metrics.keys()),
        'Standard': list(standard_metrics.values()),
        'Agentic': list(agentic_metrics.values()),
        'Improvement': [(agentic_metrics[k] - standard_metrics[k]) / standard_metrics[k] * 100 if standard_metrics[k] != 0 else float('inf') for k in standard_metrics.keys()]
    })
    
    metrics_df.to_csv(os.path.join(output_dir, 'metrics.csv'), index=False)
    
    # Generate summary report
    with open(os.path.join(output_dir, 'summary.md'), 'w') as f:
        f.write('# Agentic vs. Non-Agentic Target Selector Benchmark\n\n')
        
        f.write('## Summary Metrics\n\n')
        f.write('| Metric | Standard | Agentic | Improvement |\n')
        f.write('|--------|----------|---------|-------------|\n')
        
        for _, row in metrics_df.iterrows():
            metric = row['Metric']
            standard = row['Standard']
            agentic = row['Agentic']
            improvement = row['Improvement']
            
            # Format based on metric type
            if 'rate' in metric or 'success' in metric:
                f.write(f'| {metric} | {standard:.2f}% | {agentic:.2f}% | {improvement:+.2f}% |\n')
            elif 'time' in metric:
                f.write(f'| {metric} | {standard:.2f}s | {agentic:.2f}s | {improvement:+.2f}% |\n')
            else:
                f.write(f'| {metric} | {standard} | {agentic} | {improvement:+.2f}% |\n')
        
        f.write('\n## Conclusion\n\n')
        
        # Generate conclusion based on metrics
        if agentic_metrics['success_rate'] > standard_metrics['success_rate']:
            f.write(f"The agentic Target Selector shows a {(agentic_metrics['success_rate'] - standard_metrics['success_rate']):.2f}% higher success rate than the standard version. ")
        else:
            f.write(f"The standard Target Selector shows a {(standard_metrics['success_rate'] - agentic_metrics['success_rate']):.2f}% higher success rate than the agentic version. ")
            
        if agentic_metrics['sequence_retrieval_rate'] > standard_metrics['sequence_retrieval_rate']:
            f.write(f"It also retrieves {(agentic_metrics['sequence_retrieval_rate'] - standard_metrics['sequence_retrieval_rate']):.2f}% more sequences. ")
        else:
            f.write(f"However, it retrieves {(standard_metrics['sequence_retrieval_rate'] - agentic_metrics['sequence_retrieval_rate']):.2f}% fewer sequences. ")
            
        if agentic_metrics['avg_execution_time'] < standard_metrics['avg_execution_time']:
            f.write(f"The agentic version is {((standard_metrics['avg_execution_time'] - agentic_metrics['avg_execution_time']) / standard_metrics['avg_execution_time'] * 100):.2f}% faster.")
        else:
            f.write(f"The agentic version is {((agentic_metrics['avg_execution_time'] - standard_metrics['avg_execution_time']) / standard_metrics['avg_execution_time'] * 100):.2f}% slower.")
    
    print(f"Benchmark report generated in {output_dir}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Benchmark agentic vs. non-agentic Target Selector')
    parser.add_argument('--queries', type=str, default='tests/test_queries.txt',
                        help='Path to test queries file')
    parser.add_argument('--output-dir', type=str, default='test_results/benchmark',
                        help='Directory to save benchmark results')
    args = parser.parse_args()
    
    # Load test queries
    queries = load_test_queries(args.queries)
    print(f"Loaded {len(queries)} test queries")
    
    # Run benchmarks
    print("\nRunning benchmark with standard (non-agentic) Target Selector...")
    standard_results = run_benchmark(queries, use_agentic=False)
    
    print("\nRunning benchmark with agentic Target Selector...")
    agentic_results = run_benchmark(queries, use_agentic=True)
    
    # Generate report
    generate_report(standard_results, agentic_results, args.output_dir)


if __name__ == '__main__':
    main()
