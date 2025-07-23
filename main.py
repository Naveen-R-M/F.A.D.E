"""
Main entry point for the F.A.D.E system.
"""

import argparse
import os
import logging
import json
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Setup logging
logging.basicConfig(
    level=os.getenv("LOG_LEVEL", "INFO"),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/fade.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("fade")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="F.A.D.E: Framework for Agentic Drug Exploration")
    
    parser.add_argument(
        "--query", 
        type=str,
        help="Natural language query for drug discovery"
    )
    
    parser.add_argument(
        "--config", 
        type=str,
        default="configs/default.yaml",
        help="Path to configuration file"
    )
    
    parser.add_argument(
        "--batch-file", 
        type=str,
        help="File containing multiple queries to process in batch"
    )
    
    parser.add_argument(
        "--output-dir", 
        type=str,
        default="data/outputs",
        help="Directory for output results"
    )
    
    parser.add_argument(
        "--verbose", 
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser.parse_args()

def run_pipeline(query, config_path, output_dir):
    """
    Run the F.A.D.E pipeline for a single query.
    
    Args:
        query (str): Natural language query
        config_path (str): Path to config file
        output_dir (str): Directory to store results
    
    Returns:
        dict: Results of the pipeline
    """
    logger.info(f"Starting pipeline for query: {query}")
    
    # TODO: Implement the full pipeline
    # This is a placeholder for the actual implementation
    
    # 1. Import the workflow orchestrator
    from workflows.pipeline_runner import run_workflow
    
    # 2. Run the workflow
    results = run_workflow(query, config_path, output_dir)
    
    # 3. Return results
    return results

def main():
    """Main function to run the F.A.D.E system."""
    args = parse_arguments()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.batch_file:
        # Process multiple queries from a file
        with open(args.batch_file, 'r') as f:
            queries = [line.strip() for line in f if line.strip()]
        
        results = []
        for i, query in enumerate(queries):
            logger.info(f"Processing query {i+1}/{len(queries)}")
            result = run_pipeline(query, args.config, args.output_dir)
            results.append(result)
            
        # Save batch results
        with open(os.path.join(args.output_dir, "batch_results.json"), 'w') as f:
            json.dump(results, f, indent=2)
            
    elif args.query:
        # Process a single query
        result = run_pipeline(args.query, args.config, args.output_dir)
        
        # Save result
        with open(os.path.join(args.output_dir, "result.json"), 'w') as f:
            json.dump(result, f, indent=2)
    else:
        logger.error("No query provided. Use --query or --batch-file")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
