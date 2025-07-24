#!/usr/bin/env python3
"""
Main entry point for F.A.D.E framework.
"""

import os
import sys
import json
import argparse
import logging
from typing import Any, Dict, Optional

from dotenv import load_dotenv
from agents.target_selector import TargetSelector


def setup_logging(log_level: str = "INFO") -> None:
    """
    Set up logging for the application.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")
    
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("fade.log")
        ]
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="F.A.D.E: Framework for Agentic Drug Exploration")
    
    parser.add_argument(
        "--query", "-q",
        type=str,
        help="Natural language query describing the drug discovery goal"
    )
    
    parser.add_argument(
        "--batch-file",
        type=str,
        help="Path to a file containing multiple queries (one per line)"
    )
    
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default="results",
        help="Directory where results will be saved"
    )
    
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level"
    )
    
    return parser.parse_args()


def process_query(query: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    Process a natural language query through the F.A.D.E pipeline.
    
    Args:
        query: Natural language query describing the drug discovery goal.
        output_dir: Optional directory where results will be saved.
        
    Returns:
        Dictionary containing the results.
    """
    # Initialize the Target Selector agent
    target_selector = TargetSelector()
    
    # Process the query
    results = target_selector.process(query)
    
    # Save results if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        results_path = os.path.join(output_dir, "results.json")
        
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)
    
    return results


def main() -> None:
    """Main entry point."""
    # Load environment variables from .env file
    load_dotenv()
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set up logging
    setup_logging(args.log_level)
    
    # Process query or batch file
    if args.query:
        # Process a single query
        logging.info("Processing query: %s", args.query)
        results = process_query(args.query, args.output_dir)
        logging.info("Processing completed. Results saved to: %s", args.output_dir)
        
    elif args.batch_file:
        # Process multiple queries from a batch file
        if not os.path.exists(args.batch_file):
            logging.error("Batch file not found: %s", args.batch_file)
            sys.exit(1)
            
        logging.info("Processing batch file: %s", args.batch_file)
        
        with open(args.batch_file, "r") as f:
            queries = [line.strip() for line in f if line.strip()]
            
        for i, query in enumerate(queries):
            query_output_dir = os.path.join(args.output_dir, f"query_{i+1}")
            logging.info("Processing query %d/%d: %s", i+1, len(queries), query)
            process_query(query, query_output_dir)
            
        logging.info("Batch processing completed. Results saved to: %s", args.output_dir)
        
    else:
        # No query or batch file provided
        logging.error("No query or batch file provided. Use --query or --batch-file.")
        sys.exit(1)


if __name__ == "__main__":
    main()
