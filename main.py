#!/usr/bin/env python3
"""
Main entry point for F.A.D.E
"""

import os
import sys
import json
import argparse
from typing import Any, Dict, Optional

from dotenv import load_dotenv
from utils.logging import setup_logging, get_logger
from utils.summary import SummaryLogger
from agents.target_selector import TargetSelector
from agents.structure_predictor.structure_predictor import StructurePredictor


def setup_logging_with_output_dir(log_level: str = "INFO", output_dir: Optional[str] = None) -> None:
    """
    Set up logging for the application with output directory integration.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        output_dir: Optional output directory where results will be saved.
                   If provided, logs will also be stored in this directory.
    """
    # Determine logs directory based on output directory
    logs_dir = "logs"
    if output_dir:
        logs_dir = os.path.join(output_dir, "logs")
    
    # Set up structured logging
    setup_logging(log_level=log_level, logs_dir=logs_dir)
    
    # Get the main logger
    logger = get_logger("fade.main")
    logger.info(f"Logging initialized with logs directory: {logs_dir}")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="F.A.D.E: Fully Agentic Drug Engine")
    
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
    
    parser.add_argument(
        "--skip-structure-prediction",
        action="store_true",
        help="Skip the structure prediction step"
    )
    
    return parser.parse_args()


def process_query(query: str, output_dir: Optional[str] = None, skip_structure_prediction: bool = False) -> Dict[str, Any]:
    """
    Process a natural language query through the F.A.D.E pipeline.
    
    Args:
        query: Natural language query describing the drug discovery goal.
        output_dir: Optional directory where results will be saved.
        skip_structure_prediction: If True, skip the structure prediction step.
        
    Returns:
        Dictionary containing the results.
    """
    # Get logger
    logger = get_logger("fade.main")
    
    # Initialize summary logger if output directory is provided
    summary_logger = None
    if output_dir:
        summary_logger = SummaryLogger(output_dir)
        summary_logger.log_query(query)
    
    # Initialize the Target Selector agent
    target_selector = TargetSelector()
    
    # Process the query with Target Selector
    logger.info(f"Processing query with Target Selector: {query}")
    target_selector_results = target_selector.process(query)
    
    # Log summary of target selector results
    if summary_logger:
        summary_logger.log_agent_result("target_selector", target_selector_results)
    
    # Initialize results dictionary
    results = {
        "target_selector": target_selector_results
    }
    
    # Continue with Structure Predictor if not skipped
    if not skip_structure_prediction:
        # Initialize the Structure Predictor agent
        structure_predictor = StructurePredictor()
        
        # Process results with Structure Predictor
        logger.info("Processing with Structure Predictor")
        structure_predictor_input = {
            "sequences": target_selector_results.get("sequences", {}),
            "job_configs": target_selector_results.get("config_files", {})
        }
        structure_predictor_results = structure_predictor.process(structure_predictor_input)
        
        # Add structure predictor results to the overall results
        results["structure_predictor"] = structure_predictor_results
        
        # Log summary of structure predictor results
        if summary_logger:
            summary_logger.log_agent_result("structure_predictor", structure_predictor_results)
    
    # Save results if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        results_path = os.path.join(output_dir, "results.json")
        
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)
        logger.info(f"Results saved to: {results_path}")
        
        # Log final summary
        if summary_logger:
            summary_logger.log_final_summary(results)
    
    return results


def main() -> None:
    """Main entry point."""
    # Load environment variables from .env file
    load_dotenv()
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set up logging with output directory integration
    setup_logging_with_output_dir(args.log_level, args.output_dir)
    
    # Get logger
    logger = get_logger("fade.main")
    
    # Process query or batch file
    if args.query:
        # Process a single query
        logger.info(f"Processing query: {args.query}")
        results = process_query(args.query, args.output_dir, args.skip_structure_prediction)
        logger.info(f"Processing completed. Results saved to: {args.output_dir}")
        
    elif args.batch_file:
        # Process multiple queries from a batch file
        if not os.path.exists(args.batch_file):
            logger.error(f"Batch file not found: {args.batch_file}")
            sys.exit(1)
            
        logger.info(f"Processing batch file: {args.batch_file}")
        
        with open(args.batch_file, "r") as f:
            queries = [line.strip() for line in f if line.strip() and not line.startswith("#")]
            
        for i, query in enumerate(queries):
            query_output_dir = os.path.join(args.output_dir, f"query_{i+1}")
            logger.info(f"Processing query {i+1}/{len(queries)}: {query}")
            process_query(query, query_output_dir, args.skip_structure_prediction)
            
        logger.info(f"Batch processing completed. Results saved to: {args.output_dir}")
        
    else:
        # No query or batch file provided
        logger.error("No query or batch file provided. Use --query or --batch-file.")
        sys.exit(1)


if __name__ == "__main__":
    main()
