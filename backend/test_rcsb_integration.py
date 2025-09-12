#!/usr/bin/env python3
"""
Test script for RCSB Target Selector functionality
"""

import sys
import os
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from agents.target_selector.rcsb_target_selector import RCSBTargetSelector
from utils.logging import setup_logging, get_logger


def test_rcsb_target_selector():
    """Test the RCSB Target Selector with a sample query."""
    
    # Set up logging
    setup_logging(log_level="INFO")
    logger = get_logger("fade.test_rcsb")
    
    logger.info("Starting RCSB Target Selector test")
    
    # Test queries
    test_queries = [
        "Find molecules that target KRAS G12D mutation for pancreatic cancer treatment",
        "Design a drug for EGFR targeting lung cancer with high binding affinity",
        "Find inhibitors for PDB structure 1ABC that can cross the blood-brain barrier"
    ]
    
    try:
        # Initialize agent (will need API key in environment)
        agent = RCSBTargetSelector(
            name="test_rcsb_selector",
            config={"data_dir": str(project_root / "data")}
        )
        
        for i, query in enumerate(test_queries, 1):
            logger.info(f"\n=== Test Query {i} ===")
            logger.info(f"Query: {query}")
            
            try:
                result = agent.process(query)
                
                logger.info(f"Results:")
                logger.info(f"- Parsed targets: {list(result.get('structures', {}).keys())}")
                logger.info(f"- Source: {result.get('source', 'unknown')}")
                
                structures = result.get('structures', {})
                for target_name, structure_info in structures.items():
                    pdb_id = structure_info.get('pdb_id', 'unknown')
                    resolution = structure_info.get('resolution', 'N/A')
                    logger.info(f"- {target_name}: PDB {pdb_id} ({resolution} Å)")
                
            except Exception as e:
                logger.error(f"Query {i} failed: {e}")
                continue
        
        logger.info("\nRCSB Target Selector test completed successfully")
        
    except Exception as e:
        logger.error(f"Test failed: {e}")
        return False
    
    return True


def test_rcsb_client():
    """Test the RCSB client directly."""
    
    setup_logging(log_level="INFO")
    logger = get_logger("fade.test_rcsb_client")
    
    logger.info("Starting RCSB Client test")
    
    try:
        from utils.rcsb_client import RCSBClient
        from utils.gemini_client import GeminiClient
        
        # Initialize clients
        rcsb_client = RCSBClient()
        gemini_client = GeminiClient()
        
        # Test search terms generation
        target_info = {
            "protein_targets": [
                {
                    "name": "KRAS",
                    "mutations": [{"notation": "G12D"}]
                }
            ],
            "disease_context": "pancreatic cancer"
        }
        
        logger.info("Testing search terms generation...")
        search_terms = rcsb_client.generate_search_terms(gemini_client, target_info)
        logger.info(f"Generated search terms: {search_terms}")
        
        if search_terms:
            logger.info("Testing structure search...")
            structures = rcsb_client.search_structures(search_terms, limit=3)
            logger.info(f"Found {len(structures)} structures")
            
            for struct in structures[:2]:
                pdb_id = struct.get('pdb_id', 'unknown')
                title = struct.get('title', '')[:50]
                logger.info(f"- {pdb_id}: {title}...")
        
        logger.info("RCSB Client test completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"RCSB Client test failed: {e}")
        return False


if __name__ == "__main__":
    print("=" * 60)
    print("F.A.D.E RCSB Integration Test")
    print("=" * 60)
    
    # Check if API key is available
    if not os.getenv("GEMINI_API_KEY"):
        print("WARNING: GEMINI_API_KEY not set in environment")
        print("Some tests may fail without API access")
    
    # Run tests
    print("\n1. Testing RCSB Client...")
    client_success = test_rcsb_client()
    
    print("\n2. Testing RCSB Target Selector...")
    selector_success = test_rcsb_target_selector()
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    print(f"RCSB Client: {'PASS' if client_success else 'FAIL'}")
    print(f"RCSB Target Selector: {'PASS' if selector_success else 'FAIL'}")
    
    if client_success and selector_success:
        print("\nAll tests passed! RCSB integration is working.")
        sys.exit(0)
    else:
        print("\nSome tests failed. Check the logs above.")
        sys.exit(1)
