"""
Main entry point for F.A.D.E drug discovery pipeline using LangGraph.
"""

import sys
import argparse
import json
from pathlib import Path

from fade.workflows.drug_discovery import run_drug_discovery_pipeline
from fade.config import config
from fade.utils import setup_logging, get_logger

logger = get_logger("main")


def main():
    """Main entry point for CLI usage."""
    parser = argparse.ArgumentParser(
        description="F.A.D.E - Fully Agentic Drug Engine (LangGraph Implementation)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py "Find inhibitors for KRAS G12C"
  python main.py "Design molecules targeting EGFR T790M mutation"
  python main.py "Create blood-brain barrier penetrant drugs for tau protein"
        """
    )
    
    parser.add_argument(
        "query",
        nargs="?",
        default="Find inhibitors for KRAS G12C",
        help="Natural language query describing the drug discovery task"
    )
    
    parser.add_argument(
        "--user-id",
        help="Optional user identifier for tracking"
    )
    
    parser.add_argument(
        "--output-format",
        choices=["text", "json"],
        default="text",
        help="Output format for results"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    parser.add_argument(
        "--no-rich",
        action="store_true",
        help="Disable rich console output"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(level=log_level, use_rich=not args.no_rich)
    
    # Validate configuration
    try:
        config.validate()
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        logger.info("Please check your .env file")
        sys.exit(1)
    
    # Run the LangGraph pipeline
    logger.info(f"Starting F.A.D.E pipeline with LangGraph")
    logger.info(f"Query: {args.query}")
    
    final_state = run_drug_discovery_pipeline(args.query, args.user_id)
    
    # Output results
    if args.output_format == "json":
        # Convert to JSON-serializable format
        status = "completed"
        if final_state.get("error"):
            status = "needs_refinement" if final_state.get("error_type") == "no_structures_found" else "failed"
        
        output = {
            "run_id": final_state["run_id"],
            "query": final_state["query"],
            "status": status,
            "error": final_state.get("error"),
            "error_type": final_state.get("error_type"),
            "current_step": final_state.get("current_step"),
            "suggested_queries": final_state.get("suggested_queries"),
            "guidance": final_state.get("guidance"),
            "target": {
                "name": final_state.get("target_info", {}).get("protein_name"),
                "uniprot_id": final_state.get("target_info", {}).get("uniprot_id"),
                "gene_name": final_state.get("target_info", {}).get("gene_name"),
                "mutations": final_state.get("target_info", {}).get("mutations")
            } if final_state.get("target_info") else None,

            "messages": [msg.content if hasattr(msg, 'content') else str(msg) 
                        for msg in final_state.get("messages", [])]
        }
        print(json.dumps(output, indent=2, default=str))
        
    else:  # text output
        print("\n" + "="*80)
        print("F.A.D.E DRUG DISCOVERY RESULTS")
        print("="*80)
        
        if final_state.get("error"):
            # Check if this is guidance rather than a hard error
            if final_state.get("error_type") in ["no_structures_found", "uniprot_not_found"]:
                # Both RCSB and UniProt guidance use similar format
                print("\n" + "="*80)
                if final_state.get("error_type") == "uniprot_not_found":
                    # UniProt guidance is already well-formatted
                    print(final_state['error'])
                else:
                    # RCSB guidance
                    print("💡 QUERY REFINEMENT NEEDED")
                    print("="*80)
                    print(f"\n{final_state['error']}")
                    
                    # Show suggested queries if available
                    if final_state.get("suggested_queries"):
                        print("\n📝 Try one of these queries:")
                        for i, suggestion in enumerate(final_state["suggested_queries"], 1):
                            print(f"   {i}. python main.py \"{suggestion}\"")
                    
                    # Show detailed guidance if available
                    if final_state.get("guidance"):
                        guidance = final_state["guidance"]
                        if guidance.get("refinement") and guidance["refinement"].get("confidence"):
                            confidence = guidance["refinement"]["confidence"]
                            print(f"\n🎯 Confidence in suggestions: {confidence:.0%}")
            else:
                print(f"\n❌ Pipeline Failed: {final_state['error']}")
                print(f"   Step: {final_state.get('current_step', 'unknown')}")
        
        elif final_state.get("target_info"):
            target = final_state["target_info"]
            print(f"\n✓ Target Identified:")
            print(f"  • Protein: {target.get('protein_name', 'Unknown')}")
            print(f"  • UniProt ID: {target.get('uniprot_id', 'N/A')}")
            print(f"  • Gene: {target.get('gene_name', 'N/A')}")
            
            if target.get("mutations"):
                print(f"  • Mutations: {', '.join(target['mutations'])}")
            
            if target.get("sequence_length"):
                print(f"  • Sequence Length: {target['sequence_length']} amino acids")
            

            if target.get("existing_structures"):
                structures = target["existing_structures"]
                print(f"\n✓ Existing Structures: {len(structures)} PDB entries")
                for struct in structures[:3]:
                    print(f"  • {struct['pdb_id']}: {struct.get('title', 'N/A')[:60]}...")
                    if struct.get("resolution"):
                        print(f"    Resolution: {struct['resolution']} Å")
            
            # Show pocket information
            if final_state.get("pockets"):
                pockets = final_state["pockets"]
                print(f"\n✓ Binding Pockets Detected: {len(pockets)}")
                for pocket in pockets[:3]:
                    print(f"  • {pocket['pocket_id']}:")
                    print(f"    Druggability: {pocket.get('druggability_score', 0):.2f}")
                    print(f"    Volume: {pocket.get('volume', 0):.1f} Ų")
                    if pocket.get("description"):
                        print(f"    Description: {pocket['description']}")
            
            if final_state.get("selected_pocket"):
                pocket = final_state["selected_pocket"]
                print(f"\n✓ Selected Target Pocket: {pocket['pocket_id']}")
                if final_state.get("pocket_selection_rationale"):
                    print(f"  Rationale: {final_state['pocket_selection_rationale']}")
        
        else:
            print("\n⚠ No target information found")
            print(f"Current step: {final_state.get('current_step', 'unknown')}")
        
        # Show workflow messages
        if args.verbose and final_state.get("messages"):
            print("\n" + "-"*40)
            print("Workflow Messages:")
            for msg in final_state["messages"]:
                content = msg.content if hasattr(msg, 'content') else str(msg)
                print(f"  • {content}")
    
    # Exit with appropriate code
    sys.exit(0 if not final_state.get("error") else 1)


if __name__ == "__main__":
    main()
