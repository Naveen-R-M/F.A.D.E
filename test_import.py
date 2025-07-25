#!/usr/bin/env python3
"""
Test script to verify imports are working correctly.
"""

try:
    # Test importing the TargetSelector
    from agents.target_selector import TargetSelector
    print("✅ Successfully imported TargetSelector")
    
    # Test importing the components
    from agents.target_selector import ErrorAnalyzer, QueryReformulator, SequenceValidator, SearchStrategySelector
    print("✅ Successfully imported agentic components")
    
    # Test importing the AgenticMixin
    from agents.base.agentic_mixin import AgenticMixin
    print("✅ Successfully imported AgenticMixin")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    
print("\nImport test completed!")
