#!/bin/bash

# Run unit tests for F.A.D.E

# Set up environment
cd $(dirname $0)/..
export PYTHONPATH=$PWD:$PYTHONPATH

# Run all unit tests
if [ "$1" == "--all" ]; then
    python -m unittest discover -s tests/unit -p "test_*.py"
    exit $?
fi

# Run specific test if provided
if [ ! -z "$1" ]; then
    python -m unittest tests/unit/$1
    exit $?
fi

# Default: Run only the agentic tests
python -m unittest tests/unit/test_agentic_mixin.py tests/unit/test_error_analyzer.py
