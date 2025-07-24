#!/bin/bash

# Script to run unit tests for F.A.D.E

# Set up Python path
export PYTHONPATH="$PYTHONPATH:$(dirname $(dirname $(realpath $0)))"

# Run all unit tests
echo "Running all unit tests..."
python -m unittest discover -s tests/unit -p "test_*.py" -v

# Exit with the status of the tests
exit $?
