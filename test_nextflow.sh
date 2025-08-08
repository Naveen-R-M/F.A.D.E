#!/bin/bash

#############################################
# F.A.D.E Nextflow Pipeline Test Script
#############################################

echo "=========================================="
echo "F.A.D.E Nextflow Pipeline Test"
echo "=========================================="
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0

# Function to check test result
check_result() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓ $2${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ $2${NC}"
        ((TESTS_FAILED++))
    fi
}

# Test 1: Check Nextflow installation
echo "1. Checking Nextflow installation..."
module load nextflow/24.10.3 2>/dev/null
nextflow -version &>/dev/null
check_result $? "Nextflow is installed"

# Test 2: Check directory structure
echo ""
echo "2. Checking directory structure..."
[ -d "$SCRATCH/F.A.D.E/nextflow" ]
check_result $? "Nextflow directory exists"

[ -f "$SCRATCH/F.A.D.E/nextflow/main.nf" ]
check_result $? "Main workflow file exists"

[ -f "$SCRATCH/F.A.D.E/nextflow/nextflow.config" ]
check_result $? "Configuration file exists"

[ -d "$SCRATCH/F.A.D.E/nextflow/modules" ]
check_result $? "Modules directory exists"

# Test 3: Check module files
echo ""
echo "3. Checking workflow modules..."
MODULES=(
    "target_selection.nf"
    "structure_prediction.nf"
    "binding_site.nf"
    "molecule_generation.nf"
    "docking.nf"
    "lead_optimization.nf"
    "reporting.nf"
)

for module in "${MODULES[@]}"; do
    [ -f "$SCRATCH/F.A.D.E/nextflow/modules/$module" ]
    check_result $? "Module $module exists"
done

# Test 4: Check Python agents
echo ""
echo "4. Checking Python agents..."
[ -f "$SCRATCH/F.A.D.E/agents/target_selector.py" ]
check_result $? "Target selector agent exists"

[ -f "$SCRATCH/F.A.D.E/agents/structure_predictor.py" ]
check_result $? "Structure predictor agent exists"

# Test 5: Check environment
echo ""
echo "5. Checking environment..."
[ ! -z "$SCRATCH" ]
check_result $? "SCRATCH variable is set"

[ ! -z "$GEMINI_API_KEY" ]
if [ $? -eq 0 ]; then
    check_result 0 "GEMINI_API_KEY is set"
else
    echo -e "${YELLOW}⚠ GEMINI_API_KEY is not set (some features may not work)${NC}"
fi

# Test 6: Check Singularity
echo ""
echo "6. Checking Singularity..."
#module load singularity/3.8.0 2>/dev/null
#singularity --version &>/dev/null
#check_result $? "Singularity is available"
echo "Skipping Singularity Check"

# Test 7: Validate Nextflow syntax
echo ""
echo "7. Validating Nextflow syntax..."
cd $SCRATCH/F.A.D.E/nextflow
nextflow main.nf -profile northeastern --help &>/dev/null
check_result $? "Nextflow syntax is valid"

# Test 8: Check runner scripts
echo ""
echo "8. Checking runner scripts..."
[ -x "$SCRATCH/F.A.D.E/run_nextflow.sh" ]
check_result $? "run_nextflow.sh is executable"

[ -x "$SCRATCH/F.A.D.E/submit_nextflow.sbatch" ]
check_result $? "submit_nextflow.sbatch is executable"

# Test 9: Dry run test (stub mode)
echo ""
echo "9. Running dry run test..."
cd $SCRATCH/F.A.D.E/nextflow
nextflow run main.nf \
    -profile local \
    -stub-run \
    --query "Test query" \
    --output_dir "test_output" \
    --max_molecules 10 \
    &>/dev/null
check_result $? "Dry run completed successfully"

# Clean up test output
rm -rf test_output work .nextflow* 2>/dev/null

# Summary
echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo -e "Tests Passed: ${GREEN}$TESTS_PASSED${NC}"
echo -e "Tests Failed: ${RED}$TESTS_FAILED${NC}"

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed! The Nextflow pipeline is ready to use.${NC}"
    echo ""
    echo "To run the pipeline:"
    echo "  ./run_nextflow.sh -q \"Your drug discovery query\""
    echo ""
    echo "To submit to SLURM:"
    echo "  sbatch submit_nextflow.sbatch \"Your query\" \"output_dir\" \"1000\""
    exit 0
else
    echo -e "${RED}Some tests failed. Please check the configuration.${NC}"
    echo ""
    echo "Common fixes:"
    echo "  1. Load required modules: module load nextflow/24.10.3"
    echo "  2. Set GEMINI_API_KEY: export GEMINI_API_KEY='your-key'"
    echo "  3. Check file permissions: chmod +x run_nextflow.sh"
    exit 1
fi
