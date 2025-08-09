#!/bin/bash

# Setup script to resolve Java conflicts for Nextflow

echo "🔧 Setting up Nextflow environment..."

# Step 1: Clear any existing Java environment
echo "Step 1: Clearing existing Java environment..."
unset JAVA_HOME
unset JAVA_PATH
unset JDK_HOME

# Step 2: Purge all modules to start clean
echo "Step 2: Purging all modules..."
module purge

# Step 3: Load required modules in correct order
echo "Step 3: Loading modules..."
module load miniconda3/24.11.1
module load nextflow/24.10.3

# Step 4: Verify Java setup
echo "Step 4: Verifying Java setup..."
echo "JAVA_HOME: $JAVA_HOME"
echo "Java version:"
java -version

# Step 5: Test Nextflow
echo "Step 5: Testing Nextflow..."
nextflow -version

echo ""
echo "✅ Environment setup complete!"
echo ""
echo "You can now run:"
echo "  cd $SCRATCH/F.A.D.E"
echo "  ./run_nextflow_integrated.sh --dry-run 'Test query'"
