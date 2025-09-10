# F.A.D.E Testing

This directory contains scripts and files for testing the F.A.D.E (Fully Agentic Drug Engine) system.

## Available Test Files

1. **test_queries.txt**: A file containing test queries for batch processing.
2. **run_single_test.sh**: A bash script that runs a single test query using the conda environment.
3. **run_batch_tests.sh**: A bash script that runs several predefined test queries using the conda environment.

## How to Use

### 1. Running a Single Query

To run a single query, use the `run_single_test.sh` script:

```bash
cd $SCRATCH/F.A.D.E
./tests/run_single_test.sh "Find molecules targeting KRAS G12D with good BBB permeability"
```

This script will:
- Load the miniconda module
- Activate the fade conda environment
- Run the query
- Save results to a timestamped directory

### 2. Running Multiple Predefined Test Queries

To run the predefined test queries:

```bash
cd $SCRATCH/F.A.D.E
./tests/run_batch_tests.sh
```

This will run several test queries and save the results in a timestamped directory under `test_results/`.

### 3. Batch Processing with test_queries.txt

To process all queries in the `test_queries.txt` file, you need to activate the conda environment first:

```bash
cd $SCRATCH/F.A.D.E
module load miniconda3/24.11.1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $SCRATCH/conda-envs/fade
python main.py --batch-file tests/test_queries.txt --output-dir test_results/batch_test
```

## Adding New Test Queries

You can add new test queries to:

1. The `test_queries.txt` file (one query per line)
2. The `run_batch_tests.sh` script by adding new `run_test_query` function calls
