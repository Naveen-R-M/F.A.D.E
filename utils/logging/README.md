# F.A.D.E Logging System

This directory contains the logging configuration for the F.A.D.E (Fully Agentic Drug Engine) system. The logging system is designed to provide structured and organized logs that are easy to find and analyze.

## Features

1. **Agent-Specific Logs**: Each agent in the system has its own log file, making it easy to track the activities of specific agents.
2. **Structured Directory Organization**: Logs are organized in a structured directory hierarchy:
   - `logs/`: Base directory for all logs
   - `logs/main/`: Main application logs
   - `logs/target_selector/`: Target Selector agent logs
   - `logs/structure_predictor/`: Structure Predictor agent logs
   - `logs/molecule_generator/`: Molecule Generator agent logs
   - `logs/evaluator/`: Evaluator agent logs
   - `logs/docking/`: Docking agent logs
   - `logs/refiner/`: Refiner agent logs
   - `logs/utils/`: Utility logs

3. **Run-Specific Logs**: Each run of the system gets a unique run ID (timestamp), and all logs for that run are grouped together.
4. **Console Output**: In addition to file logs, logs can also be output to the console for real-time monitoring.
5. **Output Directory Integration**: Logs can be stored within the output directory of a specific run, making it easy to find logs for a specific experiment.

## Usage

### Basic Usage

```python
from utils.logging import setup_logging, get_logger

# Set up logging
setup_logging(log_level="INFO", logs_dir="logs")

# Get a logger for a specific component
logger = get_logger("fade.main")
logger.info("This is a log message")

# Get a logger for an agent
agent_logger = get_logger("fade.agent.target_selector")
agent_logger.info("This is an agent log message")
```

### Output Directory Integration

```python
from utils.logging import setup_logging, get_logger

# Set up logging with output directory integration
setup_logging(log_level="INFO", logs_dir="results/my_experiment/logs")

# Get a logger
logger = get_logger("fade.main")
logger.info("Logs will be stored in the specified directory")
```

## Log Level

The logging system supports the following log levels:

- `DEBUG`: Detailed debug information
- `INFO`: General information about system operation
- `WARNING`: Warning messages about potential issues
- `ERROR`: Error messages about issues that affected system operation
- `CRITICAL`: Critical error messages about issues that halted system operation

## Log Format

Each log entry has the following format:

```
YYYY-MM-DD HH:MM:SS - logger_name - LOG_LEVEL - Log message
```

For example:

```
2025-07-24 18:30:45 - fade.agent.target_selector - INFO - Processing query: Find molecules targeting KRAS G12D
```

## Implementation Details

- The logging system is implemented in `log_config.py`
- `setup_logging()` configures the root logger and creates the necessary directory structure
- `get_logger()` returns a logger with the specified name, configured to log to the appropriate file
- Each agent uses the logger from the base agent class, which is configured to use the appropriate log file

## Finding Logs

Logs for a specific run are stored in directories named with the run ID (timestamp). The run ID is generated when `setup_logging()` is called, and all logs for that run are stored in directories with the run ID as part of the file name.

For example, logs for a run with ID `20250724_183045` would be stored in:

- `logs/main/fade_20250724_183045.log`
- `logs/target_selector/target_selector_20250724_183045.log`
- etc.

When using output directory integration, logs are stored in the specified output directory. For example, if the output directory is `results/my_experiment`, logs would be stored in `results/my_experiment/logs/`.
