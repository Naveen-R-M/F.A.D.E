# F.A.D.E Summary Logger

This module provides a high-level summary logger for F.A.D.E that generates concise, natural language summaries of the user query and the outputs from each agent in the pipeline.

## Features

- **Concise Summaries**: Generates brief, natural language summaries of agent outputs
- **Human-Readable**: Designed to be easily understood by non-technical users
- **Pipeline Overview**: Provides a complete overview of the drug discovery pipeline
- **Agent-Specific Summaries**: Each agent has its own specialized summarizer

## How It Works

The summary logger uses specialized summarizers for each agent type to convert the technical output into natural language summaries. These summaries are then compiled into a single, easily readable log file.

## Summary Log Format

The summary log file has the following format:

```
# F.A.D.E Pipeline Summary

Date: YYYY-MM-DD HH:MM:SS

## User Query

The original user query goes here.

## Target Selector Results

Concise summary of target selector results.

## Structure Predictor Results

Concise summary of structure predictor results.

...

## Pipeline Summary

Pipeline progress: X/Y agents completed

Generated Z candidate molecules

Completed at: YYYY-MM-DD HH:MM:SS
```

## Example Summaries

### Target Selector Agent

```
Identified KRAS with G12D mutation. Retrieved KRAS sequence (A0A250XVZ2), 188 amino acids from Castor canadensis. Required properties: BBB permeability. Generated configuration files for structure prediction.
```

### Structure Predictor Agent

```
Generated 3D structure for KRAS with 85% confidence. Identified binding sites: GTP binding pocket, Switch II region.
```

## Implementation Details

The summary logger uses a decorator-based registration system for agent summarizers. Each agent type has a dedicated summarizer function that knows how to extract the relevant information from the agent's results and present it in a concise, natural language format.

The `SummaryLogger` class manages the creation and updating of the summary log file. It is initialized with the output directory and provides methods for logging the user query, agent results, and a final summary.

## Usage

```python
from utils.summary import SummaryLogger

# Initialize the summary logger
summary_logger = SummaryLogger(output_dir)

# Log the user query
summary_logger.log_query("Find molecules targeting KRAS G12D with good BBB permeability")

# Log agent results
summary_logger.log_agent_result("target_selector", target_selector_results)
summary_logger.log_agent_result("structure_predictor", structure_predictor_results)
# ... and so on for each agent

# Log final summary
summary_logger.log_final_summary(all_results)
```

## Extending

To add a summarizer for a new agent type, use the `@summarize_agent_result` decorator:

```python
from utils.summary import summarize_agent_result

@summarize_agent_result("new_agent_type")
def summarize_new_agent(result):
    # Extract relevant information from the result
    # Generate a concise summary
    return summary
```
