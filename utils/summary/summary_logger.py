"""
High-level summary logger for F.A.D.E

This module provides a summary logger that generates concise, natural language summaries
of the user query and the outputs from each agent in the F.A.D.E pipeline.
"""

import os
import json
import time
import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Union

# Dictionary of agent summarizers
AGENT_SUMMARIZERS = {}

def summarize_agent_result(agent_type: str):
    """
    Decorator for registering agent result summarizers.
    
    Args:
        agent_type: Type of agent (e.g., "target_selector")
    """
    def decorator(func):
        AGENT_SUMMARIZERS[agent_type] = func
        return func
    return decorator


class SummaryLogger:
    """
    High-level summary logger for F.A.D.E.
    
    This class generates and manages concise, natural language summaries of the
    user query and the outputs from each agent in the F.A.D.E pipeline.
    """
    
    def __init__(self, output_dir: str):
        """
        Initialize the summary logger.
        
        Args:
            output_dir: Directory where results and logs are stored
        """
        self.output_dir = output_dir
        self.summary_file = os.path.join(output_dir, "summary.log")
        self.logger = logging.getLogger("fade.summary")
        
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize the summary log file
        with open(self.summary_file, "w") as f:
            f.write("# F.A.D.E Pipeline Summary\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    def log_query(self, query: str) -> None:
        """
        Log the user query to the summary file.
        
        Args:
            query: Natural language query from the user
        """
        with open(self.summary_file, "a") as f:
            f.write("## User Query\n\n")
            f.write(f"{query}\n\n")
        
        self.logger.info(f"Logged user query to summary file: {self.summary_file}")
    
    def log_agent_result(self, agent_type: str, result: Dict[str, Any]) -> None:
        """
        Log an agent's results to the summary file in natural language.
        
        Args:
            agent_type: Type of agent (e.g., "target_selector")
            result: Results dictionary from the agent
        """
        # Get the appropriate summarizer for this agent type
        summarizer = AGENT_SUMMARIZERS.get(agent_type)
        
        if summarizer:
            # Generate a natural language summary
            summary = summarizer(result)
            
            with open(self.summary_file, "a") as f:
                agent_name = agent_type.replace("_", " ").title()
                f.write(f"## {agent_name} Results\n\n")
                f.write(f"{summary}\n\n")
            
            self.logger.info(f"Logged {agent_type} results to summary file")
        else:
            self.logger.warning(f"No summarizer found for agent type: {agent_type}")
    
    def log_final_summary(self, final_results: Dict[str, Any]) -> None:
        """
        Log a final summary of the entire pipeline run.
        
        Args:
            final_results: Complete results from the pipeline
        """
        with open(self.summary_file, "a") as f:
            f.write("## Pipeline Summary\n\n")
            
            # Count the number of agents that produced results
            agents_completed = sum(1 for agent in ["target_selector", "structure_predictor", 
                                                 "molecule_generator", "evaluator", 
                                                 "docking", "refiner"] 
                                 if agent in final_results)
            
            total_agents = 6  # Total number of agents in the pipeline
            
            f.write(f"Pipeline progress: {agents_completed}/{total_agents} agents completed\n\n")
            
            if "molecules" in final_results and final_results["molecules"]:
                f.write(f"Generated {len(final_results.get('molecules', []))} candidate molecules\n\n")
            
            # Add timestamp
            f.write(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        self.logger.info(f"Logged final summary to summary file")


# Define summarizers for each agent type

@summarize_agent_result("target_selector")
def summarize_target_selector(result: Dict[str, Any]) -> str:
    """
    Generate a concise summary of the Target Selector agent results.
    
    Args:
        result: Results dictionary from the Target Selector agent
        
    Returns:
        Concise natural language summary
    """
    parsed_data = result.get("parsed_data", {})
    sequences = result.get("sequences", {})
    
    # Extract protein targets
    protein_targets = parsed_data.get("protein_targets", [])
    targets_str = ", ".join(target.get("name", "") for target in protein_targets if target.get("name"))
    
    # Extract mutations
    all_mutations = []
    for target in protein_targets:
        for mutation in target.get("mutations", []):
            original = mutation.get("original_residue", "")
            position = mutation.get("position", "")
            mutated = mutation.get("mutated_residue", "")
            if original and position and mutated:
                all_mutations.append(f"{original}{position}{mutated}")
    
    mutations_str = ", ".join(all_mutations)
    
    # Extract molecule properties
    molecule_props = parsed_data.get("molecule_properties", {})
    props = []
    if molecule_props.get("blood_brain_barrier_permeability"):
        props.append("BBB permeability")
    if molecule_props.get("lipinski_rule_of_five"):
        props.append("Lipinski's rule of five")
    if molecule_props.get("toxicity_requirements"):
        props.append("specific toxicity requirements")
    if molecule_props.get("solubility_requirements"):
        props.append("specific solubility requirements")
    
    props_str = ", ".join(props) if props else "no specific properties"
    
    # Build the summary
    summary_parts = []
    
    # Add target and mutation information
    if targets_str:
        if mutations_str:
            summary_parts.append(f"Identified {targets_str} with {mutations_str} mutation")
        else:
            summary_parts.append(f"Identified {targets_str}")
    
    # Add sequence information
    for target_name, seq_info in sequences.items():
        accession = seq_info.get("accession", "")
        length = seq_info.get("length", 0)
        organism = seq_info.get("organism", "")
        
        seq_summary = f"Retrieved {target_name} sequence"
        if accession:
            seq_summary += f" ({accession})"
        if length:
            seq_summary += f", {length} amino acids"
        if organism:
            seq_summary += f" from {organism}"
        
        summary_parts.append(seq_summary)
    
    # Add property requirements
    if props:
        summary_parts.append(f"Required properties: {props_str}")
    
    # Add configuration information
    config_files = result.get("config_files", {})
    if config_files:
        summary_parts.append("Generated configuration files for structure prediction")
    
    # Join all parts
    summary = " ".join(summary_parts) + "."
    
    return summary


@summarize_agent_result("structure_predictor")
def summarize_structure_predictor(result: Dict[str, Any]) -> str:
    """
    Generate a concise summary of the Structure Predictor agent results.
    
    Args:
        result: Results dictionary from the Structure Predictor agent
        
    Returns:
        Concise natural language summary
    """
    # Since this agent isn't implemented yet, provide a placeholder summary
    structures = result.get("structures", {})
    
    if not structures:
        return "No structures were generated."
    
    summary_parts = []
    
    for protein_name, structure_info in structures.items():
        confidence = structure_info.get("confidence", 0)
        binding_sites = structure_info.get("binding_sites", [])
        
        structure_summary = f"Generated 3D structure for {protein_name}"
        if confidence:
            confidence_pct = int(confidence * 100)
            structure_summary += f" with {confidence_pct}% confidence"
        
        if binding_sites:
            binding_sites_str = ", ".join(binding_sites)
            structure_summary += f". Identified binding sites: {binding_sites_str}"
        
        summary_parts.append(structure_summary)
    
    # Join all parts
    summary = " ".join(summary_parts) + "."
    
    return summary


@summarize_agent_result("molecule_generator")
def summarize_molecule_generator(result: Dict[str, Any]) -> str:
    """
    Generate a concise summary of the Molecule Generator agent results.
    
    Args:
        result: Results dictionary from the Molecule Generator agent
        
    Returns:
        Concise natural language summary
    """
    # Since this agent isn't implemented yet, provide a placeholder summary
    molecules = result.get("molecules", [])
    
    if not molecules:
        return "No molecules were generated."
    
    num_molecules = len(molecules)
    
    summary = f"Generated {num_molecules} candidate molecules"
    
    # Add property information if available
    properties_met = result.get("properties_met", [])
    if properties_met:
        props_str = ", ".join(properties_met)
        summary += f" with {props_str}"
    
    summary += "."
    
    return summary


@summarize_agent_result("evaluator")
def summarize_evaluator(result: Dict[str, Any]) -> str:
    """
    Generate a concise summary of the Evaluator agent results.
    
    Args:
        result: Results dictionary from the Evaluator agent
        
    Returns:
        Concise natural language summary
    """
    # Since this agent isn't implemented yet, provide a placeholder summary
    evaluated_molecules = result.get("evaluated_molecules", [])
    
    if not evaluated_molecules:
        return "No molecules were evaluated."
    
    num_molecules = len(evaluated_molecules)
    num_passed = sum(1 for m in evaluated_molecules if m.get("passed_filters", False))
    
    summary = f"Evaluated {num_molecules} molecules, {num_passed} passed all filters."
    
    return summary


@summarize_agent_result("docking")
def summarize_docking(result: Dict[str, Any]) -> str:
    """
    Generate a concise summary of the Docking agent results.
    
    Args:
        result: Results dictionary from the Docking agent
        
    Returns:
        Concise natural language summary
    """
    # Since this agent isn't implemented yet, provide a placeholder summary
    docked_molecules = result.get("docked_molecules", [])
    
    if not docked_molecules:
        return "No molecules were docked."
    
    num_molecules = len(docked_molecules)
    
    # Calculate average docking score
    scores = [m.get("docking_score", 0) for m in docked_molecules if "docking_score" in m]
    avg_score = sum(scores) / len(scores) if scores else 0
    
    summary = f"Docked {num_molecules} molecules with average score of {avg_score:.2f} kcal/mol."
    
    return summary


@summarize_agent_result("refiner")
def summarize_refiner(result: Dict[str, Any]) -> str:
    """
    Generate a concise summary of the Refiner agent results.
    
    Args:
        result: Results dictionary from the Refiner agent
        
    Returns:
        Concise natural language summary
    """
    # Since this agent isn't implemented yet, provide a placeholder summary
    refined_molecules = result.get("refined_molecules", [])
    
    if not refined_molecules:
        return "No molecules were refined."
    
    num_molecules = len(refined_molecules)
    num_improved = sum(1 for m in refined_molecules if m.get("improved", False))
    
    summary = f"Refined {num_molecules} molecules, {num_improved} showed improved properties or binding."
    
    return summary
