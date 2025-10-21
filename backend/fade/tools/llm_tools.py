"""
LLM integration tools for natural language processing in the pipeline.
"""

from typing import Dict, Any, List, Optional, Tuple
import json
import re
from langchain.chat_models import init_chat_model
from langchain.prompts import ChatPromptTemplate
from langchain.schema import BaseMessage, HumanMessage, SystemMessage
from langchain.output_parsers import PydanticOutputParser
from pydantic import BaseModel, Field

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.llm_tools")


class TargetExtractionOutput(BaseModel):
    """Structured output for target extraction from natural language."""
    protein_name: Optional[str] = Field(description="Common name of the protein target")
    gene_name: Optional[str] = Field(description="Gene symbol (e.g., KRAS, EGFR)")
    uniprot_id: Optional[str] = Field(description="UniProt accession if mentioned")
    organism: str = Field(default="Homo sapiens", description="Target organism")
    mutations: List[str] = Field(default_factory=list, description="Specific mutations mentioned (e.g., G12C, T790M)")
    binding_site: Optional[str] = Field(description="Specific binding site or pocket mentioned")
    disease_context: Optional[str] = Field(description="Disease or condition mentioned")
    desired_properties: List[str] = Field(default_factory=list, description="Desired drug properties mentioned")
    constraints: List[str] = Field(default_factory=list, description="Any constraints mentioned (e.g., BBB penetrant, oral bioavailability)")


class LLMTools:
    """Tools for using LLMs in the drug discovery pipeline."""
    
    def __init__(self):
        """Initialize LLM with configuration."""
        llm_config = config.get_llm_config()
        self.llm = init_chat_model(
            model=llm_config["model"],
            api_key=llm_config["api_key"],
            temperature=llm_config.get("temperature", 0.7),
            max_tokens=llm_config.get("max_tokens", 4000)
        )
        logger.info(f"Initialized LLM: {llm_config['model']}")
    
    def extract_target_from_query(self, query: str) -> TargetExtractionOutput:
        """
        Extract protein target information from natural language query.
        
        Args:
            query: Natural language drug discovery query
            
        Returns:
            Structured target information
        """
        parser = PydanticOutputParser(pydantic_object=TargetExtractionOutput)
        
        prompt_template = ChatPromptTemplate.from_messages([
            SystemMessage(content="""You are a expert in drug discovery and molecular biology.
Your task is to extract protein target information from natural language queries.

Extract the following information:
- Protein name (common name like "KRAS protein")
- Gene name (official gene symbol like "KRAS")
- UniProt ID (if explicitly mentioned, like "P01116")
- Organism (default to "Homo sapiens" for human)
- Specific mutations (like "G12C", "T790M", "V600E")
- Binding site or pocket (if mentioned, like "GTP binding site", "kinase domain")
- Disease context (cancer type, neurological condition, etc.)
- Desired drug properties (BBB penetrant, oral bioavailability, low toxicity)
- Any constraints or requirements

Be precise with mutation notation - use standard format like G12C (amino acid position amino acid).
For gene names, use standard HUGO gene nomenclature.

{format_instructions}"""),
            HumanMessage(content=f"Extract target information from this query:\n\n{query}")
        ])
        
        messages = prompt_template.format_messages(
            format_instructions=parser.get_format_instructions()
        )
        
        try:
            response = self.llm.invoke(messages)
            
            # Try to parse the response
            try:
                result = parser.parse(response.content)
                logger.info(f"Successfully extracted target: {result.gene_name or result.protein_name}")
                return result
            except Exception as e:
                # Fallback: try to extract JSON from the response
                logger.warning(f"Failed to parse with Pydantic, trying JSON extraction: {e}")
                json_match = re.search(r'\{.*\}', response.content, re.DOTALL)
                if json_match:
                    json_str = json_match.group()
                    data = json.loads(json_str)
                    result = TargetExtractionOutput(**data)
                    return result
                else:
                    # Final fallback with defaults
                    logger.error("Could not parse LLM response, using minimal extraction")
                    return TargetExtractionOutput(
                        protein_name=self._extract_protein_name(query),
                        gene_name=self._extract_gene_name(query)
                    )
                    
        except Exception as e:
            logger.error(f"LLM extraction failed: {e}")
            # Return minimal extraction based on patterns
            return TargetExtractionOutput(
                protein_name=self._extract_protein_name(query),
                gene_name=self._extract_gene_name(query)
            )
    
    def _extract_protein_name(self, text: str) -> Optional[str]:
        """Simple pattern-based protein name extraction."""
        # Common patterns for protein names
        patterns = [
            r'(KRAS|EGFR|BRAF|ALK|ROS1|MET|HER2|PIK3CA|PTEN|TP53|CDK4|CDK6)',
            r'([A-Z][A-Z0-9]{2,})\s+(?:protein|kinase|receptor|enzyme)',
            r'(?:target|inhibit|block|modulate)\s+([A-Z][A-Z0-9]{2,})'
        ]
        
        for pattern in patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                return match.group(1)
        return None
    
    def _extract_gene_name(self, text: str) -> Optional[str]:
        """Simple pattern-based gene name extraction."""
        # Common oncogene and drug target patterns
        gene_pattern = r'\b([A-Z][A-Z0-9]{1,}(?:[A-Z0-9]|[-/]?[A-Z0-9]+)?)\b'
        
        # Known common drug targets
        known_targets = {
            'KRAS', 'EGFR', 'BRAF', 'ALK', 'ROS1', 'MET', 'HER2', 'ERBB2',
            'PIK3CA', 'PTEN', 'TP53', 'CDK4', 'CDK6', 'VEGFR', 'FGFR',
            'BCL2', 'MCL1', 'BTK', 'JAK1', 'JAK2', 'JAK3', 'TYK2',
            'FLT3', 'IDH1', 'IDH2', 'HDAC', 'PARP', 'PD1', 'PDL1',
            'CTLA4', 'TROP2', 'NECTIN4', 'FOLR1', 'CD19', 'CD20', 'CD30'
        }
        
        words = text.upper().split()
        for word in words:
            clean_word = re.sub(r'[^\w-]', '', word)
            if clean_word in known_targets:
                return clean_word
            
        # Try generic pattern
        match = re.search(gene_pattern, text)
        if match and len(match.group(1)) >= 3:
            return match.group(1)
            
        return None
    
    def clarify_ambiguous_target(self, 
                                query: str, 
                                candidates: List[Dict[str, Any]]) -> int:
        """
        Use LLM to select the best matching target from candidates.
        
        Args:
            query: Original user query
            candidates: List of potential protein targets
            
        Returns:
            Index of the best matching candidate
        """
        if not candidates:
            return -1
        if len(candidates) == 1:
            return 0
            
        # Format candidates for LLM
        candidate_text = "\n".join([
            f"{i+1}. {c.get('gene_name', 'Unknown')} - {c.get('protein_name', 'Unknown protein')} "
            f"({c.get('organism', 'Unknown organism')})"
            for i, c in enumerate(candidates)
        ])
        
        prompt = ChatPromptTemplate.from_messages([
            SystemMessage(content="""You are a drug discovery expert. 
Select the most appropriate protein target based on the user's query.
Consider the context, disease mentioned, and common drug targets.
Respond with just the number of the best match."""),
            HumanMessage(content=f"""Query: {query}
            
Available targets:
{candidate_text}

Which target best matches the query? Respond with just the number:""")
        ])
        
        try:
            response = self.llm.invoke(prompt.format_messages())
            
            # Extract number from response
            number_match = re.search(r'(\d+)', response.content)
            if number_match:
                index = int(number_match.group(1)) - 1
                if 0 <= index < len(candidates):
                    logger.info(f"LLM selected candidate {index + 1}: {candidates[index].get('gene_name')}")
                    return index
                    
        except Exception as e:
            logger.error(f"Failed to clarify target: {e}")
            
        # Default to first candidate
        return 0
    
    def generate_pocket_selection_rationale(self,
                                           pockets: List[Dict[str, Any]],
                                           selected_pocket: Dict[str, Any],
                                           target_info: Dict[str, Any]) -> str:
        """
        Generate a rationale for why a specific pocket was selected.
        
        Args:
            pockets: All identified pockets
            selected_pocket: The selected pocket
            target_info: Information about the target
            
        Returns:
            Natural language rationale
        """
        prompt = ChatPromptTemplate.from_messages([
            SystemMessage(content="""You are a structural biologist and drug designer.
Explain why a specific binding pocket was selected for drug design.
Consider druggability, size, hydrophobicity, and known binding sites.
Be concise but thorough."""),
            HumanMessage(content=f"""Target: {target_info.get('protein_name', 'Unknown')}
            
Selected pocket:
- Volume: {selected_pocket.get('volume', 0):.1f} ų
- Druggability score: {selected_pocket.get('druggability_score', 0):.2f}
- Surface area: {selected_pocket.get('surface_area', 0):.1f} ų
- Known binding site: {selected_pocket.get('is_known_site', False)}

Other pockets considered: {len(pockets) - 1}

Generate a brief rationale (2-3 sentences) for selecting this pocket:""")
        ])
        
        try:
            response = self.llm.invoke(prompt.format_messages())
            return response.content.strip()
        except Exception as e:
            logger.error(f"Failed to generate rationale: {e}")
            return f"Selected based on optimal druggability score ({selected_pocket.get('druggability_score', 0):.2f}) and suitable pocket volume."
