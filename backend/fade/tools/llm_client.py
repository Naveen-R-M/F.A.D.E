"""
Custom LLM client for F.A.D.E - NO FALLBACKS.
Fails immediately if LLM is not available.
"""

import os
import requests
import json
from typing import Dict, Any, Optional
from langchain_openai import ChatOpenAI
from langchain_ollama import ChatOllama
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.llm_client")


def get_llm_client():
    """
    Get configured LLM client - NO FALLBACKS.
    Fails if LLM cannot be initialized.
    
    Raises:
        ValueError: If LLM cannot be initialized
    """
    if not config.AI_MODEL:
        raise ValueError("AI_MODEL not configured")
    
    model = config.AI_MODEL.lower()
    
    # Check if we're using Ollama (local)
    # if config.AI_API_BASE and "localhost" in config.AI_API_BASE:
    #     logger.info(f"Using Ollama with model: {config.AI_MODEL}")
        
    #     # Try to use ChatOllama from langchain_community
    #     from langchain_ollama import ChatOllama
    #     logger.info("Using ChatOllama from langchain_community")
    #     return ChatOllama(
    #         model=config.AI_MODEL,
    #         base_url=config.AI_API_BASE.replace("/v1", ""),  # ChatOllama doesn't need /v1
    #         temperature=0.7,
    #         num_predict=2000
    #     )
    
    # For cloud models, use the standard configuration
    logger.info(f"Using cloud model: {config.AI_MODEL}")
    llm_config = config.get_llm_config()  # This will raise error if invalid
    
    # Use ChatOpenAI for all cloud models
    # return ChatOpenAI(**llm_config)
    return ChatGoogleGenerativeAI(**llm_config)


def test_llm():
    """
    Test if LLM is working - NO FALLBACK.
    
    Raises:
        Exception: If LLM test fails
    """
    llm = get_llm_client()  # Will raise if fails
    response = llm.invoke("Say 'Hello, F.A.D.E is working!' in 5 words or less")
    logger.info(f"LLM test successful: {response.content}")
    return True
