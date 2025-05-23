# ===================================================================================
# Project: ChatSkLearn
# File: processor/data_ingestion.py
# Description: This file Loads the data from URLs and store it in a VectorStoreDB (Defaults to CHROMA)
# Author: LALAN KUMAR
# Created: [15-04-2025]
# Updated: [30-04-2025]
# LAST MODIFIED BY: LALAN KUMAR [https://github.com/kumar8074]
# Version: 1.0.0
# ===================================================================================
# CAUTION: DO NOT RUN this file until you want to make changes to the vectorStoreDB.
# RECOMMENDATION: It is recommended to run this file multiple times and persist the vectorDB using all the embedding providers.
#                 Currently it persists the vectorDB locally which is efiicient for ChatSklearn application types.
#                 However you can modify the script and persist the vectorDB externally to third party databases.

import os
import sys
import argparse

# Add root path
current_file_path = os.path.abspath(__file__)
project_root = os.path.abspath(os.path.join(current_file_path, "../.."))
if project_root not in sys.path:
    sys.path.append(project_root)

from langchain_community.document_loaders import WebBaseLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_chroma import Chroma

from logger import logging
from config import settings
from app.core.embeddings import get_embeddings

def load_urls(file_path: str) -> list[str]:
    """Load and clean URLs from a file."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return [line.strip() for line in lines if line.strip()]

def load_documents(urls: list[str]):
    """Load documents from web pages."""
    loader = WebBaseLoader(web_paths=urls)
    return loader.load()

def split_documents(documents):
    """Split large documents into chunks."""
    splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
    return splitter.split_documents(documents)

def persist_vector_db(docs, embeddings, persist_directory):
    """Save documents to Chroma vector store."""
    vectordb = Chroma.from_documents(
        documents=docs,
        embedding=embeddings,
        persist_directory=persist_directory
    )
    logging.info(f"VectorDB persisted at: {persist_directory}")

def main(embedding_provider=None):
    """Main function to ingest data, with optional provider specification."""
    # Use the specified provider or default from settings
    provider = embedding_provider or settings.EMBEDDING_PROVIDER
    logging.info(f"Using embedding provider: {provider}")
    
    urls_file_path = "processor/successful_urls.txt"
    logging.info(f"Loading URLs from {urls_file_path}")
    urls = load_urls(urls_file_path)

    logging.info(f"Loaded {len(urls)} URLs...")
    
    logging.info(f"Loading Documents from URLs....(This might take a while have a coffee break)")
    documents = load_documents(urls)
    logging.info(f"Loaded {len(documents)} documents...")
    
    logging.info(f"Splitting documents into chunks...")
    docs = split_documents(documents)
    logging.info(f"Split into {len(docs)} chunks...")

    # Get dynamic embedding model based on provider
    embeddings = get_embeddings(provider)
    logging.info(f"Using Embedding model: {embeddings}")

    # Set output directory based on provider name
    persist_directory = f"DATA/chroma_store_{provider.lower()}"
    os.makedirs(persist_directory, exist_ok=True)
    
    logging.info(f"Persisting vectorDB to {persist_directory}...(This might take a while)")
    persist_vector_db(docs, embeddings, persist_directory)

def ingest_all_providers():
    """Ingest data for all available providers."""
    providers = [
        settings.EMBEDDING_PROVIDER_GEMINI,
        settings.EMBEDDING_PROVIDER_OPENAI,
        settings.EMBEDDING_PROVIDER_COHERE
    ]
    
    for provider in providers:
        # Check if API key exists for this provider
        key_var_name = f"{provider.upper()}_API_KEY"
        key = getattr(settings, key_var_name, None)
        
        if key:
            logging.info(f"Processing data for {provider} provider")
            try:
                main(provider)
                logging.info(f"Successfully processed data for {provider}")
            except Exception as e:
                logging.error(f"Error processing data for {provider}: {str(e)}")
        else:
            logging.warning(f"Skipping {provider} - API key not found")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ingest data for vector stores")
    parser.add_argument("--provider", type=str, 
                        help="Specific embedding provider to use (gemini, openai, cohere)")
    parser.add_argument("--all", action="store_true", 
                        help="Process data for all providers with valid API keys")
    
    args = parser.parse_args()
    
    if args.all:
        ingest_all_providers()
    elif args.provider:
        # Convert provider name to the constant from settings
        provider_map = {
            "gemini": settings.EMBEDDING_PROVIDER_GEMINI,
            "openai": settings.EMBEDDING_PROVIDER_OPENAI,
            "cohere": settings.EMBEDDING_PROVIDER_COHERE
        }
        
        provider = provider_map.get(args.provider.lower())
        if provider:
            main(provider)
        else:
            logging.error(f"Unknown provider: {args.provider}")
    else:
        # Use default provider from settings
        main()
        
# To ingest data for all providers with valid API keys:
#python processor/data_ingestion.py --all

# Or for a specific provider:
#python processor/data_ingestion.py --provider gemini
#python processor/data_ingestion.py --provider openai
#python processor/data_ingestion.py --provider cohere
