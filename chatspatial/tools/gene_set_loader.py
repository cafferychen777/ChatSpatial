"""
Gene Set Database Loader for Enrichment Analysis

This module provides functionality to load gene sets from various databases
for use with enrichment analysis tools.
"""

import logging
from typing import Dict, List, Optional, Union
import pandas as pd
import numpy as np
from pathlib import Path
import json
import requests
import gseapy as gp

logger = logging.getLogger(__name__)


class GeneSetLoader:
    """Load gene sets from various databases and sources."""
    
    def __init__(self, species: str = "human"):
        """
        Initialize gene set loader.
        
        Parameters
        ----------
        species : str
            Species for gene sets ('human' or 'mouse')
        """
        self.species = species.lower()
        self.organism = "Homo sapiens" if species == "human" else "Mus musculus"
        
    def load_msigdb(
        self,
        collection: str = "H",
        subcollection: Optional[str] = None,
        min_size: int = 10,
        max_size: int = 500
    ) -> Dict[str, List[str]]:
        """
        Load gene sets from MSigDB using gseapy.
        
        Parameters
        ----------
        collection : str
            MSigDB collection name:
            - H: hallmark gene sets
            - C1: positional gene sets
            - C2: curated gene sets (e.g., CGP, CP:KEGG, CP:REACTOME)
            - C3: motif gene sets
            - C4: computational gene sets
            - C5: GO gene sets (CC, BP, MF)
            - C6: oncogenic signatures
            - C7: immunologic signatures
            - C8: cell type signatures
        subcollection : Optional[str]
            Subcollection for specific databases (e.g., 'CP:KEGG', 'GO:BP')
        min_size : int
            Minimum gene set size
        max_size : int
            Maximum gene set size
            
        Returns
        -------
        Dict[str, List[str]]
            Dictionary of gene sets
        """
        try:
            # Get available gene sets
            gene_sets_dict = {}
            
            if collection == "H":
                # Hallmark gene sets
                gene_sets = gp.get_library_name(organism=self.organism)
                if 'MSigDB_Hallmark_2020' in gene_sets:
                    gene_sets_dict = gp.get_library('MSigDB_Hallmark_2020', organism=self.organism)
            
            elif collection == "C2" and subcollection == "CP:KEGG":
                # KEGG pathways
                if self.species == "human":
                    gene_sets_dict = gp.get_library('KEGG_2021_Human', organism=self.organism)
                else:
                    gene_sets_dict = gp.get_library('KEGG_2019_Mouse', organism=self.organism)
                    
            elif collection == "C2" and subcollection == "CP:REACTOME":
                # Reactome pathways
                gene_sets_dict = gp.get_library('Reactome_2022', organism=self.organism)
                
            elif collection == "C5":
                # GO gene sets
                if subcollection == "GO:BP" or subcollection is None:
                    gene_sets_dict.update(gp.get_library('GO_Biological_Process_2023', organism=self.organism))
                if subcollection == "GO:MF" or subcollection is None:
                    gene_sets_dict.update(gp.get_library('GO_Molecular_Function_2023', organism=self.organism))
                if subcollection == "GO:CC" or subcollection is None:
                    gene_sets_dict.update(gp.get_library('GO_Cellular_Component_2023', organism=self.organism))
                    
            elif collection == "C8":
                # Cell type signatures
                gene_sets_dict = gp.get_library('CellMarker_Augmented_2021', organism=self.organism)
            
            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets_dict.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes
                    
            logger.info(f"Loaded {len(filtered_sets)} gene sets from MSigDB {collection}")
            return filtered_sets
            
        except Exception as e:
            logger.error(f"Failed to load MSigDB gene sets: {e}")
            return {}
    
    def load_go_terms(
        self,
        aspect: str = "BP",
        min_size: int = 10,
        max_size: int = 500
    ) -> Dict[str, List[str]]:
        """
        Load GO terms using gseapy.
        
        Parameters
        ----------
        aspect : str
            GO aspect: 'BP' (biological process), 'MF' (molecular function), 'CC' (cellular component)
        min_size : int
            Minimum gene set size
        max_size : int
            Maximum gene set size
            
        Returns
        -------
        Dict[str, List[str]]
            Dictionary of GO gene sets
        """
        aspect_map = {
            "BP": "GO_Biological_Process_2023",
            "MF": "GO_Molecular_Function_2023",
            "CC": "GO_Cellular_Component_2023"
        }
        
        if aspect not in aspect_map:
            raise ValueError(f"Invalid GO aspect: {aspect}")
            
        try:
            gene_sets = gp.get_library(aspect_map[aspect], organism=self.organism)
            
            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes
                    
            logger.info(f"Loaded {len(filtered_sets)} GO {aspect} gene sets")
            return filtered_sets
            
        except Exception as e:
            logger.error(f"Failed to load GO gene sets: {e}")
            return {}
    
    def load_kegg_pathways(
        self,
        min_size: int = 10,
        max_size: int = 500
    ) -> Dict[str, List[str]]:
        """Load KEGG pathways."""
        try:
            if self.species == "human":
                gene_sets = gp.get_library('KEGG_2021_Human', organism=self.organism)
            else:
                gene_sets = gp.get_library('KEGG_2019_Mouse', organism=self.organism)
                
            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes
                    
            logger.info(f"Loaded {len(filtered_sets)} KEGG pathways")
            return filtered_sets
            
        except Exception as e:
            logger.error(f"Failed to load KEGG pathways: {e}")
            return {}
    
    def load_reactome_pathways(
        self,
        min_size: int = 10,
        max_size: int = 500
    ) -> Dict[str, List[str]]:
        """Load Reactome pathways."""
        try:
            gene_sets = gp.get_library('Reactome_2022', organism=self.organism)
            
            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes
                    
            logger.info(f"Loaded {len(filtered_sets)} Reactome pathways")
            return filtered_sets
            
        except Exception as e:
            logger.error(f"Failed to load Reactome pathways: {e}")
            return {}
    
    def load_cell_markers(
        self,
        min_size: int = 5,
        max_size: int = 200
    ) -> Dict[str, List[str]]:
        """Load cell type marker gene sets."""
        try:
            gene_sets = gp.get_library('CellMarker_Augmented_2021', organism=self.organism)
            
            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes
                    
            logger.info(f"Loaded {len(filtered_sets)} cell type marker sets")
            return filtered_sets
            
        except Exception as e:
            logger.error(f"Failed to load cell markers: {e}")
            return {}


async def load_gene_sets(
    database: str,
    species: str = "human",
    min_genes: int = 10,
    max_genes: int = 500,
    context = None
) -> Dict[str, List[str]]:
    """
    Load gene sets from specified database.
    
    Parameters
    ----------
    database : str
        Database name:
        - GO_Biological_Process, GO_Molecular_Function, GO_Cellular_Component
        - KEGG_Pathways
        - Reactome_Pathways
        - MSigDB_Hallmark
        - Cell_Type_Markers
    species : str
        Species ('human' or 'mouse')
    min_genes : int
        Minimum gene set size
    max_genes : int
        Maximum gene set size
    context : Optional
        MCP context for logging
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary of gene sets
    """
    loader = GeneSetLoader(species=species)
    
    database_map = {
        "GO_Biological_Process": lambda: loader.load_go_terms("BP", min_genes, max_genes),
        "GO_Molecular_Function": lambda: loader.load_go_terms("MF", min_genes, max_genes),
        "GO_Cellular_Component": lambda: loader.load_go_terms("CC", min_genes, max_genes),
        "KEGG_Pathways": lambda: loader.load_kegg_pathways(min_genes, max_genes),
        "Reactome_Pathways": lambda: loader.load_reactome_pathways(min_genes, max_genes),
        "MSigDB_Hallmark": lambda: loader.load_msigdb("H", None, min_genes, max_genes),
        "Cell_Type_Markers": lambda: loader.load_cell_markers(min_genes, max_genes),
    }
    
    if database not in database_map:
        raise ValueError(f"Unknown database: {database}. Available: {list(database_map.keys())}")
    
    if context:
        await context.info(f"Loading gene sets from {database} for {species}")
    
    gene_sets = database_map[database]()
    
    if context:
        await context.info(f"Loaded {len(gene_sets)} gene sets from {database}")
    
    return gene_sets