# Omics Integration & Over-Representation Explorer App - (VIZZHY)

## Overview

The App is a user-friendly **Streamlit** web application designed to integrate and analyze **genomic**, **transcriptomic**, and **proteomic** data. It performs over-representation analysis (ORA) using several biological databases, helping researchers uncover biological insights and connections between genes, proteins, metabolites, enzymes, pathways, diseases, and transcription factors.

The app allows users to upload datasets, perform enrichment analysis, visualize the relationships between the data, and view a summary of associations.

## Features
- **Omics Data Integration**: Upload and combine genomic, transcriptomic, and proteomic datasets.
- **Over-Representation Analysis (ORA)**: Analyze gene sets against various biological databases such as KEGG, Reactome, and more.
- **Interactive Network Visualization**: Visualize enriched terms in an interactive network.
- **Summary Table**: Generate a summary table of associations between genes, proteins, metabolites, enzymes, pathways, diseases, and more.

## Databases Used for Over-Representation Analysis:
- **HMDB_Metabolites** (Metabolites)
- **DisGeNET** (Diseases)
- **OMIM_Disease** (Diseases)
- **KEGG_2021_Human** (Pathways)
- **Reactome_2016** (Pathways)
- **TRANSFAC_and_JASPAR_PWMs** (Transcription Factors)
- **GO_Biological_Process_2021** (Biological Processes)
- **KEA_2015** (Enzymes)

## How It Works  link: https://8dbmointegrator.streamlit.app/
1. **Upload Omics Data**: Upload your **genomic**, **transcriptomic**, and **proteomic** data in CSV format via the sidebar.
2. **Over-Representation Analysis (ORA)**: The app performs ORA for common genes across all three omic layers against the selected databases.
3. **Explore Results**: The app generates a list of enriched terms and visualizes these terms and their connections through a network.
4. **Interactive Network**: Interact with a network visualization to explore connections between genes, proteins, metabolites, and more.
5. **Summary Table**: A table is generated showing the associations of each gene with its corresponding biological elements.

