# Multi-Omics Statistical Integration with UMAP App- VIZZHY

## Overview

This application integrates multi-omics data (Genomics, Transcriptomics, and Proteomics) and applies clustering and dimensionality reduction (via UMAP) to uncover meaningful insights across different biological layers. The app also performs enrichment analysis, visualizes clusters, and generates interactive networks to explore gene, protein, pathway, metabolite, and disease associations.

## Features

- **Multi-Omics Data Upload:** Users can upload Genomics, Transcriptomics, and Proteomics CSV files for analysis.
- **Data Filtering:** Customizable filters allow users to filter genomics (CADD), transcriptomics (logFC and p-value), and proteomics (intensity) data.
- **UMAP Clustering:** The application uses UMAP for dimensionality reduction and KMeans clustering to group similar omics data points.
- **Enrichment Analysis:** Over-representation analysis (ORA) is performed for each cluster across multiple biological libraries (Pathways, Metabolites, Diseases).
- **UMAP and Network Visualization:** The app provides interactive UMAP scatter plots and cluster-based interaction networks to explore relationships between genes, proteins, and pathways.
- **Cluster Enrichment Summary:** Displays a table summarizing the top enrichment terms across clusters.
- **Network Interaction:** Visualizes interactions between genes, proteins, and enrichment terms in an interactive network.

## Usage (link: https://umapmointegration.streamlit.app/)

1. **Upload Your Data:**
   - Upload your **Genomics CSV** (with columns: `Gene`, `CADD`).
   - Upload your **Transcriptomics CSV** (with columns: `Gene`, `logFC`, `p_value`).
   - Upload your **Proteomics CSV** (with columns: `Gene`, `Protein`, `Intensity`).

2. **Set Filters:**
   - Set the minimum CADD score for genomics filtering.
   - Set the minimum absolute logFC for transcriptomics filtering.
   - Set the maximum p-value for transcriptomics filtering.
   - Set the minimum intensity for proteomics filtering.
   - Choose the number of clusters to generate in the clustering analysis.

3. **View Filtered Data:**
   The app will show a preview of the filtered data for each omics layer (Genomics, Transcriptomics, Proteomics).

4. **UMAP Clustering:**
   The data is clustered using KMeans after dimensionality reduction via UMAP. You can visualize the UMAP plot and interact with the clusters.

5. **Enrichment Analysis:**
   - Enrichment analysis will be performed for each cluster based on selected biological libraries (Pathway, Metabolite, Disease).
   - The results will be displayed in a table showing the top enrichment terms per cluster.

6. **Network Visualization:**
   An interactive network of gene-protein-term interactions for each cluster is displayed, with nodes representing genes, proteins, pathways, metabolites, and diseases. 

7. **Cluster Summary:**
   - The app provides a summary table showing gene interactions across clusters, including associated proteins, pathways, metabolites, and diseases.

## Key Features and Visualizations

- **UMAP Visualization with Clusters:**
   An interactive scatter plot showing the UMAP representation of the integrated data with clusters highlighted.
   
- **Cluster Boundaries:**
   Circular boundaries are drawn around each cluster on the UMAP plot to visualize the extent of the clusters.

- **Cluster Enrichment Bar Plot:**
   A bar chart displaying the frequency of enrichment terms (Pathway, Metabolite, Disease) across all clusters.

- **Interactive Network:**
   A dynamic network showing gene-protein and gene-term relationships. You can hover over nodes to see more details and explore the associations visually.

## Code Description

### 1. Data Upload and Filtering
- Users upload their omics data (Genomics, Transcriptomics, Proteomics) through file uploaders.
- Filters are applied based on user input values for CADD score (Genomics), logFC (Transcriptomics), p-value (Transcriptomics), and Intensity (Proteomics).

### 2. UMAP and Clustering
- Data is scaled using `StandardScaler` and reduced to 2D using `UMAP`.
- KMeans clustering is applied on the UMAP-reduced data to identify clusters.

### 3. Enrichment Analysis
- The app performs enrichment analysis for each cluster using the `gseapy.enrichr` method for various libraries (Reactome Pathway, HMDB Metabolites, DisGeNET Diseases).
- The top terms for each enrichment type are displayed.

### 4. Network Visualization
- An interactive network is generated using `pyvis.network`, showing relationships between genes, proteins, and enrichment terms.
- The network allows users to interact with nodes and explore the connections dynamically.

### 5. Final Outputs
- **UMAP Visualization**: A scatter plot showing the distribution of data points with cluster colors.
- **Enrichment Summary**: A table summarizing enrichment results for each cluster.
- **Network**: An interactive graph showing interactions between genes, proteins, pathways, metabolites, and diseases.
- **Cluster Interaction Table**: A tabular view of gene-protein-enrichment relationships across clusters.

## Notes

- The app utilizes the UniProt REST API to map protein IDs to gene names in the proteomics dataset.
- The app supports interactive elements like network visualization and UMAP scatter plots to help users explore the multi-omics integration results.
