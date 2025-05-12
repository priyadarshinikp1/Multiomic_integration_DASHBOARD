# Multi-Omic Integration Streamlit Dashboard

This is a unified Streamlit dashboard that combines three interactive data visualization and analysis apps:

- **MO_UMAP App** – For dimensionality reduction and visualization using UMAP.
- **Multiomics Integrator** – Integrates and visualizes multi-omics data.
- **Enrichment DB** – Provides enrichment analysis using biological databases.


**Getting Started**
To use the Integrator App, simply visit the application hosted on Streamlit:

https://3mointegratorapps.streamlit.app/

No installation is required. The app is ready to use directly from the browser.

**Prerequisites**
A modern web browser (e.g., Chrome, Firefox, Edge).

**How to Use the App**

Upload Data
- Upload your Genomics CSV (with columns: Gene, CADD).
- Upload your Transcriptomics CSV (with columns: Gene, logFC, p_value).
- Upload your Proteomics CSV (with columns: Gene, Protein, Intensity).

**Provided Sample Input:**
- cvd_genomics- WGS variant annotated csv- PRJNA264546
- cvd_transcriptomics- Bulk rna seq differential expression results- PRJNA394884
- cvd_proteomics- Protein intensity- 

**Note**
"Details for each app are provided separately in individual README files located within their respective app directories."
