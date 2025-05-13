import streamlit as st
from PIL import Image
from pathlib import Path

st.set_page_config(page_title="Multi-Omics App Dashboard", layout="wide")
st.title("🧬 Multi-Omics Integration Dashboard")
st.markdown("Choose an analysis tool below to get started:")

# Create a row of 3 clickable images
cols = st.columns(3)

# Define image paths and links
apps = [
    {
        "title": "MOIntegration- Network Explorer",
        "image": "https://github.com/priyadarshinikp1/Multiomic_integration_DASHBOARD/blob/a503cd6379f14e32fe58739ce28cc74296d9719c/network.png?raw=true",  # Replace with actual image URL
        "page": "Multiomics_Integrator"
    },
    {
        "title": "MOIntegration- UMAP and clustering",
        "image": "https://github.com/priyadarshinikp1/Multiomic_integration_DASHBOARD/blob/a503cd6379f14e32fe58739ce28cc74296d9719c/umap.png?raw=true",
        "page": "MO_UMAP"
    },
    {
        "title": "App 3: Omics Pathway Mapper- Enrichmnet analysis",
        "image": "https://github.com/priyadarshinikp1/Multiomic_integration_DASHBOARD/blob/a503cd6379f14e32fe58739ce28cc74296d9719c/enrich.png?raw=true",
        "page": "Enrichment_DB"
    }
]

# Render each app block
for col, app in zip(cols, apps):
    with col:
        st.markdown(f"### {app['title']}", unsafe_allow_html=True)
        # Make image clickable
        st.markdown(
            f"""
            <a href="/{app['page']}" target="_self">
                <img src="{app['image']}" style="width:100%; border-radius: 10px; border: 2px solid #DDD;" />
            </a>
            """,
            unsafe_allow_html=True
        )
