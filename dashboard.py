import streamlit as st
import subprocess
import os

st.set_page_config(page_title="Multi-Omics Integration Apps", layout="wide")

# Sidebar for navigation
st.sidebar.title("Multi-Omic Integration Dashboard")
app = st.sidebar.radio("Select an App", ["MO_UMAP", "Multiomics Integrator", "Enrichment DB"])

# Run the selected app
if app == "MO_UMAP":
    st.write("Launching MO_UMAP app...")
    # Use UTF-8 encoding to open the file
    with open("app1_MO_umap/Mo_umap.py", "r", encoding="utf-8") as file:
        exec(file.read())


elif app == "Multiomics Integrator":
    st.write("Launching Multiomics Integrator app...")
    # Use UTF-8 encoding to open the file
    with open("app2_Multiomics_Integrator_app/app.py", "r", encoding="utf-8") as file:
        exec(file.read())

elif app == "Enrichment DB":
    st.write("Launching Enrichment DB app...")
    with open("app3_int_enri_db/enri_db.py", "r", encoding="utf-8") as file:
        exec(file.read())
