import os
import tempfile
import requests
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from pyvis.network import Network
from gseapy import enrichr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap

# -----------------------------
# App Configuration
# -----------------------------

st.image("https://raw.githubusercontent.com/priyadarshinikp1/Multiomics-Integrator-app/main/logo.png", width=200)
st.title("Multi-Omics Statistical Integration with UMAP Vizzhy App")

with st.sidebar:
    st.markdown("---")
    st.markdown("**Developed by: PRIYADARSHINI**")
    st.markdown("[LinkedIn](https://www.linkedin.com/in/priyadarshini24) | [GitHub](https://github.com/priyadarshinikp1)")

# -----------------------------
# Sidebar Configuration
# -----------------------------
st.sidebar.header("⚙️ Filter Settings")
cadd_thresh = float(st.sidebar.text_input("Min CADD Score (Genomics)", value="20"))
logfc_thresh = float(st.sidebar.text_input("Min |logFC| (Transcriptomics)", value="1"))
t_pval_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="1000"))
n_clusters = st.sidebar.slider("Number of Clusters", min_value=2, max_value=10, value=4)

# -----------------------------
# File Upload Section
# -----------------------------
st.header("📁 Upload Omics Data or Use Example Files")
use_example = st.checkbox("Use Example Data- Cardiomyopathy")
genomics = transcriptomics = proteomics = None

if use_example:
    input_dir = os.path.join(os.path.dirname(__file__), '..', 'Input')
    input_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'Input'))
    genomics_path = os.path.join(input_dir, "cvd_genomics.csv")
    transcriptomics_path = os.path.join(input_dir, "cvd_transcriptomics.csv")
    proteomics_path = os.path.join(input_dir, "cvd_proteomics.csv")
    try:
        genomics = open(genomics_path, "rb")
        transcriptomics = open(transcriptomics_path, "rb")
        proteomics = open(proteomics_path, "rb")
        st.success("Loaded example files")
    except FileNotFoundError:
        st.error("Example files not found. Please check the 'Input' folder and filenames.")
else:
    genomics = st.file_uploader("Upload Genomics CSV", type="csv")
    transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
    proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# -----------------------------
# Helper Functions
# -----------------------------
def extract_uniprot_ids(protein_series):
    ids = set()
    for entry in protein_series.dropna():
        for pid in str(entry).split(";"):
            if pid.strip():
                ids.add(pid.strip())
    return ids

@st.cache_data
def map_uniprot_to_gene(uniprot_ids):
    mapping = {}
    ids = list(uniprot_ids)
    for i in range(0, len(ids), 100):
        chunk = ids[i:i+100]
        query = " OR ".join([f"accession:{id_}" for id_ in chunk])
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
        try:
            r = requests.get(url)
            if r.status_code == 200:
                lines = r.text.strip().split('\n')[1:]
                for line in lines:
                    acc, genes = line.split('\t')
                    mapping[acc] = genes.split()[0] if genes else acc
        except:
            pass
    return mapping

# -----------------------------
# Main Processing
# -----------------------------

if genomics and transcriptomics and proteomics:
    with st.spinner("Processing and integrating omics data..."):
        try:
            # Load Data
            gdf = pd.read_csv(genomics)
            tdf = pd.read_csv(transcriptomics)
            pdf = pd.read_csv(proteomics)

            # Filter Data
            gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
            tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh)]

            if logfc_thresh > 0:
                tdf_filtered = tdf_filtered[tdf_filtered['logFC'] >= logfc_thresh]
            elif logfc_thresh < 0:
                tdf_filtered = tdf_filtered[tdf_filtered['logFC'] < logfc_thresh]

            pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

            # Preview
            st.subheader("🔍 Filtered Data Preview")
            top_n = st.slider("Show Top N Rows", 5, 50, 10)
            st.write("**Genomics**")
            st.dataframe(gdf_filtered.head(top_n))
            st.write("**Transcriptomics**")
            st.dataframe(tdf_filtered.head(top_n))
            st.write("**Proteomics**")
            st.dataframe(pdf_filtered.head(top_n))

            # 🔬 Common Genes Across All Omics
            common_genes_all_omics = set(gdf_filtered['Gene']) & set(tdf_filtered['Gene']) & set(pdf_filtered['Gene'])
            common_genes_all_omics_df = pd.DataFrame(list(common_genes_all_omics), columns=['Gene'])

            st.subheader("🔬 Common Genes Across All Omics Data")
            if not common_genes_all_omics:
                st.warning("No common genes found across all three omics datasets.")
            else:
                st.dataframe(common_genes_all_omics_df)

            # Integration
            union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])
            unique_uniprot_ids = extract_uniprot_ids(pdf_filtered['Protein'])
            uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)

            # Expand Proteins to Genes
            expanded_rows = []
            for _, row in pdf_filtered.iterrows():
                for pid in str(row['Protein']).split(';'):
                    pid = pid.strip()
                    gene = uniprot_gene_map.get(pid)
                    if gene:
                        expanded_rows.append({'Protein': pid, 'Gene': gene})
            expanded_protein_df = pd.DataFrame(expanded_rows)
            protein_gene_map = dict(zip(expanded_protein_df['Protein'], expanded_protein_df['Gene']))

            # Merge and cluster
            merged_df = pd.merge(gdf_filtered, tdf_filtered, on='Gene')
            merged_df = pd.merge(merged_df, pdf_filtered, on='Gene')
            scaler = StandardScaler()
            features = scaler.fit_transform(merged_df[['CADD', 'logFC', 'Intensity']])
            reducer = umap.UMAP(n_components=2, random_state=42)
            umap_coords = reducer.fit_transform(features)

            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            clusters = kmeans.fit_predict(umap_coords)
            merged_df['Cluster'] = clusters
            merged_df['UMAP1'], merged_df['UMAP2'] = umap_coords[:, 0], umap_coords[:, 1]

            # Enrichment per cluster
            enrichment_dict = {}
            all_enrichments = []
            libraries = {"Pathway": "Reactome_2016", "Metabolite": "HMDB_Metabolites", "Disease": "DisGeNET"}
            gene_cluster_map = {}

            for c in sorted(merged_df['Cluster'].unique()):
                cluster_genes = merged_df[merged_df['Cluster'] == c]['Gene'].dropna().unique().tolist()
                enrichments = {"Cluster": c, "Pathway": "", "Metabolite": "", "Disease": ""}

                for label, lib in libraries.items():
                    try:
                        enr = enrichr(gene_list=cluster_genes, gene_sets=lib, outdir=None)
                        if not enr.results.empty:
                            top_terms = "; ".join(enr.results.sort_values('P-value').head(3)['Term'])
                            enrichments[label] = top_terms
                            for term in enr.results['Term'].head(3):
                                all_enrichments.append((label, term))
                    except Exception as e:
                        st.warning(f"Enrichment failed for Cluster {c}, {label}: {e}")
                        enrichments[label] = ""

                enrichment_dict[c] = enrichments

                for gene in cluster_genes:
                    gene_cluster_map.setdefault((c, gene), {"Protein": [], "Pathway": [], "Metabolite": [], "Disease": []})
                    for prot, gname in protein_gene_map.items():
                        if gname == gene:
                            gene_cluster_map[(c, gene)]["Protein"].append(prot)
                    for label in ["Pathway", "Metabolite", "Disease"]:
                        for term in enrichments[label].split(';'):
                            if term.strip():
                                gene_cluster_map[(c, gene)][label].append(term.strip())

            # Summary Table
            st.subheader("📚 Enrichment Summary by Cluster")
            enrichment_summary_df = pd.DataFrame([
                {'Cluster': c, 'Metabolite': v['Metabolite'], 'Pathway': v['Pathway'], 'Disease': v['Disease']}
                for c, v in enrichment_dict.items()
            ])
            st.dataframe(enrichment_summary_df, use_container_width=True)

            # Bar Plot
            st.subheader("📊 Frequent Enrichment Terms Across Clusters")
            enrichment_df_all = pd.DataFrame(all_enrichments, columns=['Type', 'Term'])
            bar_data = enrichment_df_all.groupby(['Type', 'Term']).size().reset_index(name='Count')
            fig_bar = px.bar(bar_data, x='Term', y='Count', color='Type', title="Enrichment Term Frequency")
            st.plotly_chart(fig_bar)

            # UMAP
            st.subheader("🔬 UMAP Visualization with Clusters")
            fig_umap = px.scatter(
                merged_df, x='UMAP1', y='UMAP2', color=merged_df['Cluster'].astype(str),
                hover_name='Gene', hover_data=['CADD', 'logFC', 'Intensity'],
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            for cluster_id in merged_df['Cluster'].unique():
                cluster_data = merged_df[merged_df['Cluster'] == cluster_id]
                x_center, y_center = cluster_data[['UMAP1', 'UMAP2']].mean()
                radius = np.linalg.norm(cluster_data[['UMAP1', 'UMAP2']] - [x_center, y_center], axis=1).max()
                fig_umap.add_shape(type="circle", x0=x_center - radius, y0=y_center - radius,
                                   x1=x_center + radius, y1=y_center + radius,
                                   line=dict(color="black", width=2, dash="dash"),
                                   fillcolor="rgba(173,216,230,0.1)")
                fig_umap.add_annotation(x=x_center, y=y_center, text=f"Cluster {cluster_id}", showarrow=False, font=dict(size=14))
            st.plotly_chart(fig_umap, use_container_width=True)

            # Network
            st.subheader("🧠 Cluster-Based Interaction Network")
            net = Network(height='800px', width='100%', directed=False)
            net.force_atlas_2based()
            color_map = {"Gene": "gray", "Protein": "gold", "Pathway": "skyblue", "Metabolite": "lightgreen", "Disease": "lightcoral"}

            for (c, gene), info in gene_cluster_map.items():
                net.add_node(gene, label=gene, color=color_map['Gene'])
                for prot in info['Protein']:
                    net.add_node(prot, label=prot, color=color_map['Protein'])
                    net.add_edge(gene, prot, color='gray')
                for etype in ['Pathway', 'Metabolite', 'Disease']:
                    for term in info[etype]:
                        net.add_node(term, label=term, color=color_map[etype])
                        net.add_edge(gene, term, color='gray')

            with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
                net.save_graph(tmp_file.name)
                try:
                    html_content = open(tmp_file.name, 'r', encoding='utf-8').read()
                    st.components.v1.html(html_content, height=800)
                except Exception as e:
                    st.error(f"Failed to load network: {e}")

            st.subheader("🗂️ Network Legend")
            for k, v in color_map.items():
                st.markdown(f"<span style='display:inline-block; width:20px; height:20px; background-color:{v}; margin-right:10px; border-radius:50%;'></span> **{k}**", unsafe_allow_html=True)

            # Final Table
            st.subheader("📄 Interaction Table from Network")
            rows = []
            for (c, gene), info in gene_cluster_map.items():
                row = {
                    'Cluster': c,
                    'Gene': gene,
                    'Protein': "; ".join(info['Protein']),
                    'Metabolite': "; ".join(info['Metabolite']),
                    'Pathway': "; ".join(info['Pathway']),
                    'Disease': "; ".join(info['Disease'])
                }
                rows.append(row)
            st.dataframe(pd.DataFrame(rows))

        except Exception as e:
            st.error(f"An error occurred: {e}")
