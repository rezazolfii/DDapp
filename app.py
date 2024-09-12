
import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import streamlit as st
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, silhouette_score
import plotly.express as px
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from padelpy import padeldescriptor
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MiniBatchKMeans
import matplotlib.pyplot as plt
from generator import GeneratingFingerprints

st.markdown("""
    # Drug Design Using Agglomerative Clustering
    ## *Marc, Luciano and Reza*
""")
def set_bg_hack(main_bg):
    '''
    A function to set the background image of the main page
    '''
    # set bg name
    main_bg_ext = "png"
    st.markdown(
         f"""
         <style>
         .stApp {{
             background: url(data:image/{main_bg_ext};base64,{base64.b64encode(open(main_bg, "rb").read()).decode()});
             background-repeat: no-repeat;
             background-position: center top;
             background-size: cover;
             background-attachment: scroll;
             background-color: #A19897; /* Set background color to white */
         }}
         </style>
         """,
         unsafe_allow_html=True,
     )
set_bg_hack('background6.jpg')

# Streamlit Interface
st.title('Fingerprint Generation and Clustering')

# Input target name
target_name = st.text_input("Enter target name:")

# Initialize the GeneratingFingerprints class
generator = GeneratingFingerprints()

if target_name:
    try:
        # Fetch targets
        df_targets = generator.fetch_user_target(target_name)
        st.write("Targets found:", df_targets[["target_chembl_id", "pref_name", "target_type", "score"]])

        # Allow the user to select target indices
        selected_indices = st.multiselect("Select target indices", df_targets.index.tolist())

        if selected_indices:
            df_proteins = generator.extract_user_target(selected_indices)
            st.write("Proteins extracted:", df_proteins)

            # Generate Fingerprints
            if st.button("Generate Fingerprints"):
                df_fingerprints = generator.extracting_fingerprints()
                st.session_state['fingerprints'] = df_fingerprints  # Store fingerprints in session state
                st.write("Fingerprints generated!")

            # Check if fingerprints are already generated and stored in session state
            if 'fingerprints' in st.session_state:
                if st.button("Cluster and Predict"):
                    # Use fingerprints from session state
                    df_fingerprints = st.session_state['fingerprints']
                    generator.df_final = df_fingerprints  # Set it in the object

                    # Get the results from fit_and_predict
                    cluster_results, cluster_plot, num_clusters = generator.fit_and_predict()

                    # Save the results in session state
                    st.session_state['cluster_results'] = cluster_results
                    st.session_state['cluster_plot'] = cluster_plot
                    st.session_state['num_clusters'] = num_clusters

                    # Display the cluster results and the plot
                    st.write("Cluster results:", cluster_results)
                    st.pyplot(cluster_plot)  # Display the plot
                    st.write(f"Number of Clusters: {num_clusters}")


    except ValueError as e:
        st.error(str(e))
