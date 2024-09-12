
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


<<<<<<< HEAD
    except ValueError as e:
        st.error(str(e))
=======
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
             background-color: #a19897; /* Set background color to white */
         }}
         </style>
         """,
         unsafe_allow_html=True,
     )

set_bg_hack('background6.jpg')

# st.markdown("""
#     # TARGET OF INTREST
# """)

#user input target

title = st.text_input('# TARGET OF INTREST')
csv_file = 'hdac_targets.csv'
options_df = pd.read_csv(csv_file)
#Returning the df with targets and rest
options_df = options_df[['organism', 'pref_name', 'score', 'target_chembl_id', 'target_type']]
#Hide index
# hdf = options_df.assign(hack='').set_index('hack')

# #filter daframe by index
options_list = options_df.index.tolist()


if st.button('Enter'):
#Read csv file of target
     st.write("DataFrame:")
     st.dataframe(options_df)




selected_rows = st.multiselect(
    'Select one or more rows:',
    options_df.index,
    format_func=lambda x: options_df.loc[x, 'pref_name']
)

if st.button('Submit'):
    filtered_df = options_df.loc[selected_rows]
    st.write("Filtered DataFrame:")
    st.dataframe(filtered_df)



#Read csv file of target
csv_file = 'hdac.csv'
options_df2 = pd.read_csv(csv_file)

if st.button('Submit Filtered DataFrame'):
    st.write('Bioactivity Molecule DataFrame:')
    st.dataframe(options_df2, height=500)  # Set an appropriate height for the scroll bar option


morgan_fp = st.checkbox('Morgan Fingerprint Descriptors')
padel_desc = st.checkbox('Get Padel Descriptors')

if morgan_fp or padel_desc:
    if st.button('Submit Descriptors'):
        if morgan_fp:
            st.write('Morgan Fingerprint Descriptors selected.')
        elif padel_desc:
            st.write('Padel Descriptors selected.')
        else:
            padel_desc and morgan_fp
            st.write('Hold On This Might Take A Sec')


# Placeholder for descriptor calculations
if morgan_fp:
    st.write('Calculating Morgan Fingerprint Descriptors...')
    # Add Morgan Fingerprint Descriptor calculation here
if padel_desc:
    st.write('Calculating Padel Descriptors...')
    # Add Padel Descriptor calculation here

#Read csv file of target
csv_file = 'df_merged_final.csv'
options_df2 = pd.read_csv(csv_file)
# dai = options_df2.isnull().sum()
# dai[dai>0]
options_df3 = options_df2.drop(columns=['BCUT2D_MWHI', 'BCUT2D_CHGHI', 'BCUT2D_CHGLO', 'BCUT2D_LOGPHI', 'BCUT2D_LOGPLOW', 'BCUT2D_MRHI', 'BCUT2D_MWLOW', 'BCUT2D_MRLOW'])
options_df3.dropna(inplace=True)
options_df3 = options_df3[~options_df3.isin([np.inf, -np.inf]).any(axis=1)]
if st.button('Molecular Fingerprints DataFrame'):
    st.write('Molecular SMILES Features DataFrame:')
    st.dataframe(options_df2, height=500)  # Set an appropriate height for the scroll bar option

#Read csv file of target
csv_file_tsne = 'X_TSNE.csv'
X_TSNE = pd.read_csv(csv_file_tsne)

if st.button('Submit for ML'):

    # Assuming descriptors are stored in df_descriptors

    # Placeholder for descriptor DataFrame (Replace with actual descriptors data)
    df_descriptors = options_df3.copy()  # Replace with actual descriptors data

    if not df_descriptors.empty:
        st.write('Calculating KMeans and t-SNE...')

        # KMeans Clustering
        kmeans = KMeans(n_clusters=15)  # Number of clusters can be parameterized
        kmeans_labels = kmeans.fit_predict(X_TSNE)
        # t-SNE dimensionality reduction
        tsne = TSNE(n_components=2)  # Number of dimensions can be parameterized
        tsne_results = tsne.fit_transform(X_TSNE)

        # Create a DataFrame with t-SNE results and KMeans labels
        tsne_df = pd.DataFrame(tsne_results, columns=['Dim1', 'Dim2'])
        tsne_df['KMeans_Labels'] = kmeans_labels
        tsne_df['KMeans_Labels'] =tsne_df['KMeans_Labels'].astype("category")
        # Calculate Silhouette Score
        silhouette_score = silhouette_score(df_descriptors, kmeans_labels)
        st.write(f'Silhouette Score: {silhouette_score}')

        # Plotly 3D Scatter Plot
        fig = px.scatter(tsne_df, x='Dim1', y='Dim2', color='KMeans_Labels')
        st.plotly_chart(fig)
    else:
        st.write('Descriptor DataFrame is empty. Please check the data.')












# # Save to CSV
# df_merged_final.to_csv('df_merged_final.csv', index=False)
# st.write('df_merged_final.csv with descriptors is ready for download')

# # Provide download link
# with open('df_merged_final.csv', 'rb') as f:
#     csv_data = f.read()
#     b64 = base64.b64encode(csv_data).decode()
#     href = f'<a href="data:file/csv;base64,{b64}" download="df_merged_final.csv">Download df_merged_final.csv</a>'
#     st.markdown(href, unsafe_allow_html=True)



# f 'df_merged_final.csv' in locals():
#     if st.button('Return Descriptors CSV'):
#         st.write('Returning Descriptors DataFrame as CSV:')
#         st.dataframe(filtered_df, height=500)
#         csv = filtered_df.to_csv(index=False)
#         b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
#         href = f'<a href="data:file/csv;base64,{b64}" download="descriptors.csv">Download CSV File</a>'
#         st.markdown(href, unsafe_allow_html=True)
>>>>>>> 5fbeea226168a9ac98a9d2da95ca4978eca177ab
