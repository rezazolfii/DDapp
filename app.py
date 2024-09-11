
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


#user input target
title = st.text_input('TARGET OF INTREST')

#Read csv file of target
csv_file = 'hdac_targets.csv'
options_df = pd.open_csv(csv_file)

#Returning the df with targets and rest
options_df = options_df[['organism', 'pref_name', 'score', 'target_chembl_id', 'target_type']]
#Hide index
# hdf = options_df.assign(hack='').set_index('hack')


# #filter daframe by index
options_list = options_df.index.tolist()
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
