import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle

import streamlit as st
import pandas as pd


csv_file = 'hdac_targets.csv'
options_df = pd.read_csv(csv_file)
#Add index
options_df = options_df[['organism', 'pref_name', 'score', 'target_chembl_id', 'target_type']]
options_df = options_df.reset_index()


#filter daframe by index
options_list = options_df['index'].tolist()
st.write("DataFrame:")
st.dataframe(options_df)

selection = st.sidebar.selectbox('Select an option:', options_list)
if selection:
    # do something with the selected option
    print(f"User selected: {selection}")
