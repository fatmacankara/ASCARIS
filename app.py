import streamlit as st

import pandas as pd 
from os import path
import sys
sys.path.append('code/')
#sys.path.append('ASCARIS/code/') 
#import pdb_featureVector
#import alphafold_featureVector
#import argparse
st.write('Hello')
st.write('ljl')
st.write('aaaaaaa')
source = st.selectbox('Source',[1,2])
impute = st.selectbox('Impute',[True, False])
input_data = st.text_input('Input')
