import streamlit as st
source = st.selectbox('Source',[1,2])
impute = st.selectbox('Impute',[True, False])
input_data = st.text_input('Input')
