# # Generic packages
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import pickle
import random
import os

# Graphical packages
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns

# Dash related packages
import dash
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import dash_bio as dashbio

# Functions
from src.charts import *

current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)

app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], assets_folder="../assets")

def load_disease_data(disease):
    """Load disease specific data files"""
    if disease == 1:  # Kawasaki
        sigs = pd.read_csv('../data/preprocessed/sigs_kawasaki.csv', sep=",")
        grapher = pd.read_csv('../data/preprocessed/grapher_kawasaki.csv', sep=",")
        go = pd.read_csv('../data/preprocessed/go_kawasaki.csv', sep=",")
    elif disease == 2:  # Sars-cov2
        sigs = pd.read_csv('../data/preprocessed/sigs_sars-cov2.csv', sep=",")
        grapher = pd.read_csv('../data/preprocessed/grapher_sars-cov2.csv', sep=",")
        go = pd.read_csv('../data/preprocessed/go_sars-cov2.csv', sep=",")
    
    with open('../models/dds_object.pkl', 'rb') as f:
        dds = pickle.load(f)
    
    return sigs, grapher, go, dds

# Initial data load
sigs, grapher, go, dds = load_disease_data(1)  # Load Kawasaki data by default
sigs['Significance'] = np.abs(np.log10(sigs['pvalue']))
sigs['sorter'] = sigs['Significance']*sigs['log2FoldChange']

# Rest of your existing code remains the same until the update_graphs function

@app.callback(
    [Output('graph-1', 'figure'),
     Output('graph-2', 'figure'),
     Output('graph-3', 'figure'),
     Output('graph-4', 'figure')],
    [Input('disease_input', 'value'),
     Input('num_genes_input', 'value'),
     Input('patient_selection_input', 'value'),
     Input('log2fold_input', 'value'),
     Input('significance_input', 'value')]
)

def update_graphs(disease_input, num_genes_input, patient_selection_input, log2fold_input, significance_input):
    # Load disease-specific data
    sigs, grapher, go, dds = load_disease_data(disease_input)
    sigs['Significance'] = np.abs(np.log10(sigs['pvalue']))
    sigs['sorter'] = sigs['Significance']*sigs['log2FoldChange']
    
    # Data Loading
    if disease_input == 1:
        disease_title = 'Kawasaki Disease'
    elif disease_input == 2:
        disease_title = 'SARS-CoV-2'

    # Update labels based on thresholds
    conditions = [
        (sigs['Significance'] <= significance_input) & ((sigs['log2FoldChange'] <= log2fold_input[1]) | (sigs['log2FoldChange'] >= log2fold_input[0])),
        (sigs['Significance'] > significance_input) & (sigs['log2FoldChange'] > log2fold_input[1]),
        (sigs['Significance'] > significance_input) & (sigs['log2FoldChange'] < log2fold_input[0])
    ]
    sigs['label'] = np.select(conditions, ['Not significant', 'Upregulated', 'Downregulated'], default='Not significant')
    
    num_genes = int(num_genes_input/2)
    label_df = pd.concat(
            (sigs.sort_values('sorter')[-num_genes:],
            sigs.sort_values('sorter')[0 : num_genes]))
    list_annotation = list(label_df['Symbol'])

    # Volcano plot
    figure1 = volcano_plot(sigs, list_annotation, log2fold_input, significance_input, disease_title)
    
    # GO Enrichment plot
    figure2 = create_go_plot(go, disease_title)
    
    # PCA plot
    figure3 = PCA_plot(dds, grapher, patient_selection_input, disease_title)
    
    # Clustermap plot
    figure4 = clustermap_plot(grapher, list_annotation, patient_selection_input)
    
    return figure1, figure2, figure3, figure4

if __name__ == "__main__":
    app.run_server(debug=True)