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

# Biocomputing related packages
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
from sanbomics.plots import volcano
import scanpy as sc # PCA library

# Functions
from charts import *

current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)

app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], assets_folder="../assets")

sigs = pd.read_csv('../data/preprocessed/sigs.csv', sep=",")
grapher = pd.read_csv('../data/preprocessed/grapher.csv', sep=",")
with open('../models/dds_object.pkl', 'rb') as f:
    dds = pickle.load(f)
sigs['Significance'] = np.abs(np.log10(sigs['pvalue']))
sigs['sorter'] = sigs['Significance']*sigs['log2FoldChange']

patient_list = ['KD1','KD2','KD3','KD4','KD5','KD6','KD7','KD8','KD9','KD10','KD11','KD12','KD13','KD14','KD15','KD16','KD17','KD18','KD19','KD20',
                'control1','control2','control3','control4','control5','control6','control7','control8','control9','control10','control11','control12','control13'
                ,'control14','control15','control16']
patient_list_init = ['KD1','KD2','KD3','KD4','KD5','KD6','KD7','KD8','KD9','KD10','KD11','KD12','KD13','KD14','KD15','KD16','KD17','KD18','KD19','KD20',
                'control1','control2','control3','control4','control5','control6','control7','control8','control9','control10']

# Initializing 
color_map = {
    'Not significant': '#D3D3D3',  # tomato color
    'Upregulated': '#1E90FF',  # tomato color
    'Downregulated': '#FF6347',  # dodger blue color
}

fig = go.Figure(
    go.Scattergl(
        x = np.random.randn(1000),
        y = np.random.randn(1000),
        mode='markers',
        marker=dict(color=random.sample(['#ecf0f1']*500 + ["#3498db"]*500, 1000), line_width=1)
    )
)
fig.update_layout(plot_bgcolor='#FFFFFF', #width=500, height=400,
                  xaxis_visible=False, yaxis_visible=False, showlegend=False, margin=dict(l=0,r=0,t=0,b=0))

####################################  LAYOUT  ####################################
app.layout = dbc.Container([
    # Sidebar
    html.Div([
        # Sidebar - Header
        html.Div([
            html.H1("Welcome to the ""OMIX"" buddy 🤖!"),
            html.P("A tool to help pharma R&D teams to perform Differential Expression analysis with transcriptomics data.")
        ], style={"vertical-alignment": "top", "height": 260}),

        # Sidebar - Diseases selection
        html.H2('Select a disease to analyze:'),
        html.Div([
            html.Div(
                dbc.RadioItems(
                    id='disease_input',
                    className='btn-group',
                    inputClassName='btn-check',
                    labelClassName="btn btn-outline-light",
                    labelCheckedClassName="btn btn-light",
                    options=[
                        {"label": "Kawasaki", "value": 1},
                        {"label": "Parkinson", "value": 2},
                        {"label": "Melanoma", "value": 3}
                    ],
                    value=1,
                    style={'width': '100%'}
                ), style={'width': 206}
            )
        ], style={'margin-left': 15, 'margin-right': 15, 'display': 'flex'}),

        # Sidebar - Parameters selection
        html.Br(),
        html.Details([
            html.Summary('Parameters Tab(Click To Show!!!)'),
            html.Div([
                html.Div([
                    html.H2('Number of genes to investigate:'),
                    dcc.Slider(2, 14, 2,
                        id='num_genes_input',
                        value=8,
                    ),
                ]),
                html.Div([
                    html.H2('Patients to investigate:'),
                    dcc.Dropdown(
                        id='patient_selection_input',
                        options=patient_list,
                        value=patient_list_init,
                        multi=True,
                        maxHeight=250,
                        optionHeight=30
                    ),
                ]),
            ], style={'margin-left': 15, 'margin-right': 15, 'margin-top': 30}),
        ]),

        # Sidebar Image
        html.Div(
            html.Img(src='assets/Data_boost_logo.png',
                     style={'margin-left': 15, 'margin-right': 15, 'margin-top': 30, 'width': 150})
        ),
    ], style={'width': 330, 'margin-left': 15, 'margin-right': 15, 'margin-top': 20, 'margin-bottom': 20}
    ),

    # Body
    dbc.Col([
        dbc.Row(
            html.H2('Genes level analysis')
        ),
        dbc.Row([
            dbc.Col(dcc.Graph('graph-1', config={'displaylogo': False})),
            dbc.Col([
                dcc.Graph('graph-2', config={'displaylogo': False}), 
                dcc.RangeSlider(min=-6, max=6, step=1, value=[-2, 2], id='log2fold_input'),
                dcc.Slider(min=0, max=40, step=5, value=15, id='significance_input'),
                ])
        ]),
        dbc.Row(
            html.H2('Patient level analysis')
        ),
        dbc.Row([
            dbc.Col(dcc.Graph('graph-3', config={'displaylogo': False})),
            dbc.Col(dcc.Graph('graph-4', config={'displaylogo': False}))
        ])
    ], style={'width': 1210, 'margin-left': 15, 'margin-top': 15, 'margin-right': 20, 'margin-bottom': 20})  # Ajuster la largeur de la colonne pour les graphiques
], fluid=True, style={'display': 'flex'}, className='dashboard-container')

####################################  I/O  ####################################
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

####################################  Charts Update ####################################
def update_graphs(selected_value, num_genes_input, patient_selection_input, log2fold_input, significance_input):
    
    # Data Loading
    if selected_value == 1 or selected_value == 2 or selected_value == 3:
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

    #########################################  Volcano plot  #########################################
    figure2 = volcano_plot(sigs, list_annotation, log2fold_input, significance_input)

    #########################################  PCA plot  #########################################
    figure3 = PCA_plot(dds, grapher, patient_selection_input)

    #########################################  Clustermap plot  #########################################
    figure4 = clustermap_plot(grapher, list_annotation, patient_selection_input)
    
    return fig, figure2, figure3, figure4

if __name__ == "__main__":
    app.run_server(debug=True)