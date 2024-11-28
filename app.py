# # Generic packages
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import pickle
import random
import os
from pathlib import Path

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

# Initialisation du chemin courant avec pathlib
current_dir = Path(__file__).resolve()
data_dir = current_dir.parent / 'data' / 'preprocessed'
models_dir = current_dir.parent / 'models'
assets_dir = current_dir.parent / 'assets'

app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], assets_folder= str(current_dir.parent / 'assets'))
server = app.server

# Initialize disease data dictionaries
disease_data = {
    1: {
        'sigs': pd.read_csv(str(data_dir / 'sigs_kawasaki.csv'), sep=","),
        'grapher': pd.read_csv(str(data_dir / 'grapher_kawasaki.csv'), sep=","),
        'go': pd.read_csv(str(data_dir / 'go_kawasaki.csv'), sep=","),
        'title': 'Kawasaki Disease',
        'patient_list':['KD1','KD2','KD3','KD4','KD5','KD6','KD7','KD8','KD9','KD10','KD11','KD12','KD13','KD14','KD15','KD16','KD17','KD18','KD19','KD20',
                'control1','control2','control3','control4','control5','control6','control7','control8','control9','control10','control11','control12','control13'
                ,'control14','control15','control16'],
        'patient_list_init':['KD1','KD2','KD3','KD4','KD5','KD6','KD7','KD8','KD9','KD10','KD11','KD12','KD13','KD14','KD15','KD16','KD17','KD18','KD19','KD20',
                'control1','control2','control3','control4','control5','control6','control7','control8','control9','control10']
    },
    2: {
        'sigs': pd.read_csv(str(data_dir / 'sigs_sars-cov2.csv'), sep=","),
        'grapher': pd.read_csv(str(data_dir / 'grapher_sars-cov2.csv'), sep=","),
        'go': pd.read_csv(str(data_dir / 'go_sars-cov2.csv'), sep=","),
        'title': 'SARS-CoV-2',
        'patient_list':['0H_1_INF','0H_2_INF','0H_3_INF','0H_4_MOCK','0H_5_MOCK','0H_6_MOCK','1H_1_INF','1H_2_INF','1H_3_INF','1H_4_MOCK','1H_5_MOCK','1H_6_MOCK'],
        'patient_list_init':['0H_1_INF','0H_2_INF','0H_3_INF','0H_4_MOCK','0H_5_MOCK','0H_6_MOCK','1H_1_INF','1H_2_INF','1H_3_INF','1H_4_MOCK','1H_5_MOCK','1H_6_MOCK']
    }
}

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
            html.P("A tool to help pharma R&D teams to perform D.E.A., i.e. Differential Expression Analysis, with transcriptomics data.")
        ], style={"vertical-alignment": "top", "height": 210}),

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
                        {"label": "SARS cov2", "value": 2},
                    ],
                    value=1,
                    style={'width': '100%'}
                ), style={'width': 206}
            )
        ], style={'margin-left': 15, 'margin-right': 15, 'display': 'flex'}),

        # Sidebar - Parameters selection
        html.Details([
            html.Summary('Parameters Tab'),
            html.Div([
                html.Div([
                    html.H2('Significance threshold:'),
                    dcc.RangeSlider(min=-6, max=6, step=1, value=[-2, 2], id='log2fold_input'),
                    html.H2('Up-Down regulation threshold:'),
                    dcc.Slider(min=0, max=40, step=5, value=15, id='significance_input'),
                ]),
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
                        options=[],
                        value=[],
                        multi=True,
                        maxHeight=250,
                        optionHeight=30
                    ),
                ])
            ], style={'margin-left': 15, 'margin-right': 15, 'margin-top': 10}),
        ], style={'margin-top': 15}),

    ], style={'width': 330, 'margin-left': 15, 'margin-right': 15, 'margin-top': 20, 'margin-bottom': 20}
    ),

    # Body
    dbc.Col([
        dbc.Row(html.H3('Genes level analysis'), style={'margin-top': 15}),
        dbc.Row([
            dbc.Col(dcc.Graph('graph-1', config={'displaylogo': False})),
            dbc.Col(dcc.Graph('graph-2', config={'displaylogo': False}))
        ], style={'max-height': '490px'}),
        dbc.Row(html.H3('Patients level analysis'), style={'margin-top': 15}),
        dbc.Row([
            dbc.Col(dcc.Graph('graph-3', config={'displaylogo': False})),
            dbc.Col(dcc.Graph('graph-4', config={'displaylogo': False}))
        ], style={'max-height': '500px'})
    ], style={'width': 1210, 'margin-left': 0, 'margin-top': 0, 'margin-right': 20, 'margin-bottom': 20})  # Ajuster la largeur de la colonne pour les graphiques
], fluid=True, style={'display': 'flex'}, className='dashboard-container')


####################################  CALLBACKS List  ####################################
# Callback to update patient list dropdown based on disease selection
@app.callback(
    [Output('patient_selection_input', 'options'),
     Output('patient_selection_input', 'value')],
    [Input('disease_input', 'value')]
)
def update_patient_list(disease_input):
    current_data = disease_data[disease_input]
    options = [{'label': patient, 'value': patient} for patient in current_data['patient_list']]
    return options, current_data['patient_list_init']


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
def update_graphs(disease_input, num_genes_input, patient_selection_input, log2fold_input, significance_input):

    current_data = disease_data[disease_input]
    sigs = current_data['sigs']
    grapher = current_data['grapher']
    go_df = current_data['go']
    disease_title = current_data['title']

    # Replace zero p-values with the smallest non-zero p-value
    min_non_zero_pvalue = sigs['pvalue'][sigs['pvalue'] > 0].min()
    sigs.loc[sigs['pvalue'] == 0, 'pvalue'] = min_non_zero_pvalue

    # Calculate significance and sorter columns
    sigs['Significance'] = np.abs(np.log10(sigs['pvalue']))
    sigs['sorter'] = sigs['Significance'] * sigs['log2FoldChange']

    # Apply significance conditions
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

    #########################################  GO plot  #########################################
    figure1 = volcano_plot(sigs, list_annotation, log2fold_input, significance_input, disease_title)
    
    #########################################  Volcano plot  #########################################
    figure2 = create_go_plot(go_df, disease_title)

    #########################################  PCA plot  #########################################
    figure3 = PCA_plot(grapher, patient_selection_input, disease_title) 

    #########################################  Clustermap plot  #########################################
    figure4 = clustermap_plot(grapher, list_annotation, patient_selection_input)
    
    return figure1, figure2, figure3, figure4

if __name__ == "__main__":
    app.run_server(debug=True)