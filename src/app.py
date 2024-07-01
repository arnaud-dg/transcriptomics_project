# Importing necessary libraries
# # Generic packages
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import pickle

# Graphical packages
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns

# Dash related packages
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import dash_bio as dashbio

# Biocomputing related packages
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
from sanbomics.plots import volcano
import scanpy as sc # PCA library

# Data Loading
sigs = pd.read_csv('../data/preprocessed/sigs.csv', sep=",")
grapher = pd.read_csv('../data/preprocessed/grapher.csv', sep=",")
with open('../models/dds_object.pkl', 'rb') as f:
    dds = pickle.load(f)

patient_list = ['KD1','KD2','KD3','KD4','KD5','KD6','KD7','KD8','KD9','KD10','KD11','KD12','KD13','KD14','KD15','KD16','KD17','KD18','KD19','KD20',
                'control1','control2','control3','control4','control5','control6','control7','control8','control9','control10','control11','control12',
                'control13','control14','control15','control16','control17','control18','control19','control20']

# Initializing 
color_map = {
    'Not significant': '#D3D3D3',  # tomato color
    'Upregulated': '#1E90FF',  # tomato color
    'Downregulated': '#FF6347',  # dodger blue color
}

sigs['Significance'] = np.abs(np.log10(sigs['pvalue']))
sigs['sorter'] = sigs['Significance']*sigs['log2FoldChange']
conditions = [
    (sigs['Significance'] <= 20) & ((sigs['log2FoldChange'] <= 2) | (sigs['log2FoldChange'] >= -2)),
    (sigs['Significance'] > 20) & (sigs['log2FoldChange'] > 2),
    (sigs['Significance'] > 20) & (sigs['log2FoldChange'] < -2)
]
sigs['label'] = np.select(conditions, ['Not significant', 'Upregulated', 'Downregulated'], default='Not significant')

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.MINTY, dbc.icons.FONT_AWESOME])

# Define the layout of the app
app.layout = dbc.Container([
    dbc.Row([
        #########################################  Sidebar  #########################################
        dbc.Col([
            html.H1([
                html.Span("Welcome"),
                html.Br(),
                html.Span("to the ""omics"" buddy 🤖!")
            ]),
            dcc.Dropdown(
                id='disease_dropdown',
                options=[{'label': disease, 'value': disease} for disease in ["Kawasaki Disease", "Parkinson Disease", "Melanoma"]],
                value='Kawasaki Disease'
            ),
            dcc.Slider(2, 14, 2,
                id='num_genes',
                value=8,
            ),
            dcc.Dropdown(
                id='patient_selection',
                options=[{'label': patient, 'value': patient} for patient in list(grapher.drop('ensembl_gene_id', axis=1).columns)],
                value=patient_list,
                multi=True
            ),
        ], style={'width': '15%'}),  # Ajuster la largeur de la colonne pour le menu de navigation

        #########################################  Charts x4  #########################################
        dbc.Col([
            dbc.Row([
                dbc.Col(dcc.Graph(id='graph-1', config={'displaylogo': False})),
                dbc.Col(dcc.Graph(id='graph-2', config={'displaylogo': False}))
            ]),

            dbc.Row([
                dbc.Col(dcc.Graph(id='graph-3', config={'displaylogo': False})),
                dbc.Col(dcc.Graph(id='graph-4', config={'displaylogo': False}))
            ])
        ], style={'width': '85%'})
    ])
])

# Define the callback to update the graphs
@app.callback(
    [Output('graph-1', 'figure'),
     Output('graph-2', 'figure'),
     Output('graph-3', 'figure'),
     Output('graph-4', 'figure')],
    [Input('disease_dropdown', 'value'),
     Input('num_genes', 'value'),
     Input('patient_selection', 'value')]
)

def update_graphs(selected_value, num_genes, patient_list):
    
    num_genes = int(num_genes/2)

    label_df = pd.concat(
            (sigs.sort_values('sorter')[-num_genes:],
            sigs.sort_values('sorter')[0 : num_genes]))
    list_annotation = list(label_df['Symbol'])

    #########################################  Volcano plot  #########################################
    volcano = px.scatter(sigs, 
                    x='log2FoldChange', y='Significance', 
                     color='label', color_discrete_map=color_map)

    volcano.update_layout(
        title={'text': 'Volcano plot - Kawasaki Disease','x': 0.5,'xanchor': 'center'},
        xaxis_title={'text': 'log<sub>2</sub> Fold Change','standoff': 20},
        yaxis_title={'text': 'Significance -log<sub>10</sub>(pValue)','standoff': 20},
        template="seaborn",
        # width=800, height=800,
        legend=dict(
            x=1, y=1,
            xanchor='right', yanchor='top',
            bgcolor='rgba(255, 255, 255, 0)'  # Set a semi-transparent white background for the legend
        ),
        legend_title_text='Genes Differential Expression'
    )

    volcano.add_hline(y=1.5, line_width=1, line_dash="dash", line_color="grey")
    volcano.add_vline(x=1, line_width=1, line_dash="dash", line_color="grey")
    volcano.add_vline(x=-1, line_width=1, line_dash="dash", line_color="grey")

    # Add annotations for each Symbol in the DataFrame
    for symbol in list_annotation:
        row = sigs[sigs['Symbol'] == symbol].iloc[0]
        if row['log2FoldChange'] > 0:
            font_color='blue'
            shift_x = 4
            alignment = 'left'
        else:
            font_color='red'
            shift_x = -4
            alignment = 'right'
        volcano.add_annotation(
            x=row['log2FoldChange'], y=row['Significance'],
            text=row['Symbol'], font=dict(size=11, color=font_color),
            xshift=shift_x,
            xanchor=alignment, yanchor='middle',
            showarrow=False, opacity=0.7
            )

    #########################################  Clustermap plot  #########################################
    filtered_grapher = grapher[grapher['ensembl_gene_id'].isin(list_annotation)]
    columns = list(filtered_grapher[patient_list].columns.values)
    rows = list(filtered_grapher['ensembl_gene_id'])
    clustermap=dashbio.Clustergram(
        data=filtered_grapher[patient_list],
        column_labels=columns,
        row_labels=rows,
        color_map='RdBu',
        color_threshold={
            'row': 40,
            'col': 10
        },
        height=600, width=800,
        display_ratio=[0.05, 0.1],
        tick_font={'size':10}
    )

    for trace in clustermap['data']:
        if trace['type'] == 'heatmap':
            trace['showscale'] = False

    fig1 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 1')
    fig2 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 1')
    fig3 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 1')
    fig4 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 2')
    
    return fig4, volcano, fig3, clustermap

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False, port=8050)
