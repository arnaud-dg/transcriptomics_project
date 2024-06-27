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
app = dash.Dash(__name__)

# Define the layout of the app
app.layout = html.Div([
    html.Div([
        html.H2('Navigation menu'),
        dcc.Dropdown(
            id='disease_dropdown',
            options=["Kawasaki Disease", "Parkinson Disease", "Melanoma"],
            value='Kawasaki Disease'
        ),
        dcc.Slider(2, 14, 2,
            id='num_genes',
            value=8,
        ),
        html.Div(id='slider-output-container')
    ], style={'width': '15%', 'display': 'inline-block', 'vertical-align': 'top'}),
    
    html.Div([
        html.Div([
            dcc.Graph(id='graph-1')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            dcc.Graph(id='graph-2')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            dcc.Graph(id='graph-3')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            dcc.Graph(id='graph-4')
        ], style={'width': '48%', 'display': 'inline-block'})
    ], style={'width': '85%', 'display': 'inline-block', 'vertical-align': 'top'})
])

# Define the callback to update the graphs
@app.callback(
    [Output('graph-1', 'figure'),
     Output('graph-2', 'figure'),
     Output('graph-3', 'figure'),
     Output('graph-4', 'figure')],
    [Input('disease_dropdown', 'value'),
     Input('num_genes', 'value')]
)

def update_graphs(selected_value, num_genes):
    
    
    num_genes = int(num_genes/2)

    label_df = pd.concat(
            (sigs.sort_values('sorter')[-num_genes:],
            sigs.sort_values('sorter')[0 : num_genes]))
    list_annotation = list(label_df['Symbol'])

    # Volcano plot
    volcano = px.scatter(sigs, 
                    x='log2FoldChange', y='Significance', 
                     color='label', color_discrete_map=color_map)

    volcano.update_layout(
        title={'text': 'Volcano plot - Kawasaki Disease Differential Expression','x': 0.5,'xanchor': 'center'},
        xaxis_title={'text': 'log<sub>2</sub> Fold Change','standoff': 20},
        yaxis_title={'text': 'Significance -log<sub>10</sub>(pValue)','standoff': 20},
        template="seaborn",
        # width=800, height=800,
        legend=dict(
            x=1, y=1,
            xanchor='right', yanchor='top',
            bgcolor='rgba(255, 255, 255, 0)'  # Set a semi-transparent white background for the legend
        ),
        legend_title_text='Genes DE'
    )

    volcano.add_hline(y=1.5, line_width=1, line_dash="dash", line_color="grey")
    volcano.add_vline(x=1, line_width=1, line_dash="dash", line_color="grey")
    volcano.add_vline(x=-1, line_width=1, line_dash="dash", line_color="grey")

    # Add annotations for each Symbol in the DataFrame
    for symbol in list_annotation:
        row = sigs[sigs['Symbol'] == symbol].iloc[0]
        volcano.add_annotation(
            x=row['log2FoldChange'], y=row['Significance'],
            text=row['Symbol'], font=dict(size=12),
            xanchor='left', yanchor='middle',
            showarrow=False, opacity=0.7
        )

    columns = list(grapher.columns.values)
    rows = list(grapher.index)
    clustermap=dashbio.Clustergram(
        data=grapher,
        column_labels=columns,
        row_labels=rows,
        color_map='RdBu',
        color_threshold={
            'row': 40,
            'col': 7
        },
        # height=800, width=800
    )

    # for trace in clustermap['data']:
    #     if trace['type'] == 'heatmap':
    #         trace['showscale'] = False

    fig1 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 1')
    fig2 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 1')
    fig3 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 1')
    fig4 = px.scatter(sigs, x='baseMean', y='log2FoldChange', title='Scatter Plot 2')
    
    return volcano, fig2, fig3, fig4

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)
