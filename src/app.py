# Generic packages
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats

# Graphical packages
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns

# Dash related packages
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, callback
import dash_bio as dashbio

# Biocomputing related packages
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
from sanbomics.plots import volcano
import scanpy as sc # PCA library

# Data Loading
data = pd.read_csv('../data/raw/GSE178491_KD.csv', sep=",")

# Data Preprocessing
data = data.set_index('ensembl_gene_id')
data = data.drop('genename', axis=1)
data = data[data.sum(axis=1) > 0]
data = data.T
for col in data.select_dtypes(include=['float64', 'int64']).columns:
    data[col] = data[col].astype(int)

# Metadata creation
sample_list = list(data.index)
condition_list = ['Kawa' if item.startswith('KD') else 'Ctrl' for item in sample_list]
metadata = pd.DataFrame(zip(sample_list, condition_list), columns = ['Sample', 'Condition'])
metadata = metadata.set_index('Sample')

# Creation of the dds object
dds = DeseqDataSet(counts=data, 
                   metadata=metadata,
                   design_factors="Condition")
dds.deseq2()
stat_res = DeseqStats(dds, contrast = ['Condition', 'Kawa', 'Ctrl'])
stat_res.summary()
res = stat_res.results_df

# Customization of the dds object
mapper = id_map(species = 'human')
res['Symbol'] = res.index.map(mapper.mapper)
res['-log10(pValue)'] = -np.log10(res['pvalue'])
sigs = res[(res.baseMean >= 10) & (res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5) & (res.Symbol.isna() == False)]
dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
dds_sigs = dds[:, sigs.index]
grapher = pd.DataFrame(dds_sigs.layers['log1p'].T, 
                       index=dds_sigs.var_names, columns=dds_sigs.obs_names)

# Initialize the Dash app
app = dash.Dash(__name__)

# Define the layout of the app
app.layout = html.Div([
    html.Div([
        html.H2('Sidebar'),
        dcc.Dropdown(
            id='dropdown',
            options=["Kawasaki Disease", "Parkinson Disease", "Melanoma"],
            value='Kawasaki Disease'
        )
    ], style={'width': '20%', 'display': 'inline-block', 'vertical-align': 'top'}),
    
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
    ], style={'width': '75%', 'display': 'inline-block', 'vertical-align': 'top'})
])

# Define the callback to update the graphs
@app.callback(
    [Output('graph-1', 'figure'),
     Output('graph-2', 'figure'),
     Output('graph-3', 'figure'),
     Output('graph-4', 'figure')],
    [Input('dropdown', 'value')]
)
def update_graphs(selected_value):
    
    # Volcano plot
    fold_change_threshold = 1
    p_value_threshold = 0.05
    conditions = [
        (res['log2FoldChange'] > fold_change_threshold) & (res['pvalue'] < p_value_threshold),
        (res['log2FoldChange'] < -fold_change_threshold) & (res['pvalue'] < p_value_threshold)
    ]
    choices = ['Upregulated', 'Downregulated']
    res['Significance'] = np.select(conditions, choices, default='Not significant')
    fig1 = px.scatter(
        res, 
        x='log2FoldChange', 
        y='-log10(pValue)',
        title='Volcano Plot',
    )

    fig2 = dashbio.Clustergram(
        data=res.loc[rows].values,
        column_labels=columns,
        row_labels=rows,
        color_threshold={
            'row': 250,
            'col': 700
        },
        hidden_labels='row',
        height=800,
        width=700
    )

    fig3 = px.scatter(res, x='baseMean', y='log2FoldChange', title='Scatter Plot')
    fig4 = px.scatter(res, x='baseMean', y='log2FoldChange', title='Scatter Plot')
    
    return fig1, fig2, fig3, fig4

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)
