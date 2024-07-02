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

color_map = {
    'Not significant': '#D3D3D3',  # tomato color
    'Upregulated': '#1E90FF',  # tomato color
    'Downregulated': '#FF6347',  # dodger blue color
}

#########################################  Volcano plot  #########################################
def volcano_plot(df, list_annotation, log2fold_input, significance_input):
    figure = px.scatter(df, 
                    x='log2FoldChange', y='Significance', 
                        color='label', color_discrete_map=color_map)

    figure.update_layout(
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

    figure.add_hline(y=significance_input, line_width=1, line_dash="dash", line_color="grey")
    figure.add_vline(x=log2fold_input[1], line_width=1, line_dash="dash", line_color="grey")
    figure.add_vline(x=log2fold_input[0], line_width=1, line_dash="dash", line_color="grey")

    # Add annotations for each Symbol in the DataFrame
    for symbol in list_annotation:
        row = df[df['Symbol'] == symbol].iloc[0]
        if row['log2FoldChange'] > 0:
            font_color='blue'
            shift_x = 4
            alignment = 'left'
        else:
            font_color='red'
            shift_x = -4
            alignment = 'right'
            
        figure.add_annotation(
            x=row['log2FoldChange'], y=row['Significance'],
            text=row['Symbol'], font=dict(size=11, color=font_color),
            xshift=shift_x,
            xanchor=alignment, yanchor='middle',
            showarrow=False, opacity=0.7
            )
        
    return figure

#########################################  Volcano plot  #########################################
def clustermap_plot(grapher, list_annotation, patient_list):
    filtered_grapher = grapher[grapher['ensembl_gene_id'].isin(list_annotation)]
    columns = list(filtered_grapher[patient_list].columns.values)
    rows = list(filtered_grapher['ensembl_gene_id'])
    figure=dashbio.Clustergram(
        data=filtered_grapher[patient_list],
        column_labels=columns,
        row_labels=rows,
        color_map='RdBu',
        color_threshold={
            'row': 40,
            'col': 10
        },
        # height=600, width=800,
        display_ratio=[0.05, 0.1],
        tick_font={'size':10}
    )

    for trace in figure['data']:
        if trace['type'] == 'heatmap':
            trace['showscale'] = False

    return figure

#########################################  PCA plot  #########################################
def PCA_plot(dds, grapher, patient_list):

    positions = [i for i, x in enumerate(list(grapher.columns)) if x in patient_list]
    dds_subset = dds[positions , :] # sigs[sigs.Symbol.isin(list_interest_genes)].index
    sc.pp.pca(dds_subset)
    X_pca = dds_subset.X
    pca_coords = dds_subset.obsm['X_pca']
    explained_variance_ratio = dds_subset.uns['pca']['variance_ratio']
    # Create axis labels, with explained variance percentage
    x_label = f"PCA1 (Explained variance {explained_variance_ratio[0]*100:.2f}%)"
    y_label = f"PCA2 (Explained variance {explained_variance_ratio[1]*100:.2f}%)"

    pca_df = pd.DataFrame(pca_coords[:,:2], columns=['PCA1', 'PCA2'])
    figure = px.scatter(pca_df, x='PCA1', y='PCA2', title='Principal Component Analysis - Kawasaki disease', color=list(dds_subset.obs.values.flatten()))
    figure.update_layout(
        xaxis_title=x_label, yaxis_title=y_label,
        legend=dict(
            x=1, y=1,
            xanchor='right', yanchor='top',
            bgcolor='rgba(255, 255, 255, 0)'  # Set a semi-transparent white background for the legend
        ),
        legend_title_text='Patient Condition'
    )

    return figure