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
import plotly.graph_objects as go
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

color_map_patients = {
    'Condition_1': '#FF881F',  # tomato color
    'Condition_2': '#4F6780',  # tomato color
}

#########################################  Volcano plot  #########################################
def volcano_plot(df, list_annotation, log2fold_input, significance_input, disease_title):

    figure = px.scatter(df, x='log2FoldChange', y='Significance', color='label', color_discrete_map=color_map, 
                        hover_data={'Symbol': True})

    title_label = 'Volcano plot - ' + disease_title

    figure.update_layout(
        title={'text': title_label,'x': 0.5,'xanchor': 'center',
               'font': {'size': 16, 'color': '#555555', 'family': 'Poppins', 'weight': 'bold'}},
        xaxis_title={'text': 'Regulation, log<sub>2</sub> Fold Change','standoff': 16,
                     'font': {'family': 'Poppins', 'size': 13, 'color': 'black', 'style': 'italic'}},
        yaxis_title={'text': 'Significance, -log<sub>10</sub>(pValue)','standoff': 16,
                     'font': {'family': 'Poppins', 'size': 13, 'color': 'black', 'style': 'italic'}},
        template="seaborn",
        height=480,
        margin=dict(l=10, r=10, t=30, b=0),
        legend=dict(
            x=1, y=1,
            xanchor='right', yanchor='top',
            bgcolor='rgba(255, 255, 255, 0)'  # Set a semi-transparent white background for the legend
        ),
        legend_title_text='Genes Differential Expression',
        paper_bgcolor='rgba(0, 0, 0, 0)',  # Fond transparent
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

#########################################  Gene ontology plot  #########################################
def create_go_plot(go_data, disease_title):
    """Create GO enrichment plot using plotly"""
    # Filter for minimum number of genes
    go_data = go_data[go_data.n_genes > 1]
    
    # Get top 20 terms by p-value
    n_terms = 20
    top_GO = go_data.nsmallest(n_terms, 'p_corr')
    
    # Define colors for each GO class
    class_colors = {
        'biological_process': '#FF6B6B',
        'cellular_component': '#4ECDC4',
        'molecular_function': '#4B7BEC'
    }
    
    # Create traces for each class
    traces = []
    for go_class in class_colors:
        mask = top_GO['class'] == go_class
        if mask.any():
            class_data = top_GO[mask]
            traces.append(
                go.Scatter(
                    x=-np.log10(class_data['p_corr']),
                    y=class_data['term'],
                    mode='markers',
                    name=go_class.replace('_', ' ').title(),
                    marker=dict(
                        color=class_colors[go_class],
                        size=10
                    )
                )
            )
    
    # Create figure
    figure = go.Figure(data=traces)
    
    # Update layout
    figure.update_layout(
        title={'text': f'Top Enriched GO Terms - {disease_title}',
               'x': 0.5,
               'xanchor': 'center',
               'font': {'size': 16, 'color': '#555555', 'family': 'Poppins', 'weight': 'bold'}},
        xaxis_title={'text': '-log10(Corrected p-value)',
                    'standoff': 16,
                    'font': {'family': 'Poppins', 'size': 13, 'color': 'black', 'style': 'italic'}},
        yaxis_title={'text': 'GO Terms',
                    'standoff': 16,
                    'font': {'family': 'Poppins', 'size': 13, 'color': 'black', 'style': 'italic'}},
        height=480,
        template="seaborn",
        showlegend=True,
        legend=dict(
            x=1.05,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255, 255, 255, 0.8)'
        ),
        margin=dict(l=10, r=100, t=30, b=0),
        yaxis=dict(
            tickfont=dict(size=10),
            automargin=True
        )
    )
    
    return figure

#########################################  Clustermap plot  #########################################
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
        display_ratio=[0.1, 0.12],
        tick_font={'size':10}
    )

    figure.update_layout(
        height=470,
        width=570,
        margin=dict(l=20, r=0, t=10, b=10), 
        margin_pad=0)

    for trace in figure['data']:
        if trace['type'] == 'heatmap':
            trace['showscale'] = False

    return figure

#########################################  PCA plot  #########################################
def PCA_plot(dds, grapher, patient_list, disease_title):

    positions = [i for i, x in enumerate(list(grapher.columns)) if x in patient_list]
    dds_subset = dds[positions , :] # sigs[sigs.Symbol.isin(list_interest_genes)].index
    sc.pp.pca(dds_subset)
    X_pca = dds_subset.X
    pca_coords = dds_subset.obsm['X_pca']
    explained_variance_ratio = dds_subset.uns['pca']['variance_ratio']
    # Create axis labels, with explained variance percentage
    x_label = f"PCA1 (Explained variance {explained_variance_ratio[0]*100:.2f}%)"
    y_label = f"PCA2 (Explained variance {explained_variance_ratio[1]*100:.2f}%)"
    title_label = 'Principal Component Analysis - ' + disease_title

    pca_df = pd.DataFrame(pca_coords[:,:2], columns=['PCA1', 'PCA2'])
    figure = px.scatter(pca_df, x='PCA1', y='PCA2', color=list(dds_subset.obs.values.flatten()), color_discrete_map=color_map_patients)
    figure.update_layout(
        title={'text': title_label,'x': 0.5,'xanchor': 'center',
               'font': {'size': 16, 'color': '#555555', 'family': 'Poppins', 'weight': 'bold'}},
        xaxis_title={'text': x_label,'standoff': 16,
                     'font': {'family': 'Poppins', 'size': 13, 'color': 'black', 'style': 'italic'}},
        yaxis_title={'text': y_label,'standoff': 16,
                     'font': {'family': 'Poppins', 'size': 13, 'color': 'black', 'style': 'italic'}},
        legend=dict(
            x=1, y=1,
            xanchor='right', yanchor='top',
            bgcolor='rgba(255, 255, 255, 0)'  # Set a semi-transparent white background for the legend
        ),
        template="seaborn",
        legend_title_text='Patient Condition',
        paper_bgcolor='rgba(0, 0, 0, 0)',
        height=480,
        margin=dict(l=10, r=10, t=30, b=0),
    )

    return figure