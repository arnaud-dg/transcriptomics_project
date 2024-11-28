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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

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
    """Create GO enrichment plot using plotly with consistent styling"""
    # Filter for minimum number of genes
    go_data = go_data[go_data.n_genes > 1]
    
    # Get top 20 terms by p-value
    n_terms = 20
    top_GO = go_data.nsmallest(n_terms, 'p_corr')
    
    # Truncate long term names
    top_GO['term_display'] = top_GO['term'].apply(lambda x: x[:27] + '...' if len(x) > 30 else x)
    
    # Define colors for each GO class with updated palette
    class_colors = {
        'biological_process': '#FF881F',  # Match the orange from PCA
        'cellular_component': '#4F6780',  # Match the blue-grey from PCA
        'molecular_function': '#7A9CC6'   # Intermediate shade
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
                    y=class_data['term_display'],
                    mode='markers',
                    name=go_class.replace('_', ' ').title(),
                    marker=dict(
                        color=class_colors[go_class],
                        size=10
                    ),
                    hovertemplate=(
                        "<b>%{y}</b><br>" +
                        "p-value: %{x:.2f}<br>" +
                        "<extra></extra>"
                    )
                )
            )
    
    # Create figure
    figure = go.Figure(data=traces)
    
    # Update layout to match PCA plot style
    figure.update_layout(
        plot_bgcolor='#EAEAF2',
        paper_bgcolor='rgba(0,0,0,0)',
        title={
            'text': f'GO Term Enrichment Analysis - {disease_title}',
            'x': 0.5,
            'y': 0.95,
            'xanchor': 'center',
            'yanchor': 'middle',
            'font': {
                'size': 18,
                'color': '#555555',
                'family': 'Poppins'
            }
        },
        xaxis_title={
            'text': '-log₁₀(Adjusted p-value)',
            'standoff': 16,
            'font': {
                'family': 'Poppins',
                'size': 13,
                'color': 'black'
            }
        },
        yaxis_title={
            'text': 'GO Terms',
            'standoff': 16,
            'font': {
                'family': 'Poppins',
                'size': 13,
                'color': 'black'
            }
        },
        legend={
            'x': 0.5,
            'y': -0.15,
            'xanchor': 'center',
            'yanchor': 'top',
            'orientation': 'h',
            'bgcolor': 'rgba(255, 255, 255, 0)',
            'title_text': 'GO Categories'
        },
        height=480,
        margin=dict(l=10, r=10, t=50, b=70),
        showlegend=True
    )
    
    # Update axes with white background grid
    figure.update_xaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='white',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='white'
    )
    
    figure.update_yaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='white',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='white',
        tickfont=dict(size=10)
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

def PCA_plot(grapher, patient_list, disease_title):

    # Prepare data
    pca_data = grapher[patient_list]
    
    # Standardize the features
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(pca_data.T)
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(data_scaled)
    
    # Create DataFrame for plotting
    pca_df = pd.DataFrame(pca_coords, columns=['PC1', 'PC2'])
    if disease_title == 'Kawasaki Disease':
        pca_df['Condition'] = ['Sick Patients' if 'KD' in p else 'Control Patients' for p in patient_list]
    elif disease_title == 'SARS-CoV-2':
        pca_df['Condition'] = ['Sick Patients' if 'INF' in p else 'Control Patients' for p in patient_list]
    
    # Create axis labels with explained variance percentage
    x_label = f"PC1 (Explained variance {pca.explained_variance_ratio_[0]*100:.2f}%)"
    y_label = f"PC2 (Explained variance {pca.explained_variance_ratio_[1]*100:.2f}%)"
    title_label = 'Principal Component Analysis - ' + disease_title

    # Create the plot
    figure = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Condition',
        color_discrete_map=color_map_patients,
        hover_name=patient_list,
        hover_data={'PC1': ':.2f', 'PC2': ':.2f', 'Condition': False}
    )
    
    # Update layout with transparent background and centered title
    figure.update_layout(
        plot_bgcolor='#EAEAF2', 
        paper_bgcolor='rgba(0,0,0,0)',
        title={
            'text': title_label,
            'x': 0.5,
            'y': 0.95,
            'xanchor': 'center',
            'yanchor': 'middle',
            'font': {
                'size': 18, 
                'color': '#555555', 
                'family': 'Poppins'
            }
        },
        xaxis_title={
            'text': x_label,
            'standoff': 16,
            'font': {
                'family': 'Poppins', 
                'size': 13, 
                'color': 'black'
            }
        },
        yaxis_title={
            'text': y_label,
            'standoff': 16,
            'font': {
                'family': 'Poppins', 
                'size': 13, 
                'color': 'black'
            }
        },
        legend={
            'x': 0.5,
            'y': -0.15,
            'xanchor': 'center',
            'yanchor': 'top',
            'orientation': 'h',
            'bgcolor': 'rgba(255, 255, 255, 0)',
            'title_text': 'Patient Condition'
        },
        template="seaborn",
        height=480,
        margin=dict(l=10, r=10, t=50, b=70)  # Adjusted bottom margin for legend
    )
    
    # Update axes with white background grid
    figure.update_xaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='white',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='white'
    )
    
    figure.update_yaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='white',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='white'
    )

    return figure