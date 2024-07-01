# Importing necessary libraries
# Generic librairies
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import pickle

# Graphical librairies
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns

# Dash related librairies
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bio as dashbio

# Biocomputing related librairies
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
from sanbomics.plots import volcano
import scanpy as sc # PCA library

# Load raw data from CSV file
data = pd.read_csv('../data/raw/GSE178491_KD.csv', sep=",")

# Preprocess data: Set index, drop unnecessary columns, and filter rows
data.set_index('ensembl_gene_id', inplace=True)
data.drop('genename', axis=1, inplace=True)
data = data[data.sum(axis=1) > 0].T

# Convert data types to integer
for col in data.select_dtypes(include=['float64', 'int64']).columns:
    data[col] = data[col].astype(int)

# Create metadata dataframe for samples
sample_list = data.index.tolist()
condition_list = ['Kawa' if item.startswith('KD') else 'Ctrl' for item in sample_list]
metadata = pd.DataFrame({'Sample': sample_list, 'Condition': condition_list}).set_index('Sample')

# Create DESeq2 dataset object and perform differential expression analysis
dds = DeseqDataSet(counts=data, metadata=metadata, design_factors="Condition")
dds.deseq2()
stat_res = DeseqStats(dds, contrast=['Condition', 'Kawa', 'Ctrl'])
stat_res.summary()
res = stat_res.results_df

# Map gene symbols and filter significant results
mapper = id_map(species='human')
res['Symbol'] = res.index.map(mapper.mapper)
res['-log10(pValue)'] = -np.log10(res['pvalue'])
sigs = res[(res['baseMean'] >= 5) & (res['padj'] < 0.05) & (abs(res['log2FoldChange']) > 0.025) & (res['Symbol'].notna())]

# Identification of 20 main significance genes
sigs['Significance'] = np.abs(np.log10(sigs['pvalue']))
sigs['sorter'] = sigs['Significance']*sigs['log2FoldChange']
num_genes = 15
label_df = pd.concat((sigs.sort_values('sorter')[-num_genes:], sigs.sort_values('sorter')[0 : num_genes]))
list_annotation = list(label_df['Symbol'])

# Log-transform normalized counts and filter significant genes
dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
dds_sigs = dds[:, sigs.index]

dds_sub = dds[:, sigs[sigs.Symbol.isin(list_annotation)].index]
grapher = pd.DataFrame(dds_sub.layers['log1p'].T, index=dds_sub.var_names, columns=dds_sub.obs_names)
grapher.index = grapher.index.map(mapper.mapper)

# Save the preprocessed files
# sigs dataframe
sigs.to_csv('../data/preprocessed/sigs.csv', sep=",")
# grapher dataframe
grapher.to_csv('../data/preprocessed/grapher.csv', sep=",")
#dds object
with open('../models/dds_object.pkl', 'wb') as f:
    pickle.dump(dds_sigs, f)