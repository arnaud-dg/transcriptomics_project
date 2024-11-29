# OMIX Dashboard

## Overview
[**The OMIX Buddy**](https://the-omix-buddy-6e7889e55137.herokuapp.com/) is a web-based dashboard designed to assist pharmaceutical R&D teams in performing **Differential Expression Analysis** with transcriptomics data. The tool provides interactive visualizations and analysis capabilities for gene expression data.

---

## Features

- **Disease-specific analysis selector**
- **Interactive parameter adjustment**
- **Four main visualization components**:
  - **Volcano Plot**: Gene expression analysis
  - **GO Plot**: Gene ontology analysis
  - **PCA Plot**: Patient clustering analysis
  - **Clustermap**: Gene expression patterns

---

## Technical Stack

### Dependencies

#### Data Analysis
- pandas
- numpy
- scipy.stats
- pydeseq2
- scanpy
- sanbomics

#### Visualization
- matplotlib
- plotly
- seaborn
- dash-bio

#### Web Framework
- Dash
- Dash Bootstrap Components

---

## Installation and Setup

1. Clone the repository.
2. Install required packages with pip.

### Usage

1. Run the preprocessing script to prepare the data.
2. Launch the dashboard.
3. Access the dashboard through your web browser at http://localhost:8050.

---

## Interactive Features

### Parameter Adjustments
- Significance threshold slider
- Up/Down regulation threshold slider
- Number of genes selector
- Patient selection dropdown

### Disease Selection
- Kawasaki Disease
- Parkinson's Disease
- Small Cell Lung Cancer

---

## Data Processing
The preprocessing pipeline includes:
1. Raw data loading and cleaning.
2. Differential expression analysis with DESeq2.
3. Gene symbol mapping.
4. Statistical filtering and significance calculations.
5. Data normalization and transformation.

---

## Visualization Components

### Volcano Plot
- Highlights differential expression significance.
- Interactive labeling for upregulated and downregulated genes.
- **Color Scheme**:
  - Not significant: Light gray (#D3D3D3)
  - Upregulated: Dodger blue (#1E90FF)
  - Downregulated: Tomato red (#FF6347)

### GO Plot
- Gene Ontology analysis visualization.
- Pathway enrichment display.

### PCA Plot
- Principal Component Analysis for patient clustering.
- Displays variance explained metrics.

### Clustermap
- Hierarchical clustering of gene expression data.
- Heatmap visualization with clustering for samples and genes.

---

## Contributing
For contributions or bug reports, please open an issue or submit a pull request.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

---

Created by Arnaud Duigou
