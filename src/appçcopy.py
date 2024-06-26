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

res.to_csv()

# Customization of the dds object
mapper = id_map(species = 'human')
res['Symbol'] = res.index.map(mapper.mapper)
res.to_csv('res_kawasaki.csv', sep=",")