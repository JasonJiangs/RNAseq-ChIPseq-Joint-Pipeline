import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns
mpl.use('Agg')

# Load data
df = pd.read_csv('../data/rnaseq/gene_fpkm_matrix.csv', index_col=0)

# Log transform data
df_log2 = np.log2(df + 1)

# filter out lowly expressed genes
# Calculate average expression level for each gene
average_expression = df_log2.mean(axis=1)

# Set threshold
threshold = 11

# Remove lowly expressed genes
df_filtered = df_log2[average_expression > threshold]
# shape pf df_filtered: (n_genes, n_samples)
print(df_filtered.shape)

# # Generate a clustermap: Process finished with exit code -1073741571 (0xC00000FD)
# g = sns.clustermap(df_log2, cmap='viridis', figsize=(10, 10))
# g.savefig('clustermap.png')


# Perform k-means clustering
kmeans = KMeans(n_clusters=4, random_state=0, n_init='auto').fit(df_filtered.T)

# Add cluster labels to data frame
df_filtered.T['cluster'] = kmeans.labels_

# Plot heatmap
g = sns.clustermap(df_log2, cmap='viridis', figsize=(10, 10))
g.savefig('clustermap.png')

# Remove the cluster column after plotting
df_log2 = df_log2.drop('cluster', axis=1)
