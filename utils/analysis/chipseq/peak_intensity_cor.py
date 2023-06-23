import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('TkAgg')

# Load the peak data for two replicates
# rep1 = pd.read_csv('../data/chipseq/macs2/LNCaP_DHT_AR_1_PK_summits.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score'])
# rep2 = pd.read_csv('../data/chipseq/macs2/LNCaP_DHT_AR_2_PK_summits.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score'])
rep1 = pd.read_csv('../data/chipseq/macs2/LNCaP_DHT_AR_1_PK_peaks.narrowPeak', sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'])
rep2 = pd.read_csv('../data/chipseq/macs2/LNCaP_DHT_AR_2_PK_peaks.narrowPeak', sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'])


# Merge the two dataframes on the peak coordinates
merged = pd.merge(rep1, rep2, on=['chr', 'start'], suffixes=('_rep1', '_rep2'))
# merged2 = pd.merge(rep3, rep4, on=['chr', 'start', 'end'], suffixes=('_rep1', '_rep2'))

# Calculate Pearson and Spearman correlation
pearson_corr, _ = pearsonr(merged['score_rep1'], merged['score_rep2'])
spearman_corr, _ = spearmanr(merged['score_rep1'], merged['score_rep2'])

print(f'Pearson correlation: {pearson_corr}')
print(f'Spearman correlation: {spearman_corr}')

# Create scatter plot
plt.figure(figsize=(10, 10))
plt.scatter(merged['score_rep1'], merged['score_rep2'], alpha=0.5)
plt.title('Correlation of Peak Intensities Between Replicates')
plt.xlabel('Peak Intensity Replicate 1')
plt.ylabel('Peak Intensity Replicate 2')
# plt.show()
# save
plt.savefig('result/peak_intensity_cor.png')


# Create histogram
plt.figure(figsize=(10, 10))
plt.hist(merged['score_rep1'], bins=100, alpha=0.5, label='Replicate 1')
plt.hist(merged['score_rep2'], bins=100, alpha=0.5, label='Replicate 2')
plt.title('Peak Intensity Distribution')
plt.xlabel('Peak Intensity')
plt.ylabel('Frequency')
plt.legend()
# plt.show()
# save
plt.savefig('result/peak_intensity_cor_hist.png')

# Create boxplot
plt.figure(figsize=(10, 10))
plt.boxplot([merged['score_rep1'], merged['score_rep2']], labels=['Replicate 1', 'Replicate 2'])
plt.title('Peak Intensity Distribution')
plt.ylabel('Peak Intensity')
# plt.show()
# save
plt.savefig('result/peak_intensity_cor_boxplot.png')




