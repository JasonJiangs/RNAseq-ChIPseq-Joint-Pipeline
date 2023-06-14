import sqlite3
import pandas as pd

gtf_path = r'data/hg38_refseq_genes.gtf'

df = pd.read_csv(gtf_path, sep='\t', header=None, comment='#')

# Define the column names
df.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Filter for gene features
# df = df[df['feature'] == 'gene']

# Parse the attribute field
df['gene_id'] = df['attribute'].str.extract('gene_id "([^"]+)')
df['gene_name'] = df['attribute'].str.extract('gene_name "([^"]+)')
df['transcript_id'] = df['attribute'].str.extract('transcript_id "([^"]+)')
df['exon_number'] = df['attribute'].str.extract('exon_number "([^"]+)')
df['exon_id'] = df['attribute'].str.extract('exon_id "([^"]+)')


# Select the columns for the GeneTable and rename them
df = df[['gene_id', 'gene_name', 'chrom', 'feature', 'start', 'end', 'strand', 'transcript_id', 'exon_number', 'exon_id']]
df.columns = ['gene_id', 'gene_name', 'chrom', 'feature', 'start', 'end', 'strand', 'transcript_id', 'exon_number', 'exon_id']

# Connect to the SQLite database (it will be created if it doesn't exist)
conn = sqlite3.connect('GeneTable.db')

# Write the DataFrame to the SQLite database
df.to_sql('GeneTable', conn, if_exists='replace', index=False)