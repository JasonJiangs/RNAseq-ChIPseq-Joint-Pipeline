import pandas as pd

# Read GTF file
gtf = pd.read_csv('data/hg38_refseq_genes.gtf', sep='\t', comment='#', header=None)

# Filter for exon entries
gtf = gtf[gtf[2] == 'exon']

# Parse gene_id from attributes
gtf['gene_id'] = gtf[8].apply(lambda x: dict(item.split() for item in x.split('; '))['gene_id'])

# Group by gene_id and get the minimum start coordinate for each group, 6: strand
gtf_grouped = gtf.groupby('gene_id').agg({0: 'first', 3: 'min', 4: 'max', 6: 'first'}).reset_index()

# # save
# gtf_grouped.to_csv('gtf_grouped.csv', index=False)

# Define promoter regions based on the strand of the gene
gtf_grouped['promoter_start'] = gtf_grouped.apply(lambda row: row[3] - 2000 if row[6] == '+' else row[4] + 500, axis=1)
gtf_grouped['promoter_end'] = gtf_grouped.apply(lambda row: row[3] + 500 if row[6] == '+' else row[4] - 2000, axis=1)

# rename columns strand
gtf_grouped.rename(columns={6: 'strand'}, inplace=True)

# Write BED file
gtf_grouped[[0, 'promoter_start', 'promoter_end', 'strand']].to_csv('promoters_strand.bed', sep='\t', header=False, index=False)
