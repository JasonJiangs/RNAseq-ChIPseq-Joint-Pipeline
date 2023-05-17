import pandas as pd

# load the MACS2 output file
df = pd.read_csv('../data/chipseq/macs2/LNCaP_DHT_AR_2_PK_summits.bed', sep='\t', header=None)

# create a new column that concatenates the chromosome, start, and end positions
df['loc'] = df[0].astype(str) + ':' + df[1].astype(str) + '-' + df[2].astype(str)

# count the total number of reads
total_reads = df.shape[0]

# count the number of duplicate reads
dup = df.duplicated('loc')
# dup[0] = True
duplicate_reads_list = df[dup]
duplicate_reads = duplicate_reads_list.shape[0]

# calculate the duplicate rate
duplicate_rate = duplicate_reads / total_reads

print(f'The duplicate rate is: {duplicate_rate}')
