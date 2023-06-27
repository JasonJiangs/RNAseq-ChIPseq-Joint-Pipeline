import pandas as pd
import numpy as np

# Load the BED file into a pandas DataFrame
df = pd.read_csv('../data/chipseq/macs2/LNCaP_DHT_AR_1_PK_summits.bed', sep='\t', header=None)

# Determine how many chunks you need
num_chunks = len(df) // 5000 + 1

# Split the DataFrame into chunks
chunks = np.array_split(df, num_chunks)

# Save each chunk to a separate file
for i, chunk in enumerate(chunks):
    chunk.to_csv(f'regions_chunk_{i}.bed', sep='\t', index=False, header=False)
