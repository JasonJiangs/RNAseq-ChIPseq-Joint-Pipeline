import pandas as pd

# Load the DESeq2 output into a pd
deseq_results = pd.read_csv("data/sig001/deseq2_results_sig_001_with_loc.csv")
# print length with word
print("Number of differentially expressed genes: ", len(deseq_results))
# recognize up or down regulated genes from the DESeq2 results and assign the value to a new column
# print(deseq_results['log2FoldChange'].head())
deseq_results['up_down'] = deseq_results[' log2FoldChange'].apply(lambda x: 'up' if x > 0 else 'down')
# print the first 5 rows

# rename the firt column to gene_id
# deseq_results.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)

print(deseq_results.head())

# save
deseq_results.to_csv("data/sig001/deseq2_results_sig_001_with_loc_up_down.csv", index=False)



