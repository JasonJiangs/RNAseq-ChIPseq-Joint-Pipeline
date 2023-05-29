import pandas as pd

# Load the differentially expressed genes into a list from a csv
deseq_results = pd.read_csv("data/sig01/deseq2_results_sig_01.csv")

# unique the gene ids
deseq_results['gene_id'] = deseq_results['gene_id'].str.split('|').str[0][1:]

print()