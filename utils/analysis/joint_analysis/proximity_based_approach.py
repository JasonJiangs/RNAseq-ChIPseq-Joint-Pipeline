import pandas as pd

# Load the differentially expressed genes into a list from a csv
deseq_results = pd.read_csv("data/sig001/deseq2_results_sig_001_with_loc_up_down.csv")
print("Number of differentially expressed genes: ", len(deseq_results))

# Load the AR binding sites into dataframes
AR_enhancer = pd.read_csv("data/AR_enhancer_between_2kb_10kb.bed", sep="\t", header=None, names=['chrom', 'start', 'end'])
print("Number of AR enhancers: ", len(AR_enhancer))

AR_promoter = pd.read_csv("data/AR_promoter_binding_UD2kb.bed", sep="\t", header=None, names=['chrom', 'start', 'end'])
print("Number of AR promoters: ", len(AR_promoter))

# initialize two dataframes to hold the enhancer and promoter overlaps for up and down regulated genes respectively
enhancer_overlaps_up = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                             'enhancer_chrom', 'enhancer_start', 'enhancer_end',
                                             'start_distance', 'end_distance', 'overlap'])
enhancer_overlaps_down = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                               'enhancer_chrom', 'enhancer_start', 'enhancer_end',
                                               'start_distance', 'end_distance', 'overlap'])
promoter_overlaps_up = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                             'promoter_chrom', 'promoter_start', 'promoter_end',
                                             'start_distance', 'end_distance', 'overlap'])
promoter_overlaps_down = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                               'promoter_chrom', 'promoter_start', 'promoter_end',
                                               'start_distance', 'end_distance', 'overlap'])

# Define a function to check if two regions overlap
def is_within_distance(gene_loc, binding_loc, distance):
    gene_chrom, gene_start, gene_end = gene_loc
    binding_chrom, binding_start, binding_end = binding_loc
    return gene_chrom == binding_chrom and (gene_start - distance <= binding_start <= gene_start + distance or gene_end - distance <= binding_end <= gene_end + distance)

def is_overlapping(gene_loc, binding_loc):
    gene_chrom, gene_start, gene_end = gene_loc
    binding_chrom, binding_start, binding_end = binding_loc
    return gene_chrom == binding_chrom and ((gene_start <= binding_start <= gene_end or gene_start <= binding_end <= gene_end) or
                                            (binding_start <= gene_start <= binding_end or binding_start <= gene_end <= binding_end))

# Initialize empty dictionaries to hold genes associated with each enhancer and promoter
enhancer_associated_genes = {i: [] for i in range(len(AR_enhancer))}
promoter_associated_genes = {i: [] for i in range(len(AR_promoter))}

# Iterate over differentially expressed genes
for idx, row in deseq_results.iterrows():
    gene_loc = (row['chrom'], row['start'], row['end'])

    # Check if gene is within distance of each enhancer
    for i, enhancer_row in AR_enhancer.iterrows():
        enhancer_loc = (enhancer_row['chrom'], enhancer_row['start'], enhancer_row['end'])
        if is_overlapping(gene_loc, enhancer_loc):
            # collect the information for the overlap by up or down regulated genes
            if row['log2FoldChange'] > 0:
                enhancer_overlaps_up.loc[len(enhancer_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                       enhancer_row['chrom'], enhancer_row['start'], enhancer_row['end'],
                                                                       enhancer_row['start'] - row['start'], enhancer_row['end'] - row['end'], 'Yes']
            else:
                enhancer_overlaps_down.loc[len(enhancer_overlaps_down)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                           enhancer_row['chrom'], enhancer_row['start'], enhancer_row['end'],
                                                                           enhancer_row['start'] - row['start'], enhancer_row['end'] - row['end'], 'Yes']

        if is_within_distance(gene_loc, enhancer_loc, 10000):  # adjust distance as needed
            enhancer_associated_genes[i].append(row['gene_id'])
            # collect the information for the overlap by up or down regulated genes
            if row['log2FoldChange'] > 0:
                enhancer_overlaps_up.loc[len(enhancer_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                       enhancer_row['chrom'], enhancer_row['start'], enhancer_row['end'],
                                                                       enhancer_row['start'] - row['start'], enhancer_row['end'] - row['end'], 'No']
            else:
                enhancer_overlaps_down.loc[len(enhancer_overlaps_down)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                           enhancer_row['chrom'], enhancer_row['start'], enhancer_row['end'],
                                                                           enhancer_row['start'] - row['start'], enhancer_row['end'] - row['end'], 'No']

    # Check if gene is within distance of each promoter
    for i, promoter_row in AR_promoter.iterrows():
        promoter_loc = (promoter_row['chrom'], promoter_row['start'], promoter_row['end'])
        if is_overlapping(gene_loc, promoter_loc):
            # collect the information for the overlap by up or down regulated genes
            if row['log2FoldChange'] > 0:
                promoter_overlaps_up.loc[len(promoter_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                       promoter_row['chrom'], promoter_row['start'], promoter_row['end'],
                                                                       promoter_row['start'] - row['start'], promoter_row['end'] - row['end'], 'Yes']
            else:
                promoter_overlaps_down.loc[len(promoter_overlaps_down)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                           promoter_row['chrom'], promoter_row['start'], promoter_row['end'],
                                                                           promoter_row['start'] - row['start'], promoter_row['end'] - row['end'], 'Yes']

        if is_within_distance(gene_loc, promoter_loc, 10000):  # adjust distance as needed
            promoter_associated_genes[i].append(row['gene_id'])
            # collect the information for the overlap by up or down regulated genes
            if row['log2FoldChange'] > 0:
                promoter_overlaps_up.loc[len(promoter_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                       promoter_row['chrom'], promoter_row['start'], promoter_row['end'],
                                                                       promoter_row['start'] - row['start'], promoter_row['end'] - row['end'], 'No']
            else:
                promoter_overlaps_down.loc[len(promoter_overlaps_down)] = [row['gene_id'], row['chrom'], row['start'], row['end'],
                                                                           promoter_row['chrom'], promoter_row['start'], promoter_row['end'],
                                                                           promoter_row['start'] - row['start'], promoter_row['end'] - row['end'], 'No']

# Print the number of enhancers and promoters associated with at least one gene
print("Number of enhancers associated with at least one gene: ", sum([bool(x) for x in enhancer_associated_genes.values()]))
print("Number of promoters associated with at least one gene: ", sum([bool(x) for x in promoter_associated_genes.values()]))

# save
enhancer_overlaps_up.to_csv("result/sig_001/enhancer_overlaps_up.csv", index=False)
enhancer_overlaps_down.to_csv("result/sig_001/enhancer_overlaps_down.csv", index=False)
promoter_overlaps_up.to_csv("result/sig_001/promoter_overlaps_up.csv", index=False)
promoter_overlaps_down.to_csv("result/sig_001/promoter_overlaps_down.csv", index=False)


