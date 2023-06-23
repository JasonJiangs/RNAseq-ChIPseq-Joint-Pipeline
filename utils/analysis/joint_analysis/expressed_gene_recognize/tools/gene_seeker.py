import pandas as pd


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


def promoter_gene_detector(deseq_result_path, promoter_region_path, result_path):
    deseq_results = pd.read_csv(deseq_result_path)
    promoter_region = pd.read_csv(promoter_region_path, sep="\t", header=None, names=['chrom', 'start', 'end'])
    print("Number of promoters regions: ", len(promoter_region))
    # divide up and down regulated genes
    promoter_overlaps_up = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                                 'promoter_chrom', 'promoter_start', 'promoter_end',
                                                 'start_distance', 'end_distance', 'overlap'])
    promoter_overlaps_down = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                                   'promoter_chrom', 'promoter_start', 'promoter_end',
                                                   'start_distance', 'end_distance', 'overlap'])

    # Initialize empty dictionaries to hold genes associated with each enhancer and promoter
    promoter_associated_genes = {i: [] for i in range(len(promoter_region))}

    # Iterate over differentially expressed genes
    for idx, row in deseq_results.iterrows():
        gene_loc = (row['chrom'], row['start'], row['end'])

        # Check if gene is within distance of each promoter
        for i, promoter_row in promoter_region.iterrows():
            promoter_loc = (promoter_row['chrom'], promoter_row['start'], promoter_row['end'])
            if is_overlapping(gene_loc, promoter_loc):
                # collect the information for the overlap by up or down regulated genes
                if row['log2FoldChange'] > 0:
                    promoter_overlaps_up.loc[len(promoter_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'],
                                                                           row['end'],
                                                                           promoter_row['chrom'], promoter_row['start'],
                                                                           promoter_row['end'],
                                                                           promoter_row['start'] - row['start'],
                                                                           promoter_row['end'] - row['end'], 'Yes']
                    continue
                else:
                    promoter_overlaps_down.loc[len(promoter_overlaps_down)] = [row['gene_id'], row['chrom'],
                                                                               row['start'], row['end'],
                                                                               promoter_row['chrom'],
                                                                               promoter_row['start'],
                                                                               promoter_row['end'],
                                                                               promoter_row['start'] - row['start'],
                                                                               promoter_row['end'] - row['end'], 'Yes']
                    continue

            if is_within_distance(gene_loc, promoter_loc, 10000):  # adjust distance as needed
                promoter_associated_genes[i].append(row['gene_id'])
                # collect the information for the overlap by up or down regulated genes
                if row['log2FoldChange'] > 0:
                    promoter_overlaps_up.loc[len(promoter_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'],
                                                                           row['end'],
                                                                           promoter_row['chrom'], promoter_row['start'],
                                                                           promoter_row['end'],
                                                                           promoter_row['start'] - row['start'],
                                                                           promoter_row['end'] - row['end'], 'No']
                else:
                    promoter_overlaps_down.loc[len(promoter_overlaps_down)] = [row['gene_id'], row['chrom'],
                                                                               row['start'], row['end'],
                                                                               promoter_row['chrom'],
                                                                               promoter_row['start'],
                                                                               promoter_row['end'],
                                                                               promoter_row['start'] - row['start'],
                                                                               promoter_row['end'] - row['end'], 'No']

    # Print the number of enhancers and promoters associated with at least one gene
    print("Number of promoter regions that associated with at least one gene: ",
          sum([bool(x) for x in promoter_associated_genes.values()]))

    # save
    promoter_overlaps_up.to_csv(result_path+"/promoter_overlaps_up.csv", index=False)
    print("Number of up regulated genes associated with at least one promoter: ", len(promoter_overlaps_up))
    promoter_overlaps_down.to_csv(result_path+"/promoter_overlaps_down.csv", index=False)
    print("Number of down regulated genes associated with at least one promoter: ", len(promoter_overlaps_down))





def enhancer_gene_detector(deseq_result_path, enhancer_region_path, result_path):
    deseq_results = pd.read_csv(deseq_result_path)
    # Load the AR binding sites into dataframes
    enhancer_region = pd.read_csv(enhancer_region_path, sep="\t", header=None, names=['chrom', 'start', 'end'])
    print("Number of enhancers regions: ", len(enhancer_region))
    # divide up and down regulated genes
    enhancer_overlaps_up = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                                 'enhancer_chrom', 'enhancer_start', 'enhancer_end',
                                                 'start_distance', 'end_distance', 'overlap'])
    enhancer_overlaps_down = pd.DataFrame(columns=['gene_id', 'chrom', 'start', 'end',
                                                   'enhancer_chrom', 'enhancer_start', 'enhancer_end',
                                                   'start_distance', 'end_distance', 'overlap'])

    # Initialize empty dictionaries to hold genes associated with each enhancer and promoter
    enhancer_associated_genes = {i: [] for i in range(len(enhancer_region))}

    # Iterate over differentially expressed genes
    for idx, row in deseq_results.iterrows():
        gene_loc = (row['chrom'], row['start'], row['end'])

        # Check if gene is within distance of each enhancer
        for i, enhancer_row in enhancer_region.iterrows():
            enhancer_loc = (enhancer_row['chrom'], enhancer_row['start'], enhancer_row['end'])
            if is_overlapping(gene_loc, enhancer_loc):
                # collect the information for the overlap by up or down regulated genes
                if row['log2FoldChange'] > 0:
                    enhancer_overlaps_up.loc[len(enhancer_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'],
                                                                           row['end'],
                                                                           enhancer_row['chrom'], enhancer_row['start'],
                                                                           enhancer_row['end'],
                                                                           enhancer_row['start'] - row['start'],
                                                                           enhancer_row['end'] - row['end'], 'Yes']
                    continue
                else:
                    enhancer_overlaps_down.loc[len(enhancer_overlaps_down)] = [row['gene_id'], row['chrom'],
                                                                               row['start'], row['end'],
                                                                               enhancer_row['chrom'],
                                                                               enhancer_row['start'],
                                                                               enhancer_row['end'],
                                                                               enhancer_row['start'] - row['start'],
                                                                               enhancer_row['end'] - row['end'], 'Yes']
                    continue

            if is_within_distance(gene_loc, enhancer_loc, 10000):  # adjust distance as needed
                enhancer_associated_genes[i].append(row['gene_id'])
                # collect the information for the overlap by up or down regulated genes
                if row['log2FoldChange'] > 0:
                    enhancer_overlaps_up.loc[len(enhancer_overlaps_up)] = [row['gene_id'], row['chrom'], row['start'],
                                                                           row['end'],
                                                                           enhancer_row['chrom'], enhancer_row['start'],
                                                                           enhancer_row['end'],
                                                                           enhancer_row['start'] - row['start'],
                                                                           enhancer_row['end'] - row['end'], 'No']
                else:
                    enhancer_overlaps_down.loc[len(enhancer_overlaps_down)] = [row['gene_id'], row['chrom'],
                                                                               row['start'], row['end'],
                                                                               enhancer_row['chrom'],
                                                                               enhancer_row['start'],
                                                                               enhancer_row['end'],
                                                                               enhancer_row['start'] - row['start'],
                                                                               enhancer_row['end'] - row['end'], 'No']

    # Print the number of enhancers and promoters associated with at least one gene
    print("Number of enhancer regions that associated with at least one gene: ",
          sum([bool(x) for x in enhancer_associated_genes.values()]))

    # save
    enhancer_overlaps_up.to_csv(result_path+"/enhancer_overlaps_up.csv", index=False)
    print("Number of up regulated genes associated with at least one enhancer: ", len(enhancer_overlaps_up))
    enhancer_overlaps_down.to_csv(result_path+"/enhancer_overlaps_down.csv", index=False)
    print("Number of down regulated genes associated with at least one enhancer: ", len(enhancer_overlaps_down))



