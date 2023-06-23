import pandas as pd

def find_gene_loc(de_gene_path, save_path):
    print('Finding the gene locations ...')
    # Load the DESeq2 output into a dictionary
    deseq_results = {}
    with open(de_gene_path, "r") as file:
        next(file)  # Skip header line
        for line in file:
            fields = line.strip().split(',')
            gene_id = fields[0].split('|')[0][1:]
            deseq_results[gene_id] = fields[1:]

    # Load the gene locations into a dictionary from GTF file
    gene_locations = {}
    with open("tools/data/public/hg38_refseq_genes.gtf", "r") as file:
        for line in file:
            if not line.startswith('#'):  # Skip comment lines
                fields = line.strip().split('\t')
                chrom = fields[0]
                start = fields[3]
                end = fields[4]
                info = fields[8]
                # This might need to be adjusted based on your file format
                gene_id = info.split(';')[0].split(' ')[1].replace("\"", "")
                gene_locations[gene_id] = (chrom, start, end)

    # get keys of gene_locations as a list
    gene_ids = list(gene_locations.keys())
    # Add the gene locations to the DESeq2 results
    for gene_id in deseq_results:
        if gene_id in gene_ids:
            deseq_results[gene_id].extend(gene_locations[gene_id])

    # # Let's print the first 5 items to check
    # for i, (gene_id, info) in enumerate(deseq_results.items()):
    #     print(f"{gene_id}: {info}")
    #     if i > 5:
    #         break

    # save the results to a file
    with open(save_path, "w") as file:
        file.write("gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,chrom,start,end\n")
        for gene_id, info in deseq_results.items():
            file.write(f"{gene_id},{','.join(info)}\n")

    print('Result is saved to', save_path)


def recognize_up_down_regulated(with_loc_path, with_up_down_path):
    print('Recognizing up or down regulated genes ...')
    # Load the DESeq2 output into a pd
    deseq_results = pd.read_csv(with_loc_path)
    # print length with word
    print("Number of differentially expressed genes: ", len(deseq_results))
    # recognize up or down regulated genes from the DESeq2 results and assign the value to a new column
    # print(deseq_results['log2FoldChange'].head())
    deseq_results['up_down'] = deseq_results['log2FoldChange'].apply(lambda x: 'up' if x > 0 else 'down')
    # print the first 5 rows

    # rename the firt column to gene_id
    # deseq_results.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)

    # save
    deseq_results.to_csv(with_up_down_path, index=False)
    print('Result is saved to', with_up_down_path)


