# Load the DESeq2 output into a dictionary
deseq_results = {}
with open("data/sig1/deseq2_results_sig.csv", "r") as file:
    next(file)  # Skip header line
    for line in file:
        fields = line.strip().split(',')
        gene_id = fields[0].split('|')[0][1:]
        deseq_results[gene_id] = fields[1:]

# Load the gene locations into a dictionary from GTF file
gene_locations = {}
with open("data/hg38_refseq_genes.gtf", "r") as file:
    for line in file:
        if not line.startswith('#'):  # Skip comment lines
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            info = fields[8]
            gene_id = info.split(';')[0].split(' ')[1].replace("\"","")  # This might need to be adjusted based on your file format
            gene_locations[gene_id] = (chrom, start, end)

# get keys of gene_locations as a list
gene_ids = list(gene_locations.keys())
# Add the gene locations to the DESeq2 results
for gene_id in deseq_results:
    if gene_id in gene_ids:
        deseq_results[gene_id].extend(gene_locations[gene_id])

# Let's print the first 5 items to check
for i, (gene_id, info) in enumerate(deseq_results.items()):
    print(f"{gene_id}: {info}")
    if i > 5:
        break

# save the results to a file
with open("data/sig1/deseq2_results_sig_with_loc.csv", "w") as file:
    file.write("gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, chrom, start, end\n")
    for gene_id, info in deseq_results.items():
        file.write(f"{gene_id},{','.join(info)}\n")

# Now we can check for overlaps