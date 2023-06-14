import pandas as pd

# Load the differentially expressed genes into a list from a csv
deseq_results = pd.read_csv("data/sig1/deseq2_results_sig_with_loc.csv")
# print length with word
print("Number of differentially expressed genes: ", len(deseq_results))


# Load the AR binding sites into lists
with open("data/AR_enhancer_between_2kb_10kb.bed", "r") as file:
    AR_enhancer = file.readlines()
# print length with word
print("Number of AR enhancers: ", len(AR_enhancer))

with open("data/AR_promoter_binding_UD2kb.bed", "r") as file:
    AR_promoter = file.readlines()
# print length with word
print("Number of AR promoters: ", len(AR_promoter))

# Define a function to check overlaps
def check_overlap(gene_list, AR_list):
    overlaps = []
    for i in range(len(gene_list)):
        gene_chr = gene_list.iloc[i, 7]
        gene_start = gene_list.iloc[i, 8]
        gene_end = gene_list.iloc[i, 9]
        for j in range(len(AR_list)):
            AR_fields = AR_list[j].rstrip().split('\t')
            AR_chr = AR_fields[0]
            AR_start = int(AR_fields[1])
            AR_end = int(AR_fields[2])
            if gene_chr == AR_chr:
                if gene_end <= AR_end:
                    overlaps.append((gene_list.iloc[i, 0], AR_list[j]))
                elif gene_start >= AR_start:
                    overlaps.append((gene_list.iloc[i, 0], AR_list[j]))
                else:
                    continue

    # for gene in gene_list:
    #     gene_fields = gene.rstrip().split(',')
    #     gene_chr = gene_fields[0]
    #     gene_start = int(gene_fields[1])
    #     gene_end = int(gene_fields[2])
    #     for AR in AR_list:
    #         AR_fields = AR.rstrip().split('\t')
    #         AR_chr = AR_fields[0]
    #         AR_start = int(AR_fields[1])
    #         AR_end = int(AR_fields[2])
    #         if gene_chr == AR_chr and gene_start <= AR_end and gene_end >= AR_start:
    #             overlaps.append((gene, AR))
    return overlaps


# Check overlaps
enhancer_overlaps = check_overlap(deseq_results, AR_enhancer)
promoter_overlaps = check_overlap(deseq_results, AR_promoter)

# Print results
print("Enhancer overlaps: ", enhancer_overlaps)
print("Number of enhancer overlaps: ", len(enhancer_overlaps))
print("Promoter overlaps: ", promoter_overlaps)
print("Number of promoter overlaps: ", len(promoter_overlaps))
