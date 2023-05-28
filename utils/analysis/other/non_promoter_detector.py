import pandas as pd


def option01(gtf, promoters):
    # Delete all rows from the annotation file as long as its start or end is between promoters[1] and promoters[2]
    for i in range(len(promoters)):
        gtf = gtf[~((gtf[3] >= promoters[1][i]) & (gtf[3] <= promoters[2][i]))]
        gtf = gtf[~((gtf[4] >= promoters[1][i]) & (gtf[4] <= promoters[2][i]))]
        print('round', i)
        print('length of non-promoters', len(gtf))

    # save
    gtf[[0, 3, 4]].to_csv('result/non_promoters_2.bed', sep='\t', header=False, index=False)
    return gtf


def option02(gtf, promoters):
    # Delete all rows from the annotation file as long as its start or end is between promoters[1] and promoters[2]
    # based on the strand of the gene
    for i in range(len(promoters)):
        if promoters[3][i] == '+':
            gtf = gtf[~((gtf[6] == '+') & (gtf[3] >= promoters[1][i]) & (gtf[3] <= promoters[2][i]))]
            gtf = gtf[~((gtf[6] == '+') & (gtf[4] >= promoters[1][i]) & (gtf[4] <= promoters[2][i]))]
        else:
            gtf = gtf[~((gtf[6] == '-') & (gtf[3] >= promoters[1][i]) & (gtf[3] >= promoters[2][i]))]
            gtf = gtf[~((gtf[6] == '-') & (gtf[4] >= promoters[1][i]) & (gtf[4] >= promoters[2][i]))]
        print('round', i, ' with strand', promoters[3][i])
        print('length of non-promoters', len(gtf))

    # save
    gtf[[0, 3, 4]].to_csv('result/non_promoters_strand.bed', sep='\t', header=False, index=False)
    return gtf


# find intersecting piece of gtf
def find_intersect_piece(gtf, promoter):
    # define new gtf
    new_gtfs = pd.DataFrame()
    for i in range(len(promoter)):
        if promoter[3][i] == '+':
            # new gtf = new rows + new gtf
            new_gtfs = new_gtfs.append(gtf[((gtf[6] == '+') & (gtf[3] >= promoter[1][i]) & (gtf[3] <= promoter[2][i])) |
                           ((gtf[6] == '+') & (gtf[4] >= promoter[1][i]) & (gtf[4] <= promoter[2][i]))])
        else:
            new_gtfs = new_gtfs.append(gtf[((gtf[6] == '-') & (gtf[3] >= promoter[1][i]) & (gtf[3] >= promoter[2][i])) |
                           ((gtf[6] == '-') & (gtf[4] >= promoter[1][i]) & (gtf[4] >= promoter[2][i]))])
        print('round', i, ' with strand', promoter[3][i])
        print('length of non-promoters', len(new_gtfs))

    # save
    new_gtfs[[0, 3, 4]].to_csv('result/UD2kb/promoters_regions.bed', sep='\t', header=False, index=False)
    return new_gtfs



# Read GTF file
gtf = pd.read_csv('data/hg38_refseq_genes.gtf', sep='\t', comment='#', header=None)
# print length of gtf
print('length of annotations', len(gtf))

# get promoter regions bed file
promoters = pd.read_csv('result/promoters_regions_UD2kb.bed', sep='\t', header=None)
print('length of promoters', len(promoters))


find_intersect_piece(gtf, promoters)

