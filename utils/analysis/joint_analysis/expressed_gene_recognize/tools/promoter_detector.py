import pandas as pd

def range_detector(upstream, downstream, strand):
    # Read GTF file
    gtf = pd.read_csv('tools/data/public/hg38_refseq_genes.gtf', sep='\t', comment='#', header=None)

    # Filter for exon entries
    gtf = gtf[gtf[2] == 'exon']

    # Parse gene_id from attributes
    gtf['gene_id'] = gtf[8].apply(lambda x: dict(item.split() for item in x.split('; '))['gene_id'])

    # Group by gene_id and get the minimum start coordinate for each group, 6: strand
    gtf_grouped = gtf.groupby('gene_id').agg({0: 'first', 3: 'min', 4: 'max', 6: 'first'}).reset_index()


    # Define promoter regions based on the strand of the gene
    gtf_grouped['promoter_start'] = gtf_grouped.apply(lambda row: row[3] - upstream if row[6] == strand
                                                        else row[4] + downstream, axis=1)
    gtf_grouped['promoter_end'] = gtf_grouped.apply(lambda row: row[3] + downstream if row[6] == strand
                                                        else row[4] - upstream, axis=1)

    # rename columns strand
    gtf_grouped.rename(columns={6: 'strand'}, inplace=True)

    # Write BED file
    result_path = 'result/range/range'+str(upstream)+'_'+str(downstream)+'_'+strand+'.bed'
    gtf_grouped[[0, 'promoter_start', 'promoter_end', 'strand']].to_csv(result_path, sep='\t', header=False,
                                                                        index=False)
    print('Finish recognizing regions upstream='+str(upstream)+' downstream='+str(downstream)+' strand='+strand+' ...')
    print('Result is saved in '+result_path + '\n')


def TSS_detector(upstream, downstream, strand):
    gtf = pd.read_csv('tools/data/public/hg38_refseq_genes.gtf', sep='\t', comment='#', header=None)

    # filter for start codon entries
    gtf = gtf[gtf[2] == 'start_codon']

    # parse gene_id from attributes
    gtf['gene_id'] = gtf[8].apply(lambda x: dict(item.split() for item in x.split('; '))['gene_id'])
    gtf_grouped = gtf.reset_index()

    # group by gene_id and get the minimum start coordinate for each group, 6: strand
    # gtf_grouped = gtf.groupby('gene_id').agg({0: 'first', 3: 'min', 4: 'max', 6: 'first'}).reset_index()
    # delete index column
    del gtf_grouped['index']

    # define promoter regions based on the strand of the gene
    gtf_grouped['range_start'] = gtf_grouped.apply(lambda row: row[3] - upstream if row[6] == strand
                                                        else row[4] + downstream, axis=1)
    gtf_grouped['range_end'] = gtf_grouped.apply(lambda row: row[3] + downstream if row[6] == strand
                                                        else row[4] - upstream, axis=1)

    # rename columns strand
    gtf_grouped.rename(columns={6: 'strand'}, inplace=True)

    # write BED file
    result_path = 'result/range/TSS_range'+str(upstream)+'_'+str(downstream)+'_'+strand+'.bed'

    gtf_grouped[[0, 'range_start', 'range_end', 'strand']].to_csv(result_path, sep='\t', header=False,
                                                                        index=False)
    print('Finish recognizing regions upstream='+str(upstream)+' downstream='+str(downstream)+' strand='+strand+' ...')
    print('Result is saved in '+result_path + '\n')
