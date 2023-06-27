from tools.promoter_detector import range_detector, TSS_detector
from tools.bedtools import read_bed_file, intersect, subtract, write_bed_file
from tools.deseq2_preprocessor import find_gene_loc, recognize_up_down_regulated
from tools.gene_seeker import promoter_gene_detector, enhancer_gene_detector

def main(log2FoldChange, summit_file=None, range_file_path=None, promoter_region_path=None,
         enhancer_region_path=None, gene_path=None, result_path=None):

    # get the range of regions: UD2kb.bed and UD10kb.bed
    # decipher and find the promoter and enhancer regions from hg38_refseq_genes.gtf
    # based on the first exon of each gene
    # TSS_detector(10000, 10000, '+')
    # TSS_detector(2000, 2000, '+')

    # process with bedtools' logic to find the promoter and enhancer regions
    summits = read_bed_file('tools/data/macs2/LNCaP_DHT_AR_1_PK_summits.bed')
    range_2kb = read_bed_file('result/range/TSS_range2000_2000_+.bed')
    range_10kb = read_bed_file('result/range/TSS_range10000_10000_+.bed')
    print('\n')

    # Intersect the summits with the promoters, return complete ranges from summits file that intersect with range file
    # bedtools intersect -a A.bed -b B.bed -wa
    # AR_promoter_binding_2kb = intersect(summits, range_2kb)  # *****
    # AR_promoter_enhancer_binding_10kb = intersect(summits, range_10kb)  # *****
    print('\n')

    # Subtract the 2kb promoter binding from the 10kb promoter binding
    # AR_enhancer_between_2kb_10kb = subtract(AR_promoter_enhancer_binding_10kb, AR_promoter_binding_2kb)  # *****
    print('\n')

    # Save the results to BED files
    promoter_region_path = 'result/interaction/AR_promoter_binding_UD2kb_np.bed'
    enhancer_region_path = 'result/interaction/AR_enhancer_between_2kb_10kb.bed'
    # write_bed_file(promoter_region_path, AR_promoter_binding_2kb)
    # write_bed_file('result/interaction/AR_promoter_enhancer_binding_UD10kb_np.bed', AR_promoter_enhancer_binding_10kb)
    # write_bed_file(enhancer_region_path, AR_enhancer_between_2kb_10kb)
    # print('\n')

    # preprocess the data for DESeq2
    print('Preprocessing the data for DESeq2, the data is from log2FoldChange='+str(log2FoldChange)+' ...')
    # de_gene_path = 'tools/data/deseq2/log2FoldChange_'+str(log2FoldChange)+'/deseq2_results_sig.csv'
    # with_loc_path = 'tools/data/deseq2/log2FoldChange_'+str(log2FoldChange)+'/deseq2_results_sig_loc.csv'
    with_up_down_path = 'tools/data/deseq2/log2FoldChange_'+str(log2FoldChange)+'/deseq2_results_sig_loc_updown.csv'
    # find_gene_loc(de_gene_path, with_loc_path)
    # recognize_up_down_regulated(with_loc_path, with_up_down_path)
    # print('\n')

    # find the genes in the promoter and enhancer regions
    print('Finding the genes in the promoter and enhancer regions ...')
    print('The result will be saved in result/final/log2FoldChange_'+str(log2FoldChange))
    result_path = 'result/final/log2FoldChange_'+str(log2FoldChange)
    # promoter
    print('\nPromoter ...')
    promoter_gene_detector(with_up_down_path, promoter_region_path, result_path, distance_tolerance=5000)
    # enhancer
    print('\nEnhancer ...')
    enhancer_gene_detector(with_up_down_path, enhancer_region_path, result_path, distance_tolerance=5000)


if __name__ == '__main__':
    # choice: 1, 01, 05, 001
    log2FoldChange = '001'
    main(log2FoldChange)
