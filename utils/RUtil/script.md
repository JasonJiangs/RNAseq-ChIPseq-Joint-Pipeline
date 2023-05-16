# Script command line arguments description

This file describes the command line arguments of the scripts in this directory.

### replicate_correlation.R
Take the gene_count_matrix.csv and transcript_count_matrix.csv as input, 
and output the figures and csv files of replicate correlation analysis.

The output files are saved in a newly created result directory in the output file path.
```shell
Rscript replicate_correlation.R args01 args02
```
- args01: input file path
- args02: output file path

### vis_fpkm.R
Take the gene_count_matrix.csv and transcript_count_matrix.csv as input,
