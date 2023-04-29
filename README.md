# RNA-seq & ChIP-seq Joint Analysis Workflow and Pipeline
The repository contains the workflow and pipelines for exploring transcriptional regulatory
function with high-throughout data using RNA-seq and ChIP-seq.

## Introduction

## Workflow


## How to use
### 1. Download the repository
```
git clone
```

### 2. Install the required packages
There is no instruction for installing the required packages, only the versions of the packages are listed.
#### Python packages
```
pip install -r pyrequirements.txt
``` 
#### R packages
```
```
#### Bioinformatics tools with Version
```
FastQC -- 
Bowtie2 -- 
Samtools --
HTSeq --
DESeq2 --
CEAS --
Bedtools --
```

### 3. Configure the configuration file
The configuration file is `config.yaml`, in which you should store your datasource path and result path, 
otherwise the result will be stored in the default path. The datasource path must be defined.
```yaml
# it should take the directory
# mapping-index: mapping indexes for bowtie2, bwa, hisat2, salmon, star
# annotation: gene annotation for hisat2, stringtie, cufflinks (.gtf, .gff, and .bed)
datasource:
  rna-seq:
    dir_path: /path/to/rna-seq/
    files:
      - LNCaP_DHT_RNA_rep1
      - LNCaP_DHT_RNA_rep2
      - LNCaP_Veh_RNA_rep1
      - LNCaP_Veh_RNA_rep2
    paired-end: y
    suffix: .fastq.gz
  chip-seq:
    dir_path: /path/to/chip-seq/
    files:
      - LNCaP_DHT_AR_1
      - LNCaP_DHT_AR_2
      - LNCaP_Veh_AR_1
      - LNCaP_Veh_AR_2
    paired-end: n
    suffix: .fastq.gz
  mapping-index: /path/to/mapping-index
  annotation: /path/to/annotation

resultdestination:
  rna-seq: /path/to/rna-seq-results
  chip-seq: /path/to/chip-seq-results
  joint-analysis: /path/to/joint-analysis-results
  modeling: /path/to/modeling-results
  log: /path/to/log

rnaseq:
  bowtie2:
    -p: 16
    -t: y

chipseq:
  hisat2:
    -p: 16

joint:
  fastqc:
    -f: fastq
  # .sam file conversion, sorting, indexing
  samtools:
    view:
      -b: y
      -S: y
      -@: 16
    sort:
      -@: 16
    index: y

```

### 3. Run the workflow

#### Parameter Clarification
- `-c`: get the configuration file path.
- `-qc`: multiple quality control options.
    -  `o`: Only run the quality control workflow.
    -  `a`: Run the rest of the workflow except for the quality control workflow.
    - `oa`: Run the whole workflow.
- `-m`: choose to do some modeling analysis
    - `c`: Run the clustering analysis.
    - `r`: Run the regression analysis.
    - `cr`: Run the clustering and regression analysis.

Sample command (Whole workflow):
```
cd RNA-seq-ChIP-seq-Joint-Analysis
python main.py -c config.yaml -qc oa -m cr
```

## Result
The default result storage is the `result` folder. You can also change the storage path in the `config.yaml` file.
### Result file structure:
```bash

```

## Reference
