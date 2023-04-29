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
