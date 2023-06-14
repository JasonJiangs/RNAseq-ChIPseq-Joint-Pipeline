# RNA-seq & ChIP-seq Joint Analysis Pipeline
*Under development*

## Introduction
The repository contains the workflow and pipelines for exploring transcriptional regulatory
function with high-throughout data using RNA-seq and ChIP-seq.

It is an extensible and modular pipeline combining RNA-seq and ChIP-seq tools, with easy configuration via YAML files and command line arguments, providing users with phase-based execution flexibility and advanced modeling analysis options.

## Features
- **Modular**: The pipeline is designed to be modular, which means that you can choose to run the whole pipeline or just part of it.
- **Extensible**: The pipeline is designed to be extensible, which means that you can add your own tools and integrate them into the pipeline.
- **Easy to use**: The pipeline is designed to be easy to use, which means that you can just run the pipeline with a single command.
- **Phase-based execution**: The pipeline is designed to be phase-based execution, which means that you can choose to run the pipeline with a single phase or multiple phases.
- **Advanced modeling analysis**: The pipeline is designed to be advanced modeling analysis, which means that you can run the clustering and regression analysis with a single command.

## Workflow


## How to use
### 1. Download the repository
```
git clone https://github.com/JasonJiangs/RNAseq-ChIPseq-Joint-Pipeline.git
cd RNA-seq-ChIP-seq-Joint-Analysis
```

### 2. Install the required packages
There is no instruction for installing the required packages, only the versions of the packages are listed.
#### Python packages
We use python 3 for the pipeline.
```
pipreqs . --force
module load anaconda
pip install -r requirements.txt
``` 
#### R packages
```
```
#### Bioinformatics tools
Those tools are required in the pipeline.

Phase 1
- [x] [fastp]()
- [x] [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [x] Use fastqc before and after alignment.
  - [x] Add qc tolerance.
- [x] [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [x] check paired-end or single-end
- [x] [hisat2](http://daehwankimlab.github.io/hisat2/)
  - [x] check paired-end or single-end
- [ ] [multiqc](https://multiqc.info/)

Phase 2
- [x] [samtools](http://www.htslib.org/)
- [x] [MACS2]()
- [x] [stringtie]()
- [x] [prepDE.py](https://ccb.jhu.edu/software/stringtie/dl/prepDE.py)
- [x] [getTPM.py](https://ccb.jhu.edu/software/stringtie/dl/getTPM.py)
- [x] [getPFKM.py](https://ccb.jhu.edu/software/stringtie/dl/getFPKM.py)  
- [ ] [deeptools](https://deeptools.readthedocs.io/en/develop/)
- [ ] [bedtools](https://bedtools.readthedocs.io/en/latest/)

Phase 3


### 3. Configure the configuration file
The configuration file is `config.yaml`, in which you should store your datasource path and result path, 
otherwise the result will be stored in the default path. The datasource path must be defined.
```yaml
# Choose only to do script generation
script-only: y

# it should take the directory
# mapping-index: mapping indexes for bowtie2, bwa, hisat2, salmon, star
# annotation: gene annotation for hisat2, stringtie, cufflinks (.gtf, .gff, and .bed)
datasource:
  rna-seq:
    dir_path: /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/
    files:
      - LNCaP_DHT_RNA_rep1
      - LNCaP_DHT_RNA_rep2
      - LNCaP_Veh_RNA_rep1
      - LNCaP_Veh_RNA_rep2
    # if paired-end is y, then LNCaP_DHT_RNA_rep1 -> LNCaP_DHT_RNA_rep1_1.fastq.gz, LNCaP_DHT_RNA_rep1_2.fastq.gz
    paired-end: y
    suffix: .fastq.gz # fastq.gz supported only for now
  chip-seq:
    dir_path: /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/
    files:
      - LNCaP_DHT_AR_1
      - LNCaP_DHT_AR_2
      - LNCaP_Veh_AR_1
      - LNCaP_Veh_AR_2
    # if paired-end is n, then LNCaP_DHT_AR_1 -> LNCaP_DHT_AR_1.fastq.gz
    paired-end: n
    suffix: .fastq.gz  # fastq.gz supported only for now
  # e.g.: Storing bowtie2 index inside as a subdirectory in name of bowtie2_index.
  mapping-index:
    hisat2: /nv/vol190/zanglab/zw5j/data/index/hisat2_index/hg38
    bowtie2: /nv/vol190/zanglab/zw5j/data/index/bowtie_index/bowtie2/hg38
  annotation: 
    base: /nv/vol190/zanglab/shared/StudentTestData/annotation/

# Root directory for the results
# TODO: add '/' at the end of the path
resultdestination: /scratch/pfq7pm/test_pipeline/proj_result

# slurm configuration
slurm:
  partition: standard
  time: '72:00:00'
  mem: 64G
  cpus-per-task: 16
  A: zanglab

# qc-toleration: number of files that can fail qc before the pipeline stops
qc-toleration:
  warn: 4
  fail: 1

chipseq:
  # read alignment
  bowtie2:
    -p: 16
    -t: y
  macs2:
    -t: DHT
    -c: Veh
    -f: BAM
    -g: hs
    -bdg: y
  prepDE: y
  getTPM: y
  getFPKM: y

rnaseq:
  # read alignment
  hisat2:
    -p: 16
  stringtie:
    -e: y
    -B: y
    -G: hg38_refseq_genes.gtf   
  

joint:
  # quality control
  fastqc:
    -f: fastq
  # .sam file conversion, sorting, indexing
  samtools:
    view:
      -b: y
      -S: y
      -@: 8
    sort:
      -@: 8
    index: y

```

### 4. Run the workflow

#### Parameter Clarification
- `-c`: get the configuration file path
- `-p`: multiple phases you can choose to run.
  - `1`: Run quality control (read alignment and fastqc before and after).
  - `2`: Run quality check like mappability, replicate correlation check, and duplicate rate.
  - `3`: Run rest of the RNA-seq and ChIP-seq process and joint analysis.
- `-m`: choose to do some modeling analysis
  - `c`: Run the clustering analysis.
  - `r`: Run the regression analysis.
  - `cr`: Run the clustering and regression analysis.

#### Instructions
If you choose to run with  `-p 1` and check the quality control report is just fine. You want to keep the following process
running, you can just run with `-p 2`. The same logic when you want to keep run the phase 3 code.
Also, if you want run the whole workflow, you can just run with `-p 123`.

Sample command (full process):
```
python main.py -c config.yaml -p 123 -m cr
```
### 5. Result
You should appoint the storage path, which is `resultdestination` in the `config.yaml` file.

### 6. Generated Folder Structure
It includes result, program log file, shell scripts, and shell logs. 
```text
-bash-4.2$module load tree/1.8.0
-bash-4.2$pwd
/scratch/pfq7pm/test_pipeline
-bash-4.2$tree
.
└── proj_result
    ├── execution_log.txt
    └── result
        ├── chipseq
        │   ├── bowtie2
        │   └── macs2
        ├── joint
        │   └── qc
        │       ├── after
        │       └── before
        ├── rnaseq
        │   ├── deseq2
        │   ├── hisat2
        │   ├── htseq
        │   └── stringtie
        ├── shell_log
        ├── shell_script
        │   ├── bowtie2.sh
        │   ├── fastqc_after.sh
        │   ├── fastqc_before.sh
        │   └── hisat2.sh
        └── total_script.sh
```
Some sample result can be referred to `sample_result` folder.

## Clarification
### Script-only mode
Since the pipeline basically intake the configuration file and decode it into shell files that can run under
slurm environment. You should make sure the environment setting is correct, since the pipeline uses `module load` to
load the required packages.

Otherwise, you can just take the generated scripts as a reference and run it on your own environment.
When you config the `config.yaml` file with `script-only: y`, you should get only script in the result folder - you 
can basically take it as a script generator.


## Acknowledgement
