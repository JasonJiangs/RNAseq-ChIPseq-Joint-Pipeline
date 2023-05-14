#!/bin/bash
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH -A zanglab
#SBATCH --job-name=total_task
#SBATCH --output=total_log.out
module load gcc
module load fastqc

for i in LNCaP_DHT_AR_1 LNCaP_DHT_AR_2 LNCaP_Veh_AR_1 LNCaP_Veh_AR_2
do
fastqc -f fastq -o /scratch/pfq7pm/test_pipeline/proj_result/result/joint/qc/before/ /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i".fastq.gz
done

for i in LNCaP_DHT_RNA_rep1_1 LNCaP_DHT_RNA_rep1_2 LNCaP_DHT_RNA_rep2_1 LNCaP_DHT_RNA_rep2_2 LNCaP_Veh_RNA_rep1_1 LNCaP_Veh_RNA_rep1_2 LNCaP_Veh_RNA_rep2_1 LNCaP_Veh_RNA_rep2_2
do
fastqc -f fastq -o /scratch/pfq7pm/test_pipeline/proj_result/result/joint/qc/before/ /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i".fastq.gz
done

module load bowtie2

for i in LNCaP_DHT_AR_1 LNCaP_DHT_AR_2 LNCaP_Veh_AR_1 LNCaP_Veh_AR_2
do
bowtie2 -p 16 -x /nv/vol190/zanglab/zw5j/data/index/bowtie_index/bowtie2/hg38 -U /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i".fastq.gz -S /scratch/pfq7pm/test_pipeline/proj_result/result/chipseq/bowtie2/"$i".sam
done

module load hisat2

for i in LNCaP_DHT_RNA_rep1 LNCaP_DHT_RNA_rep2 LNCaP_Veh_RNA_rep1 LNCaP_Veh_RNA_rep2
do
hisat2 -t -p 16 -x /nv/vol190/zanglab/zw5j/data/index/hisat2_index/hg38 -1 /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i"_1.fastq.gz -2 /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i"_2.fastq.gz -S /scratch/pfq7pm/test_pipeline/proj_result/result/rnaseq/hisat2/"$i".sam
done
module load gcc
module load fastqc

for i in LNCaP_DHT_AR_1 LNCaP_DHT_AR_2 LNCaP_Veh_AR_1 LNCaP_Veh_AR_2
do
fastqc -f bam -o /scratch/pfq7pm/test_pipeline/proj_result/result/joint/qc/after/ /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i".fastq.gz
done

for i in LNCaP_DHT_RNA_rep1_1 LNCaP_DHT_RNA_rep1_2 LNCaP_DHT_RNA_rep2_1 LNCaP_DHT_RNA_rep2_2 LNCaP_Veh_RNA_rep1_1 LNCaP_Veh_RNA_rep1_2 LNCaP_Veh_RNA_rep2_1 LNCaP_Veh_RNA_rep2_2
do
fastqc -f bam -o /scratch/pfq7pm/test_pipeline/proj_result/result/joint/qc/after/ /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i".fastq.gz
done

