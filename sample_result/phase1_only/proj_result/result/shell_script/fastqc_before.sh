#!/bin/bash
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH -A zanglab
#SBATCH --job-name=fastqc_before
#SBATCH --output="/scratch/pfq7pm/test_pipeline/proj_result/result/shell_log/fastqc_before.out"

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
