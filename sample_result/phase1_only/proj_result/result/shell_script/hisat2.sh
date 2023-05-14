#!/bin/bash
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH -A zanglab
#SBATCH --job-name=hisat2_result
#SBATCH --output="/scratch/pfq7pm/test_pipeline/proj_result/result/shell_log/hisat2_result.out"

module load gcc
module load hisat2

for i in LNCaP_DHT_RNA_rep1 LNCaP_DHT_RNA_rep2 LNCaP_Veh_RNA_rep1 LNCaP_Veh_RNA_rep2
do
hisat2 -t -p 16 -x /nv/vol190/zanglab/zw5j/data/index/hisat2_index/hg38 -1 /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i"_1.fastq.gz -2 /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i"_2.fastq.gz -S /scratch/pfq7pm/test_pipeline/proj_result/result/rnaseq/hisat2/"$i".sam
done
