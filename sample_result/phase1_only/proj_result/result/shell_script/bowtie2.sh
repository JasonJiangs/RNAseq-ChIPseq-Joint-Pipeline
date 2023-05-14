#!/bin/bash
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH -A zanglab
#SBATCH --job-name=bowtie2_result
#SBATCH --output="/scratch/pfq7pm/test_pipeline/proj_result/result/shell_log/bowtie2_result.out"

module load gcc
module load bowtie2

for i in LNCaP_DHT_AR_1 LNCaP_DHT_AR_2 LNCaP_Veh_AR_1 LNCaP_Veh_AR_2
do
bowtie2 -p 16 -x /nv/vol190/zanglab/zw5j/data/index/bowtie_index/bowtie2/hg38 -U /nv/vol190/zanglab/shared/StudentTestData/rawdata_LNCaP_AR/"$i".fastq.gz -S /scratch/pfq7pm/test_pipeline/proj_result/result/chipseq/bowtie2/"$i".sam
done
