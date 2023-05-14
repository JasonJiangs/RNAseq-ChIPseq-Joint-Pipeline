# Title     : DESeq2
# Objective : Use DESeq2 to do differential expression analysis
# Created by: Jason
# Created on: 2023/5/14

# Usage: Rscript /path/deseq.R /path/transcript_count_matrix.csv /path/results.csv

library(DESeq2)

cmd_list <- commandArgs(trailingOnly = TRUE)

# "transcript_count_matrix.csv"
transcript_count_matrix_path <- cmd_list[0]
# "results.csv"
result_path <- cmd_list[1]

database <- as.matrix(read.csv(transcript_count_matrix_path, row.names="transcript_id"))
condition <- factor(c("control","control","KD","KD"))

coldata <- data.frame(row.names = colnames(database), condition)
countData <- countData[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
resordered <- res[order(res$padj),]

summary(res)

write.csv(as.data.frame(resordered), file=result_path)

