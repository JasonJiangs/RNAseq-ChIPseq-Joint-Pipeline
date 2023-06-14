#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# Load necessary libraries
library(DESeq2)

# Load the read counts data
counts <- read.csv("data/gene_count_matrix.csv", row.names = 1)

# Load the sample data
coldata <- read.csv("data/coldata.csv", row.names = 1)

# Convert the data to DESeq2 objects
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Run the differential expression analysis
dds <- DESeq(dds)

# Extract the results
res <- results(dds)

# Print the results
print(res)
#save
write.csv(as.data.frame(res), file = "result/deseq2/deseq2_results.csv")

# List significantly differentially expressed genes with |log2 fold change| > 1
resSig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.01),]
print(resSig)
#save
write.csv(as.data.frame(resSig), file = "result/deseq2/deseq2_results_sig_001.csv")

# Load necessary library
library(ggplot2)

# convert res to data.frame
res <- as.data.frame(res)

# Add a significance column to the results dataframe
res$sig <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")

# Create a volcano plot with colored points
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot", color = "Significance") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "firebrick") # line for p-value = 0.05

#save
ggsave("result/deseq2/deseq2_volcano_plot.png", width = 10, height = 8, units = "in")

# Load necessary library
library(pheatmap)

# Normalize the count data
rld <- rlog(dds)

# Extract the normalized counts
norm_counts <- assay(rld)

# Subset the normalized counts to just the differentially expressed genes
de_counts <- norm_counts[rownames(resSig),]

# Generate the heatmap
pheatmap(de_counts, annotation_col = coldata, scale = "row")

#save
png("result/deseq2/deseq2_heatmap.png", width = 800, height = 800)

# GeneOntology Enrichment Analysis
