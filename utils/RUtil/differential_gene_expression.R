# Title     : TODO
# Objective : TODO
# Created by: Jason
# Created on: 2023/5/17

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

countData <- as.matrix(read.csv("data/gene_count_matrix.csv",row.names="gene_id"))

# filter genes that have low expression
countData <- countData[rowMeans(countData)>1,]

# divide in groups: two Veh and two DHT
condition <- factor(c(rep("Replicate_DHT",2),rep("Replicate_Veh",2)))

# 
colData <- data.frame(row.names=colnames(countData), condition)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# check
head(dds)

# Run the differential expression analysis
dds <- DESeq(dds)

# normalization
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
#将结果用result()函数来获取
res <- results(dds1)
summary(dds1)

# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)

#可以直接用DESeq2的plotCounts
dds <- makeExampleDESeqDataSet()
plotCounts(dds1, gene = "HNRNPUL2|HNRNPUL2",intgroup = "condition")   # 指定某个基因

#可用ggplot对图片样式进行修改，并用ggrepel进行标注
d=data.frame(t(subset(countData,rownames(countData)=="HNRNPUL2|HNRNPUL2")))
ggplot(d, aes(x = condition, y = AT4G38770, color = condition))+
  geom_point(position=position_jitter(w=0.2,h=0))+
  geom_text_repel(aes(label=rownames(d)))+
  theme_bw()+
  ggtitle("HNRNPUL2|HNRNPUL2")+
  theme(plot.title=element_text(hjust=0.5))


df <- countData[intersect(rownames(countData),rownames(res1_total)),]    
# 在原表达矩阵中找到差异表达基因
df2<- as.matrix(df)                                                 
pheatmap(df2,
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=T,
         height=10,  
         scale = "row",
         frontsize = 10,
         angle_col=45, 
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
         clustering_method = 'single',
) 


genes<- res1
# 根据上调、下调、不变为基因添加颜色信息
genes$color <- ifelse(genes$padj<0.05 & abs(genes$log2FoldChange)>= 1,ifelse(genes$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p <- ggplot(
  # 指定数据、映射、颜色
  genes, aes(log2FoldChange, -log10(padj), col = color)) +  
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  # 辅助线
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  # 图例
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)+
          # 注释
          geom_text_repel(
            data = subset(genes, padj < 1e-100 & abs(genes$log2FoldChange) >= 10),
            aes(label = rownames(genes)),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines"))
  )
p
