file_path <- "DAVID analysis/data/enhancer+promoter/GO MF.txt"

# 读取 tsv 文件
enrichment_results <- read.delim(file_path, header = TRUE, sep = "\t")

# 加载所需的包
library(ggplot2)
library(clusterProfiler)

# 设置富集分析结果数据
pathway_data <- enrichment_results[, c("Term", "PValue", "Count", "GeneRatio")]

# 根据P值排序通路数据
pathway_data <- pathway_data[order(pathway_data$PValue), ]

# 选择前20个通路数据
top20_pathways <- pathway_data[1:20, ]

# 创建因子变量，按照P值从小到大的顺序排列
top20_pathways$Term <- factor(top20_pathways$Term, levels = top20_pathways$Term[order(top20_pathways$PValue)])

# 创建气泡图
ggplot(top20_pathways, aes(x = GeneRatio, y = reorder(Term, -PValue), size = Count, color = PValue)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(1, 3)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "GeneRatio", y = "GO MF Enrichment Term", size = "Count", color = "PValue") +
  theme_minimal()

