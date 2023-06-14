# remove.packages("rlang")
# remove.packages("tidyverse")
# remove.packages("ggplot2")

# install.packages("tidyverse")
# tidyverse自动安装依赖的rlang、ggplot2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
BiocManager::install("edgeR")
BiocManager::install("ggplot2")

library(ggplot2)
library(edgeR)

# set command line arguments, the format is: Rscript xxx.R arg1 arg2 arg3
args <- commandArgs(trailingOnly = TRUE)

input_dir_path <- args[1]
output_dir_path <- args[2]

# loop all files in input_dir_path
for (file_name in list.files(input_dir_path)) {
  # get file path
  file_path <- paste(input_dir_path, file_name, sep = "/")
  # read file
  data <- read.csv(file_path, row.names = 1)


}

# create a folder 'result' in current directory if not exist
if (!file.exists("result")) {
  dir.create("result")
}

data <- read.csv("data/transcript_count_matrix.csv", row.names = 1)

# Normalization
data_norm <- cpm(data)

# Calculate correlation
correlation_matrix <- cor(data_norm, method = "spearman")
# save
write.csv(correlation_matrix, "result/correlation_matrix.csv")

# Visualization
heatmap:heatmap(correlation_matrix)

# advanced visualization
pheatmap::pheatmap(correlation_matrix)
# save
png("result/heatmap.png", width = 800, height = 800)

# with ggplot2
# Check the structure of your data
# print(str(data_norm))

# Check the column names of your data
# print(names(as.data.frame(data_norm)))
data_norm02 <- as.data.frame(data_norm)
# names(data_norm) <- c("LNCaP_DHT_RNA_rep1_batch", "LNCaP_DHT_RNA_rep2_batch",
#                       "LNCaP_Veh_RNA_rep1_batch", "LNCaP_Veh_RNA_rep2_batch")

# Replace 'Replicate1' and 'Replicate2' with your actual column names
column1 <- "LNCaP_DHT_RNA_rep1_batch"
column2 <- "LNCaP_DHT_RNA_rep2_batch"

# Check if the columns exist
if (!(column1 %in% names(data_norm02)) | !(column2 %in% names(data_norm02))) {
  print(paste("Columns", column1, "and/or", column2, "not found in the data"))
} else {
  # Plot the data
  plot <- ggplot(data_norm02, aes_string(x=column1, y=column2)) + 
    geom_point() +
    theme_bw() +
    labs(x=paste("Expression", column1), 
         y=paste("Expression", column2), 
         title=paste("Scatter plot of", column1, "vs", column2)) +
    geom_smooth(method=lm , color="red", se=FALSE)
  
    print(plot)
  # save
    ggsave("result/scatter_plot01.png", plot, width = 8, height = 8)
}

column1 <- "LNCaP_Veh_RNA_rep1_batch"
column2 <- "LNCaP_Veh_RNA_rep2_batch"

# Check if the columns exist
if (!(column1 %in% names(data_norm02)) | !(column2 %in% names(data_norm02))) {
  print(paste("Columns", column1, "and/or", column2, "not found in the data"))
} else {
  # Plot the data
  plot <- ggplot(data_norm02, aes_string(x=column1, y=column2)) + 
    geom_point() +
    theme_bw() +
    labs(x=paste("Expression", column1), 
         y=paste("Expression", column2), 
         title=paste("Scatter plot of", column1, "vs", column2)) +
    geom_smooth(method=lm , color="red", se=FALSE)
  
    print(plot)
  # save
    ggsave("result/scatter_plot02.png", plot, width = 8, height = 8)
}



