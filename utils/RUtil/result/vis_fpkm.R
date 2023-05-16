# Title     : TODO
# Objective : TODO
# Created by: Jason
# Created on: 2023/5/15

# read data
data_fpkm <- read.csv("data/gene_fpkm_matrix.csv", row.names = 1)

# Assuming data_fpkm is your FPKM matrix
# Convert the matrix to a numeric vector
fpkm_values <- as.numeric(as.matrix(data_fpkm))

# Plot a histogram with a logarithmic scale
hist(fpkm_values, breaks=100, main="Histogram of log-transformed FPKM values", xlab="log(FPKM)", log="x")

# Alternatively, plot a density plot
plot(density(fpkm_values, na.rm=TRUE), main="Density plot of log-transformed FPKM values", xlab="log(FPKM)")

# Set a threshold for low expression
threshold <- 1

# Calculate the row means of the FPKM matrix
row_means <- rowMeans(data_fpkm, na.rm=TRUE)

# Create a logical vector indicating whether each row mean is above the threshold
keep <- row_means > threshold

# Subset the FPKM matrix to keep only the rows with above-threshold means
data_fpkm_filtered <- data_fpkm[keep, ]
