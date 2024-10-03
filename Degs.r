####unique_GC_matrix_for_new.csv obtained after gene quantification"
library(DESeq2)
setwd("/home/rohan/Documents/Sakshi/mh_rna_last/output/new_deseqanalysis")

# Loading the data in R
countData <- read.table("unique_GC_matrix_for_new.csv", header=TRUE, sep=",")
head(countData)
rownames(countData) <- countData[,1] # Rename the row names with gene names
countData <- countData[,-1] # Remove the gene column

# Loading the sample information
colData <- read.csv("/home/rohan/Downloads/phenoinfo.csv", sep=",", row.names=1)
head(colData)
# Ensure the order of samples in colData matches the columns in countData
countData <- countData[,rownames(colData)]

# Print column names of countData
colnames(countData)

# Print row names of colData
rownames(colData)

# Check if all row names in colData are in the columns of countData
all(rownames(colData) %in% colnames(countData))

# Check if all columns in countData are in row names of colData
all(colnames(countData) %in% rownames(colData))



# Check if row names of colData match the column names of countData
if (!all(rownames(colData) %in% colnames(countData))) {
   stop("Mismatch between sample names in colData and countData. Check the sample names in both files.")
}

if (!all(rownames(colData) == colnames(countData))) {
   stop("Order of samples in colData does not match the columns in countData. Check the order of samples.")
}

# Convert 'condition' and 'batch' columns to factors
colData$condition <- as.factor(colData$condition)
colData$batch <- as.factor(colData$batch)

# Check and ensure colData contains 'batch'
if (!"batch" %in% colnames(colData)) {
   stop("The 'batch' column is missing in colData.")
}


# Filter out genes with a total count of less than 50 across all samples
countData <- countData[rowSums(countData) >= 50, ]
cat("Number of genes after filtering low counts:", nrow(countData), "\n")


# Create DESeq2 dataset object with batch correction
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ batch + condition)

# Normalize counts
dds1 <- estimateSizeFactors(dds)
normalized_counts <- counts(dds1, normalized=TRUE)
head(normalized_counts)

# Save normalized counts to a CSV file
write.csv(normalized_counts, file="DEseq2_normalized_counts.csv", quote=F)


# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results for the condition of interest (Test vs. Control)
res <- results(dds, contrast = c("condition", "Test", "Control"), alpha = 0.05)
summary(res)

# Save the DESeq2 results
write.table(res, file="NEWanalysis_deseq_result.csv", sep="\t", quote=F, col.names=NA)

upregulated_genes <- res[which(res$log2FoldChange > 1 & res$padj < 0.05), ]
downregulated_genes <- res[which(res$log2FoldChange < -1 & res$padj < 0.05), ]
upregulated_genes_ordered <- upregulated_genes[order(upregulated_genes$pvalue), ]
downregulated_genes_ordered <- downregulated_genes[order(downregulated_genes$pvalue), ]



# Summary of upregulated genes
cat("Summary of Upregulated Genes:\n")
cat("Total number of upregulated genes:", nrow(upregulated_genes), "\n")
summary(upregulated_genes)
summary(downregulated_genes)
# View first few rows of the ordered upregulated genes
head(upregulated_genes_ordered)



# Example of examining the distribution of counts
boxplot(rowSums(countData), main="Total Counts per Gene", ylab="Counts")


# Perform PCA on the normalized counts data (transposed)
pca_res <- prcomp(t(normalized_counts), scale. = TRUE)

# Plot the first two principal components
plot(pca_res$x[,1], pca_res$x[,2],
     xlab="PC1", ylab="PC2",
     main="PCA of Normalized Counts Data",
     col="violet", pch=19)
abline(h=0, v=0, lty=2, col="gray")

# Add text labels for each point in the plot
text(pca_res$x[,1], pca_res$x[,2], labels=rownames(normalized_counts), pos=4, cex=0.7)


# Load the normalized count data
normalized_counts <- read.table("Batchanalysis_norm_count.csv", header=TRUE, sep="\t", row.names=1)

# Check the first few rows of normalized counts to ensure correct loading
head(normalized_counts)

# Load the sample information
colData <- read.csv("/home/rohan/Downloads/phenoinfo.csv", sep=",", row.names=1)
head(colData)
# Ensure the order of samples in colData matches the columns in countData
normalized_counts <- normalized_counts[, rownames(colData)]

# Check column names of normalized_counts
colnames(normalized_counts)

# Check row names of colData
rownames(colData)


# Convert the data to a matrix format for plotting
normalized_counts_matrix <- as.matrix(normalized_counts)

# Generate a boxplot for the normalized counts
boxplot(normalized_counts_matrix,
        las=2,  # Rotate axis labels for better readability
        main="Boxplot of Normalized Counts by Sample",
        ylab="Normalized Count",
        col="pink",
        outline=FALSE)  # Remove outliers

# Add a horizontal line at the median value of normalized counts
median_normalized_counts <- median(normalized_counts_matrix)
abline(h=median_normalized_counts, col="black", lty=2)  # Black dashed line



write.table(as.data.frame(upregulated_genes_ordered),"/home/bioinfo/Upregulated_genes_Deseq.txt", sep = "\t")
