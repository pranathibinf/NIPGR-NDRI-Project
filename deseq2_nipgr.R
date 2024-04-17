# Load required libraries
library(DESeq2)
library(openxlsx)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)

# Set the working directory to the specified folder
setwd("path_to_file/gene_count_matrix.csv")

# Load gene count matrix from a CSV file
countData <- read.table("C:/Users/Lenovo/OneDrive/Documents/gene_count_matrix.csv", header = TRUE, sep = ",", row.names = 1)
head(countData)

# Load phenotypic data from a CSV file
colData <- read.table("phenodata2.csv", header = TRUE, sep = ",")
head(colData)

# Create an example dataset with 4 conditions
dds <- makeExampleDESeqDataSet(m=4)
# Specify the model design formula as an intercept only (no conditions)
design(dds) <- formula(~ 1)

# Create a DESeqDataSet object using count data and phenotypic data with a specific model design
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~stage_condition)
# Filter to keep only rows with nonzero counts
dds <- dds[rowSums(counts(dds)) > 0]

# Run the DESeq pipeline for differential expression analysis
dds <- DESeq(dds)

# Retrieve the results of the analysis and sort by adjusted p-value
res <- results(dds)
head(res)
summary(res)
res <- res[order(res$padj),]

# Subset the results for significant entries based on adjusted p-value and log2 fold change
res_sig <- subset(res, padj<0.05)
res_lfc <- subset(res_sig, abs(log2FoldChange) > 2) 
head(res_lfc)

# Save the significant results to an Excel file
write.xlsx(as.data.frame(res_sig), file="ressig.xlsx",rowNames=T,asTable = TRUE)

# Apply a regularized log transformation to the count data
rld <- rlog(dds,blind=F)

# Identify the top 100 genes with the highest variance
topVarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T),100)
matrix <- assay(rld)[ topVarianceGenes, ]
matrix <- matrix - rowMeans(matrix)

# Plot MA plot of results
plotMA(res)

# Perform PCA analysis
dds <- estimateSizeFactors(dds)
# Log transform the normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
# Plot PCA
plotPCA(DESeqTransform(se))

# Generate a volcano plot for visualizing differential expression
EnhancedVolcano(res,lab = rownames(res),x = 'log2FoldChange',y = 'pvalue',labSize = 0, pCutoff = 0.05,FCcutoff = 1,xlim=c(-30,12))
plot(EnhancedVolcano)

# Generate a heatmap of sample distances
variance_dds <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(variance_dds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(variance_dds$Condition, variance_dds$sample, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

# Create heatmaps for the 20 most highly expressed genes
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[, c("sample", "stage_condition")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
