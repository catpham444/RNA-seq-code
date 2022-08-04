# script to perform differential gene expression analysis using DESeq2 package

# load libraries
library(DESeq2)
library(tidyverse)

# Step 1: preparing count data ----------------

# read in counts data
x <- read.csv('counts_data.csv')
counts_data <- x %>% remove_rownames %>% column_to_rownames(var="ï..")

head(counts_data)


# read in sample info
y <- read.csv('sample_info.csv')
colData<- y %>% remove_rownames %>% column_to_rownames(var="X")
head(colData)

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res
write.csv(as.data.frame(res[order(res$padj),] ), file="result.csv")


# Explore Results ----------------

summary(res)

table(res$padj < 0.05 & (res$log2FoldChange>1 | res$log2FoldChange < -1)	)

write.table(subset(res, padj < 0.05& (res$log2FoldChange>1 | res$log2FoldChange < -1)),
            file = "DESeq2result.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)


DGE.results.shrnk <- lfcShrink(dds, coef = 2, type = "apeglm")
# MA plot
plotMA(DGE.results.shrnk)

#PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dexamethasone")

#volcano plot
library(EnhancedVolcano)
vp1 <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'padj', pCutoff = 0.05,
                       title = "treated / untreated")
print(vp1)

vp2 <- EnhancedVolcano(DGE.results.shrnk,
                       lab = rownames(DGE.results.shrnk),
                       x = 'log2FoldChange',
                       y = 'padj', pCutoff = 0.05,
                       title = "with logFC shrinkage")
library(patchwork)
vp1 + vp2

