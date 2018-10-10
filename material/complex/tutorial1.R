## ---- eval=FALSE, echo=FALSE---------------------------------------------
## knitr::purl("tutorial1.Rmd")

## ------------------------------------------------------------------------
sampleNames <- c("trapnell_counts_C1_R1", "trapnell_counts_C1_R2", "trapnell_counts_C1_R3", "trapnell_counts_C2_R1", "trapnell_counts_C2_R2", "trapnell_counts_C2_R3")

sampleFiles <- c("trapnell_counts_C1_R1.tab", "trapnell_counts_C1_R2.tab", "trapnell_counts_C1_R3.tab", "trapnell_counts_C2_R1.tab", "trapnell_counts_C2_R2.tab", "trapnell_counts_C2_R3.tab")

sampleConditions <- c("C1", "C1", "C1", "C2", "C2", "C2")

## ------------------------------------------------------------------------
sampleNames
sampleFiles
sampleConditions

## ------------------------------------------------------------------------
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)

sampleTable

## ---- message=FALSE------------------------------------------------------
library("DESeq2")

## ------------------------------------------------------------------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design= ~ condition)

## ------------------------------------------------------------------------
ddsHTSeq <- DESeq(ddsHTSeq)

## ------------------------------------------------------------------------
resHTSeq <- results(ddsHTSeq)

head(resHTSeq)

## ------------------------------------------------------------------------
table(resHTSeq$padj < 0.05)

## ------------------------------------------------------------------------
table(resHTSeq$padj < 0.01)
table(resHTSeq$pvalue < 0.01)

## ------------------------------------------------------------------------
orderedRes <- resHTSeq[ order(resHTSeq$padj), ]

write.csv(as.data.frame(orderedRes), file="trapnell_C1_VS_C2.DESeq2.csv")

## ------------------------------------------------------------------------
normCounts <- counts(ddsHTSeq, normalized = TRUE)

head(normCounts)

write.csv(as.data.frame(orderedRes), file="trapnell_normCounts.DESeq2.csv")

## ------------------------------------------------------------------------
merged.results <- merge(normCounts, orderedRes, by="row.names")

head(merged.results)

## ------------------------------------------------------------------------
merged.results <- merged.results[ order(merged.results$padj), ]
head(merged.results)

## ------------------------------------------------------------------------
plotDispEsts(ddsHTSeq)

## ------------------------------------------------------------------------
hist(resHTSeq$pvalue, breaks=0:50/50, xlab="p value", main="Histogram of nominal p values")

## ------------------------------------------------------------------------
plotMA(resHTSeq)

## ------------------------------------------------------------------------
resultsNames(ddsHTSeq)
resHTSeqShrunk <- lfcShrink(ddsHTSeq, coef=2)
plotMA(resHTSeqShrunk)

## ------------------------------------------------------------------------
highlight <- which(resHTSeq$padj < 0.05)

plot(resHTSeq$log2FoldChange, -log10(resHTSeq$pvalue), xlab="log2 Fold-change", ylab="-log P-adjusted", pch=20, cex=0.5)
points(resHTSeq$log2FoldChange[ highlight ], -log10(resHTSeq$pvalue[ highlight ]), col="red", pch=20, cex=0.5)
abline(v=0, h=-log10(0.05), lty="dashed", col="grey")

## ------------------------------------------------------------------------
highlight <- which(resHTSeqShrunk$padj < 0.01)

plot(resHTSeqShrunk$log2FoldChange, -log10(resHTSeqShrunk$pvalue), xlab="shrunken log2 Fold-change", ylab="-log P-adjusted", pch=20, cex=0.5)
points(resHTSeqShrunk$log2FoldChange[ highlight ], -log10(resHTSeqShrunk$pvalue[ highlight ]), col="green", pch=20, cex=0.5)
abline(v=0, h=-log10(0.01), lty="dashed", col="grey")

## ------------------------------------------------------------------------
transformed.vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=TRUE)

plotPCA(transformed.vsd)

## ------------------------------------------------------------------------
dists <- as.matrix(dist(t(normCounts)))
heatmap(dists, main="Clustering of sample-to-sample distances", scale="none")

## ------------------------------------------------------------------------
log10_rawCounts <- log10(counts(ddsHTSeq) + 1)
  
dists <- 1 - cor(log10_rawCounts, method="pearson")
heatmap(dists, main="Clustering of sample-to-sample pearson correlations", scale="none")

## ---- fig.height=8, fig.width=5------------------------------------------
library(gplots)

diffgenes <- rownames(resHTSeq)[ which(resHTSeq$padj < 0.05) ]
diffcounts <- normCounts[ diffgenes, ]

heatmap.2(diffcounts, 
          labRow = "", 
          trace = "none", density.info = "none",
          scale = "row",
          distfun = function(x) as.dist(1 - cor(t(x))))

## ------------------------------------------------------------------------
library(pheatmap)

# select the 20 most differentially expressed genes
select <- row.names(orderedRes[1:20, ])

# transform the counts to log10
log10_normCounts <- log10(normCounts + 1)

# get the values for the selected genes
values <- log10_normCounts[ select, ]

pheatmap(values,
         scale = "none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_row = 8,
         annotation_names_col = FALSE,
         gaps_col = c(3,6),
         display_numbers = TRUE,
         number_format = "%.2f",         
         height=12,
         width=6)

## ------------------------------------------------------------------------
sessionInfo()

