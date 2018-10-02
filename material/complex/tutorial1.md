---
title: "Using DESeq2 in R"
author: Daniel Neves and Daniel Sobral
date: "October 10, 2018"
output: 
  html_document: 
    keep_md: yes
---

This document demonstrates how to use *DESeq2* and *edgeR* in the *R environment* to perform a differential expression analysis using the the Trapnell datasets as an example. We will first need to tell R what samples are going to be analysed, then run the *DESeq2* pipeline plot the results of the analysis.

Finally, we will repeat the analysis using another common differential expression package called *edgeR*. 

# Setting up the environment

First we need to make sure that R is running on the same directory where we placed the counts files (the files called trapnell_counts_C1_R1.tab, trapnell_counts_C1_R2.tab, etc...). To do this either type `setwd("path/to/directory")` in the R console, or navigate to the counts directory using the *Files* panel in RStudio and select "Set As Working Directory".

# Setting up the count data and metadata

In this example, instead of loading the sample counts ourselves, we are going to let *DESeq2* handle that for us. For this, we just need to tell *DESeq2* what files correspond to each sample. We start by setting variables to hold the list of samples we are going to analyze. We create a list of sample names, a list of sample files (where the counts are), and a list of experimental conditions, telling which samples correspond to each condition. Type the following lines in the R console and press Enter.


```r
sampleNames <- c("trapnell_counts_C1_R1", "trapnell_counts_C1_R2", "trapnell_counts_C1_R3", "trapnell_counts_C2_R1", "trapnell_counts_C2_R2", "trapnell_counts_C2_R3")

sampleFiles <- c("trapnell_counts_C1_R1.tab", "trapnell_counts_C1_R2.tab", "trapnell_counts_C1_R3.tab", "trapnell_counts_C2_R1.tab", "trapnell_counts_C2_R2.tab", "trapnell_counts_C2_R3.tab")

sampleConditions <- c("C1", "C1", "C1", "C2", "C2", "C2")
```

We can confirm the values in these variables by simply typing a variable name in the R console and pressing Enter.


```r
sampleNames
```

```
## [1] "trapnell_counts_C1_R1" "trapnell_counts_C1_R2" "trapnell_counts_C1_R3"
## [4] "trapnell_counts_C2_R1" "trapnell_counts_C2_R2" "trapnell_counts_C2_R3"
```

```r
sampleFiles
```

```
## [1] "trapnell_counts_C1_R1.tab" "trapnell_counts_C1_R2.tab"
## [3] "trapnell_counts_C1_R3.tab" "trapnell_counts_C2_R1.tab"
## [5] "trapnell_counts_C2_R2.tab" "trapnell_counts_C2_R3.tab"
```

```r
sampleConditions
```

```
## [1] "C1" "C1" "C1" "C2" "C2" "C2"
```

For convenience, we place this information in a table variable that we call `sampleTable`.


```r
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)

sampleTable
```

```
##              sampleName                  fileName condition
## 1 trapnell_counts_C1_R1 trapnell_counts_C1_R1.tab        C1
## 2 trapnell_counts_C1_R2 trapnell_counts_C1_R2.tab        C1
## 3 trapnell_counts_C1_R3 trapnell_counts_C1_R3.tab        C1
## 4 trapnell_counts_C2_R1 trapnell_counts_C2_R1.tab        C2
## 5 trapnell_counts_C2_R2 trapnell_counts_C2_R2.tab        C2
## 6 trapnell_counts_C2_R3 trapnell_counts_C2_R3.tab        C2
```

# Running a differential expression test with DESeq2

With the sample table prepared, we are ready to run **DESeq2**. First need to import it into the R environment. This is done with the `library` command.


```r
library("DESeq2")
```

Then, we prepare a special structure to tell *DESeq2* what samples we are going to analyse (our sample table), and what comparison we are goind to make. Here we use the `condition` column (C1 or C2) as the experimental variable.


```r
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design= ~ condition)
```

We can run the whole DESeq2 pipeline with a single command. This will perform normalization of the raw counts, estimate variances, and perform the differential expression tests.


```r
ddsHTSeq <- DESeq(ddsHTSeq)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

We can then extract the results in the form of a table using the `results` function. The `head` function will print the first lines of this table on the console.


```r
resHTSeq <- results(ddsHTSeq)

head(resHTSeq)
```

```
## log2 fold change (MLE): condition C2 vs C1 
## Wald test p-value: condition C2 vs C1 
## DataFrame with 6 rows and 6 columns
##               baseMean log2FoldChange      lfcSE       stat     pvalue
##              <numeric>      <numeric>  <numeric>  <numeric>  <numeric>
## FBgn0000003    0.00000             NA         NA         NA         NA
## FBgn0000008  578.61667    -0.01363018 0.09778363 -0.1393912 0.88914103
## FBgn0000014   98.25848     0.26279172 0.22066004  1.1909348 0.23367917
## FBgn0000015   65.03895     0.12821306 0.26373025  0.4861523 0.62685920
## FBgn0000017 2539.89350     0.11317799 0.05904858  1.9166927 0.05527698
## FBgn0000018  327.26607    -0.06538546 0.12276560 -0.5326041 0.59430768
##                  padj
##             <numeric>
## FBgn0000003        NA
## FBgn0000008 0.9992327
## FBgn0000014 0.9515670
## FBgn0000015        NA
## FBgn0000017 0.6446372
## FBgn0000018 0.9992327
```

We can ask how many genes are differentially expressed (using a cutoff of 0.05) with this command.


```r
table(resHTSeq$padj < 0.05)
```

```
## 
## FALSE  TRUE 
##  6322   267
```


Finally, we sort this table by p-value (smaller p-values on top), and save it to a file so that we can later import it into Excel. 


```r
orderedRes <- resHTSeq[ order(resHTSeq$padj), ]

write.csv(as.data.frame(orderedRes), file="trapnell_C1_VS_C2.DESeq2.csv")
```

We can also retrieve and save a table of normalized counts.


```r
normCounts <- counts(ddsHTSeq, normalized = TRUE)

head(normCounts)
```

```
##             trapnell_counts_C1_R1 trapnell_counts_C1_R2
## FBgn0000003               0.00000               0.00000
## FBgn0000008             601.61191             601.95507
## FBgn0000014              88.44375              78.64252
## FBgn0000015              70.94939              65.04998
## FBgn0000017            2593.05424            2334.03223
## FBgn0000018             318.78628             333.01708
##             trapnell_counts_C1_R3 trapnell_counts_C2_R1
## FBgn0000003                0.0000               0.00000
## FBgn0000008              540.4924             544.48136
## FBgn0000014              100.9178              89.20881
## FBgn0000015               50.4589              56.39637
## FBgn0000017             2393.8864            2614.74099
## FBgn0000018              352.2419             311.71814
##             trapnell_counts_C2_R2 trapnell_counts_C2_R3
## FBgn0000003               0.00000               0.00000
## FBgn0000008             622.53341             560.62587
## FBgn0000014             127.98795             104.35007
## FBgn0000015              72.69716              74.68191
## FBgn0000017            2672.38848            2631.25864
## FBgn0000018             294.88425             352.94877
```

```r
write.csv(as.data.frame(orderedRes), file="trapnell_normCounts.DESeq2.csv")
```

# Visualizing results

*DESeq2* provides several functions to visualize the results, while additional plots can be made using the extensive R graphics cappabilities. Visualization can help to better understand the results, and catch potential problems in the data and analysis.

We can plot the *DESeq2* dispersion re-estimation procedure by typing:


```r
plotDispEsts(ddsHTSeq)
```

![](tutorial1_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

As a sanity check, we can inspect the distribution of p-values using the `hist` function.


```r
hist(resHTSeq$pvalue, breaks=0:50/50, xlab="p value", main="Histogram of nominal p values")
```

![](tutorial1_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

Two common visualizations for differential expression analyses are the *MA-plot*, that displays the relationship between a genes' mean expression and its fold-change between experimental conditions, and the *Volcano plot*, that displays the relationship between fold-change and evidence of differential expression (represented as -log p-value).

To display an *MA-plot* type the following in the R console. 


```r
plotMA(resHTSeq)
```

![](tutorial1_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

*DESeq2* doesn't provide a function to display a *Volcano plot*, but we can create one using R's base plot functions. In red we highlight genes differentially expressed with Padj < 0.05.


```r
plot(resHTSeq$log2FoldChange, -log10(resHTSeq$pvalue), xlab="log2 Fold-change", ylab="-log P-value", pch=20, cex=0.5)
points(resHTSeq$log2FoldChange[ resHTSeq$padj<0.05 ], -log10(resHTSeq$pvalue[ resHTSeq$padj<0.05 ]), col="red", pch=20, cex=0.5)
abline(v=0, h=-log10(0.05), lty="dashed", col="grey")
```

![](tutorial1_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

*DESeq2* provides a function to make a Principal Component Analysis (PCA) of the count data. The *DESeq2* [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations) recommend using transformed counts as input to the PCA routines, as these transformations remove the dependence of the sample-to-sample variance on the genes' mean expression. We do this with the `varianceStabilizingTransformation` function.


```r
vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)

plotPCA(vsd)
```

![](tutorial1_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

Another common visualization of high-throughput datasets is a clustered heatmap of sample-to-sample distances (or correlations). This visualization groups togheter the samples that are more similar to each other. As expected, we see that the Trapnell samples group according to condition (C1 or C2). Here we use the transformed counts defined above.


```r
dists <- dist(t(assay(vsd)))

# headmap of distances
heatmap(as.matrix(dists), main="Clustering of euclidean distances", scale="none")
```

![](tutorial1_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Here we plot the relative expression of all differentially expressed genes in the 6 samples. This figure is useful to visualize the differences in expression between samples. 


```r
library(gplots)
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:IRanges':
## 
##     space
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     space
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
diffgenes <- rownames(resHTSeq)[ which(resHTSeq$padj < 0.05) ]
diffcounts <- normCounts[ diffgenes, ]

heatmap.2(diffcounts, 
          labRow = "", 
          trace = "none", density.info = "none",
          scale = "row",
          distfun = function(x) as.dist(1 - cor(t(x))))
```

![](tutorial1_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

The following commands are used to plot a heatmap of the 20 most differentially expressed genes. For this, we use the ordered results table to determine which genes are most differentially expressed, and then plot the values from the normalized counts table (transformed to log10).


```r
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
```

![](tutorial1_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

# Differential expression with edgeR

Another commonly used package for differential expression analysis is *edgeR*. Here we repeat the analysis using *edgeR* and compare the results with those of *DESeq2*.

To initialize *edgeR* we need to provide a table containing the count data for all samples in individual columns. For this we load individual samples into R, and merge them into a single table with the following commands.


```r
sampleFiles <- c("trapnell_counts_C1_R1.tab", "trapnell_counts_C1_R2.tab", "trapnell_counts_C1_R3.tab", "trapnell_counts_C2_R1.tab", "trapnell_counts_C2_R2.tab", "trapnell_counts_C2_R3.tab")

tabs <- lapply(sampleFiles, function(x) read.table(x, col.names = c("Gene", x)))
countdata <- Reduce(f = function(x, y) merge(x, y, by="Gene"), x = tabs)

head(countdata)
```

```
##          Gene trapnell_counts_C1_R1.tab trapnell_counts_C1_R2.tab
## 1 FBgn0000003                         0                         0
## 2 FBgn0000008                       619                       620
## 3 FBgn0000014                        91                        81
## 4 FBgn0000015                        73                        67
## 5 FBgn0000017                      2668                      2404
## 6 FBgn0000018                       328                       343
##   trapnell_counts_C1_R3.tab trapnell_counts_C2_R1.tab
## 1                         0                         0
## 2                       557                       531
## 3                       104                        87
## 4                        52                        55
## 5                      2467                      2550
## 6                       363                       304
##   trapnell_counts_C2_R2.tab trapnell_counts_C2_R3.tab
## 1                         0                         0
## 2                       608                       548
## 3                       125                       102
## 4                        71                        73
## 5                      2610                      2572
## 6                       288                       345
```

```r
rownames(countdata) <- as.character(countdata$Gene)
countdata$Gene<-NULL
```

We load *edgeR* into R. 


```r
library(edgeR)
```

```
## Loading required package: limma
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:DESeq2':
## 
##     plotMA
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

We initialize *edgeR* by providing the counts table, gene names, and experimental conditions.


```r
mygroups <- c("C1","C1","C1","C2","C2","C2")

y <- DGEList(counts=countdata, genes=rownames(countdata), group = mygroups)
```

To run *edgeR* in classic mode, we need to perform 3 steps: calculate normalization factors, estimade dispersions, and finally perform the exact test for differential expression.


```r
y <- calcNormFactors(y)
y <- estimateDisp(y)
```

```
## Design matrix not provided. Switch to the classic mode.
```

```r
et <- exactTest(y)
```

We extract the results using the function `topTags`.


```r
result_edgeR <- as.data.frame(topTags(et, n=nrow(countdata)))

table(result_edgeR$FDR < 0.05)
```

```
## 
## FALSE  TRUE 
## 17224   266
```

```r
plot(result_edgeR$logFC, -log10(result_edgeR$FDR), col=ifelse(result_edgeR$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
```

![](tutorial1_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
hist(result_edgeR$PValue, breaks=20, xlab="P-Value", col="royalblue", ylab="Frequency", main="P-value distribution")
```

![](tutorial1_files/figure-html/unnamed-chunk-23-2.png)<!-- -->

# Comparison between the two methods

To compare the result of the two methods, we first merge both results tables. 


```r
comp_table <- merge(as.data.frame(resHTSeq), result_edgeR, by="row.names")

head(comp_table)
```

```
##     Row.names   baseMean log2FoldChange      lfcSE       stat     pvalue
## 1 FBgn0000003    0.00000             NA         NA         NA         NA
## 2 FBgn0000008  578.61667    -0.01363018 0.09778363 -0.1393912 0.88914103
## 3 FBgn0000014   98.25848     0.26279172 0.22066004  1.1909348 0.23367917
## 4 FBgn0000015   65.03895     0.12821306 0.26373025  0.4861523 0.62685920
## 5 FBgn0000017 2539.89350     0.11317799 0.05904858  1.9166927 0.05527698
## 6 FBgn0000018  327.26607    -0.06538546 0.12276560 -0.5326041 0.59430768
##        padj       genes       logFC    logCPM     PValue FDR
## 1        NA FBgn0000003  0.00000000 -2.457703 1.00000000   1
## 2 0.9992327 FBgn0000008 -0.01016611  5.728385 0.92504901   1
## 3 0.9515670 FBgn0000014  0.26626896  3.192089 0.13105356   1
## 4        NA FBgn0000015  0.13150517  2.612451 0.51949013   1
## 5 0.6446372 FBgn0000017  0.11677540  7.858463 0.06306125   1
## 6 0.9992327 FBgn0000018 -0.06179114  4.910233 0.60647902   1
```

We can then ask for a table comparing differentially expressed genes. Only 8 genes were classified differently in the two tests.


```r
table("DESeq2" = comp_table$padj < 0.05, "edgeR" = comp_table$FDR < 0.05)
```

```
##        edgeR
## DESeq2  FALSE TRUE
##   FALSE  6320    2
##   TRUE      6  261
```

# Independent filtering

Notice from the edgeR p-value distribution that we seem to have alot of genes with FDR close to 1. This may indicate the presence of many genes not being expressed, or with very low expression. 

It is often recommended to remove these genes from the analysis. Doing so will reduce the number of statistical tests we are making, and has an impact on the p-value adjustments.

We can use the command below to remove from the table of counts all genes that are not expressed.


```r
w <- which(rowSums(countdata) > 0)
countdata <- countdata[ w, ]
```

Run the code below to repeat the edgeR analysis with this counts table.


```r
y <- DGEList(counts=countdata, genes=rownames(countdata), group = mygroups)
y <- calcNormFactors(y)
y <- estimateDisp(y)
```

```
## Design matrix not provided. Switch to the classic mode.
```

```r
et <- exactTest(y)

result_edgeR_2 <- as.data.frame(topTags(et, n=nrow(countdata)))

table(result_edgeR_2$FDR < 0.05)
```

```
## 
## FALSE  TRUE 
##  9937   271
```

```r
hist(result_edgeR_2$PValue, breaks=20, xlab="P-Value", col="royalblue", ylab="Frequency", main="P-value distribution")
```

![](tutorial1_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


