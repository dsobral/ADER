![IGClogo](https://github.com/dsobral/ADER17S/raw/master/Logo_IGC_2014.png "IGC")
![Acknowledgments](https://github.com/dsobral/ADER/raw/master/material/images/acknowledgments.jpg "Acknowledgments")


# ADER18S #

## Oct 08 - Oct 12 2018

# Analysis of Differential Expression with RNAseq #

A hands-on training course at Instituto Gulbenkian de Ciência (5-days)

Official course page of the Gulbenkian Training Programme in Bioinformatics - GTPB

http://gtpb.igc.gulbenkian.pt/bicourses/2018/ADER18S/index.html



# Overview

This introductory course covers practical aspects of the analysis of differential gene expression by RNAseq, from planning the gathering of sequence data to the generation of tables of differentially expressed gene lists and visualization of results. For this edition of the course, we will also explore some specificities of single-cell RNAseq data analysis. Towards the end, we will cover some of the initial steps of secondary analysis, such as functional enrichment of the obtained gene lists. Participants will first start learning the concepts using small example datasets, and then will apply the learned concepts in the training room using real sized examples. At the end of the course, participants should be able to autonomously apply most of the learned methods to their own data.

# Target Audiences

Life Scientists who want to be able to use NGS data (RNAseq) to infer genes differentially expressed between different conditions. Computational researchers that wish to get acquainted with the concepts and methodologies used in RNAseq are also welcome. 


# Pre-requisites

Familiarity with elementary statistics and a few basics of scripting in R will be helpful.

Please have a look at the following resources and gauge your ability to use R in statitics at the basic level: [Coursera videos](http://blog.revolutionanalytics.com/2012/12/coursera-videos.html); [Introduction to r](http://bitesizebio.com/webinar/20600/beginners-introduction-to-r-statistical-software)

Basic Unix command line skills, such as being able to navigate in a directory tree and copy files. See, for example, ["Session 1" of the Software Carpentry training for a Unix introduction](http://bioinformatics-core-shared-training.github.io/shell-novice/). 


# Learning Objectives

Course participants will go through a series of experiences that utimately lead to create enhanced capabilities to:

1. List broad characteristics of NGS technologies and choose adequate sequencing for your biological question
2. Have a broad overview of the steps in the analysis of RNA-Seq differential expression experiments
3. Assess the general quality of the raw data from the sequencing facility
4. Do simple processing operations in the raw data to improve its quality
5. Generate alignments against a reference genome
6. Assess the general quality of the alignments and detect possible problems
7. Generate tables of counts using the alignment and a reference gene annotation
8. Generate lists of differentially expressed genes, at least for a simple pairwise comparison
9. Understand specificies of differential gene expression in the case of single-cell RNAseq
10. Perform simple functional enrichment analysis and understand the concepts behind them

For this, we are providing small example datasets and exercises that participants can use to learn. 


## Learning outcomes (LO) for each unit:

### LO 1 - Plan your experiment using NGS technologies:

####	LO 1.1 - [List possibilities and limitations of NGS sequencing technologies](material/L01_L02.md#LO1)
		What choices do you have when sending your samples to the sequencing facility

####	LO 1.2 - [Choose adequate sequencing for your biological question](material/L01_L02.md#LO1)
		How do the sequencing choices influence the kind of questions you can answer

### LO 2 - [List steps in the analysis of RNAseq differential expression experiments](material/L01_L02.md#LO2)
		What are the steps in RNAseq data analysis

### LO 3 - Assess the general quality of the raw data from the sequencing facility

#### 	LO 3.1 - [Interpret what are fastq files and what is their content](material/L03.md#LO3.1)
		What information is in fastq files, and how is it organized

#### 	LO 3.2 - [Use software like FastQC to process fastq files and produce QC reports](material/L03.md#LO3.2)
		Detect low quality bases in the QC reports
		Detect sequence bias and possible presence of adaptors and other contaminants

### LO 4 - Do simple processing operations in the raw data to improve its quality

#### 	LO 4.1 - [Use  trimmomatic to remove low quality bases from your reads](material/L04.md#LO4.1)
		Use trimmomatic to filter/trim low quality bases from your reads

#### 	LO 4.2 - [Use trimmomatic to remove adaptors and other artefactual sequences from your reads](material/L04.md#LO4.2)
		Remove adaptors (eg. illumina adaptors) and unwanted sequences (eg. polyA tails) from your reads
		Check results using FastQC on filtered data

### LO 5 - Generate alignments of processed reads against a reference genome

#### 	LO 5.1 - [What is a reference genome, versioning and where to obtain genomes](material/L05.md#LO5.1)
		Are genomes constant?
		Obtain genome fasta from Ensembl

#### 	LO 5.2 - [Alignment software: hisat2](material/L05.md#LO5.2)
		What are the requisites for using burrows-wheeler approaches?	
		Prepare a reference genome to use with hisat2

#### 	LO 5.3 - [Run an alignment: the SAM/BAM alignment format](material/L05.md#LO5.3)
		Run hisat2 in an example dataset
		What is the SAM/BAM format
	
### LO 6 - Assess the general quality of the alignments and detect possible problems

#### 	LO 6.1 - [What is a reference gene annotation, versioning and where to obtain](material/L06.md#LO6.1)
		What is the GFF/GTF format
		Obtain genome GTF from Ensembl
			
#### 	LO 6.2 - [Visualizing alignments in IGV for single genes](material/L06.md#LO6.2)

#### 	LO 6.3 - [Use Qualimap to assess quality of alignments](material/L06.md#LO6.3)
		Interpret general alignment statistics such as percentage of aligned reads
		Check the reports to assess RNA integrity and diversity

### LO 7 - Generate tables of counts

#### 	LO 7.1 - [The process of generating gene counts from genome aligments](material/L07.md#LO7.1)
		What parameters we need to consider when counting

#### 	LO 7.2 - [Use featurecounts to generate table of gene counts](material/L07.md#LO7.2)
		Interpret results from featurecounts

#### 	LO 7.3 - [Using Salmon to generate counts only with the transcriptome](material/L07.md#LO7.3)
		Interpret results from salmon

### LO 8 - Generate lists of differentially expressed genes, at least for a simple pairwise comparison

#### 	LO 8.1 - [Execute a pairwise differential expression analysis](material/L08.md#LO8.1)
		Use Galaxy to produce differentially expressed genes with DESeq2

#### 	LO 8.2 - [Interpretation and visualization of results](material/L08.md#LO8.2)
		PCA plots comparing all samples: detection of outliers, and batch effects
		Heatmaps and other plots

#### 	LO 8.3 - [Use more complex settings than simple pairwise comparisons](material/L08.md#LO8.3)
		Account for batch effects and paired data 
				
#### 	LO 8.4 - [Gain control over your analysis using R and Rstudio](material/L08.md#LO8.4)
		Use R in Rstudio to make a pairwise comparison using DESeq2 and edgeR
		Use edgeR to perform more complex analysis such as ANOVA-like all versus all comparisons

### LO 9 - Understand specificies of differential gene expression in single-cell RNAseq

#### 	LO 9.1 - Generate a count matrix for a single-cell RNAseq dataset (Chromium - 10x Genomics)
		Specificities of single-cell RNAseq, using the Chromium system as example
		Diferences in raw data preprocessing and counting

#### 	LO 9.2 - Generate groups of cells by clustering gene expression
		Quality checking and filtering of the count table
		Interpreting PCA plots and dimensionality reduction

#### 	LO 9.3 - Obtain marker genes for the different groups of cells
		Indentify genes that distinguish the different groups

### LO 10 - Perform simple functional enrichment analysis and understand the concepts involved
		
#### 	LO 10.1 - [How to extract meaning from a list of genes](material/L10.md#LO10.1)
		What are functional annotations, what types exist, and where to get them
       
#### 	LO 10.2 - [Understand the concept of functional enrichment analysis, and the statistics involved](material/L10.md#LO10.2)
		What is enrichment analysis and how is it performed
		How to define sample and population sets
		Why do we need multiple test corrections

#### 	LO 10.3 - [Interpret the results of functional enrichment analysis](material/L10.md#LO10.3)
		What can we get from enrichment analysis results
		Using functional enrichment analysis with real lists of genes


## Detailed Program

### Monday, October 8th

+ 09:30 - 10:00 Introduction to the course and self presentation of the participants
+ 10:00 - 11:00 [The High Throughput Sequencing Workflow. Designing your experiment for Differential Expression using RNA-Seq](material/L01_L02.md#LO1). [Steps in the analysis of RNA-Seq differential expression experiments](material/L01_L02.md#LO2).
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 [Interpret what are fastq files and what is their content](material/L03.md#LO3.1). [Use software like FastQC to process fastq files and produce quality reports (QC)](material/L03.md#LO3.2). 
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 [Remove low quality bases](material/L04.md#LO4.1), [Remove adaptors and other artefactual sequences from your reads](material/L04.md#LO4.2).
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 [What is a reference genome, versioning and where to obtain genomes](material/L05.md#LO5.1). [Alignment software: hisat2; bwa; salmon](material/L05.md#LO5.2). [Run an alignment: the SAM/BAM alignment format](material/L05.md#LO5.3).

### Tuesday, October 9th

+ 09:30 - 10:00 Morning Wrap-up (what have we done so far?)
+ 10:00 - 11:00 [What is a reference gene annotation, versioning and where to obtain](material/L06.md#LO6.1). [Visualizing alignments in IGV for single genes](material/L06.md#LO6.2).
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 [Use Qualimap to assess quality of alignments](material/L06.md#LO6.3).
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 [The process of generating gene counts from genome aligments](material/L07.md#LO7.1). [Use featurecounts to generate tables of gene counts](material/L07.md#LO7.2). [Use Salmon to generate counts using only the transcriptome](material/L07.md#LO7.3).
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 [Execute a pairwise differential expression analysis](material/L08.md#LO8.1)


### Wednesday, October 10th

+ 09:30 - 10:00 Morning Wrap-up (what have we done so far?)
+ 10:00 - 11:00 [Interpretation and visualization of results](material/L08.md#LO8.2). 
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 [Interpretation and visualization of results](material/L08.md#LO8.2).
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 [More complex settings: batch effects and paired samples](material/L08.md#LO8.3).
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 [Gain control over your analysis using R and Rstudio](material/L08.md#LO8.4).


### Thursday, October 11th

+ 09:30 - 10:00 Morning wrap-up (what have we done so far?)
+ 10:00 - 11:00 [Specificities of single-cell RNAseq, using the Chromium system as example](material/L09.md#LO9.1).
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 [Generate a count matrix for a single-cell RNAseq dataset](material/L09.md#LO9.1).
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 [Generate groups of cells by clustering gene expression](material/L09.md#LO9.2).
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 [Obtain marker genes for the different groups of cells](material/L09.md#LO9.3).

### Friday, October 11th

+ 09:30 - 10:00 Morning wrap-up (what have we done so far?)
+ 10:00 - 11:00 [How to extract meaning from a list of genes](material/L10.md#LO10.1).
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 [Understand the concept of functional enrichment analysis, and the statistics involved](material/L10.md#LO10.2).
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 [Interpret the results of functional enrichment analysis](material/L10.md#LO10.3).
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 Final wrap-up Session.

