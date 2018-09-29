# <a id="LO1">Learning Outcome 1: Plan your experiment using NGS technologies</a>

A good source of information for this part is [RNA-seqlopedia](http://rnaseq.uoregon.edu).

## The High Throughput Sequencing Workflow

[Sanger sequencing](https://en.wikipedia.org/wiki/Sanger_sequencing) brought about a technological revolution, as it allowed to directly read DNA molecules with relative ease and affordability. The [Human Genome Project](https://en.wikipedia.org/wiki/Human_Genome_Project) motivated further progress, leading to automated DNA-sequencing machines capable of sequencing up to 384 samples in a single batch using capillary electrophoresis. Further advances enabled the development of high throughput sequencing (HTS), also known as next generation sequencing (NGS) platforms.

At the moment, the high throughput sequencing technology most often used (by far) is Illumina. Similarly to the Sanger method, it is also based on the addition of nucleotides specifically modified to block DNA strand elongation, where each nucleotide is marked with a different color. Unlike the Sanger method, where a single DNA molecule is "read" at a time, modern illumina machines allow reading up to millions of DNA molecules simultaneously.  

## Commmon steps in most high throughput sequencing workflows:
<p>

  * Extraction and purification of the DNA template (even RNA must usually be converted to cDNA)

  * Fragmentation of the DNA template (into a size range that can be accommodated by the machine)

  * Attachment of sequencing tags (to enable reading by the machine)

  * Amplification of signal (usually trough PCR, often already in the machine)

  * Reading of signal and conversion into nucleotide bases

</p>
<br/>

![NGS Workflow](images/NGSworkflow.jpg)

## Paired-end sequencing

Many sequencing machines can read both ends of a fragment. This is called paired-end sequencing.

![Adaptor](images/paired-end.jpg)
<br/>
<br/>

## Common parameters to consider when sequencing

When sending your samples to a sequencing facility, these are the most frequent parameters to consider:
<p>

  * Single versus Paired-end sequencing

  * Read Length

  * Coverage (number of reads)

</p>
<br/>

The following links are a good source of information regarding illumina sequencing technology:
* [Illumina Sequencing by Synthesis](https://www.youtube.com/watch?&v=fCd6B5HRaZ8).
* [Elaine Mardis talk on NGS](https://www.youtube.com/watch?v=v1DbcJD4Ry0).

<br/>

## Considerations when obtaining your RNA.

The first step in a transcriptomic experiment is to obtain the RNA. After isolating total RNA from cells, one can directly sequence it. Nonetheless, the majority of the RNA in a cell is ribosomal RNA, which may need to be removed using specific kits. Moreover, total RNA also contains unprocessed immature transcripts and RNA targeted for degradation (at different stages of processing).

Therefore, unless one is interested in non-coding RNAs or other aspects related to transcription, it is usually better to apply protocols that extract the mature mRNAs (usually through the PolyA tails). Since most people are interested in coding-genes, it is more common to use mRNA-specific protocols.

Some protocols can also keep strand information. In this case, the reads have the same (or the reverse) strand as the transcribed RNA. This is particularly relevant when sequencing total RNA, noticeably to distinguish real transcripts from transcriptional activity resulting from stalled promoters or enhancers. It can also be useful to distinguish between overlapping genes.

Finally, we also need to consider the amount of material available. Are we dealing with samples with a lot of RNA (eg. cell cultures), or short amounts (eg. small tissue samples, single-cell) that are prone to amplification artifacts and presence of contaminant sequences?

## Designing your experiment for differential expression using RNAseq.

Longer read length, paired-end sequencing and strand-specific library preparation are particularly relevant to reveal gene structure. For example, on a non-model organism for which there is no genome sequenced, or the genes are poorly annotated. They are also relevant when alternative splicing is a factor to take into consideration. Discovering gene structure is a complex process and it would be the subject of an entire course on its own.

For this course, we will focus on the analysis of differential gene expression between conditions, on organisms for which gene annotation is available. Under these conditions, long reads, paired-end, and stranded library preparation methods are not as important. Therefore, for this type of experiments, we can safely go for the cheaper single-end sequencing and shorter read lengths (eg. 50bp or 76bp).

To infer genes differentially expressed between conditions, we need to obtain accurate measures of gene expression variance between the conditions. For this, we need replicates containing as much of the expected biological variance as possible. Chosing the number of replicates and depth of sequencing (number of reads) depends on the experiment. For highly controlled conditions (such as cell cultures), 2-3 replicates could be enough. In terms of coverage, 10-40 million reads should be enough to capture most "reasonably" expressed genes. Nonetheless, to be able to more accurately estimate how much is needed, one should always start from [small pilot datasets](http://scotty.genetics.utah.edu/scotty.php), although in practice this is rarely done.

At IGC we mainly use two library preparation methods (both unstranded): [Smart-seq](https://www.takarabio.com/products/next-generation-sequencing/single-cell-rna-seq/ultra-low-input-mrna-seq-for-illumina) and [QuantSeq](https://www.lexogen.com/quantseq-3mrna-sequencing/). QuantSeq is adequate for "normal" bulk samples (with many cells), it only sequences the ends of the transcripts, requiring less reads per sample (because only a small portion of the transcript is sequenced). Since it only sequences a small portion of the transcript, it can only be used for differential gene expression analysis. Smart-Seq, on the other hand, sequences full length cDNAs, and it can be used with bulk samples, as well as with samples with very low numbers of cells, including even single-cell. Specific analysis techniques are necessary for samples with very low cell numbers, which we will cover later in the course.


# <a id="LO2">Learning Outcome 2: List steps in the analysis of RNAseq differential expression experiments</a>

Steps in the analysis of RNA-Seq:
<p>

  * QC of Raw Data; ([Learning Outcome 3](#LO3))

  * Preprocessing of Raw Data (if needed); ([Learning Outcome 4](#LO4))

  * Alignment of "clean" reads to reference genome ([Learning Outcome 5](#LO5))

  * QC of Aligments ([Learning Outcome 6](#LO6))

  * Generate table of counts of genes/transcripts ([Learning Outcome 7](#LO7))

  * Differential Analysis tests ([Learning Outcome 8](#LO8))

  * Post-analysis: Functional Enrichment ([Learning Outcome 10](#L10))
</p>
<br/>


# <a id="LO3">Learning Outcome 3: Assess the general quality of the raw data from the sequencing facility</a>

## <a id="LO3.1">LO 3.1 - Interpret what are fastq files and what is their content</a>

Most high-throughput sequencing (HTS) machines output [fastq files](https://en.wikipedia.org/wiki/FASTQ_format), the “de facto” current standard in HTS. Fastq files are simply text files, where each block of information (a sequenced DNA fragment, or read) in this format is encoded as 4 lines:

	@read_identifier
	read_sequence
	+ separator line
	base_qualities

Each base has a quality character associated with it, representing how confidently the machine identified (called) the base. The probability of error per base is given as a [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score), calculated from an integer value (Q) derived from the quality character associated to the base. Useful reference values of Q include:
* Q=10 - 90% accuracy
* Q=20 - 99% accuracy
* Q=30 - 99.9% accuracy
* Q=40 - 99.99% accuracy

Although there's theoretically no limit, Q usually goes up to around 40 in recent illumina machines.

You can see a few fastq files in the folder fastq_examples:
* sample_quality_and_adaptors.fastq.gz
* sample_adaptors.fastq.gz
* 20150821.A-2_BGVR_P218_R1.sample.fastq.gz
* 20150821.A-2_BGVR_P218_R2.sample.fastq.gz

Since each fastq can have several million reads, they can become very big. Therefore, it is usual to keep them in a compressed format such as gzip. Most recent software dealing with NGS data can directly read compressed fastq files.

**QUESTION**: Uncompress the file sample_adaptors.fastq.gz and open it in a text editor. What's the identifier and sequence of the first read of that fastq file?
<details><summary>Click Here to see the answer</summary><p>
The identifier is HWI-M01876:76:000000000-AF16W:1:1101:10853:1000. The sequence is NTGTACTTCATCCGAAACTCGTGCTCATCTCTGCTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTGAT.
</p></details>
<br/>

**QUESTION**: Uncompress the files 20150821.A-2_BGVR_P218_R1.sample.fastq.gz and 20150821.A-2_BGVR_P218_R2.sample.fastq.gz and open them in a text editor. What's the identifier and sequence of the first read of both fastq files?
<details><summary>Click Here to see the answer</summary><p>
The read identifier is the same for both files (HWI-D00418:83:C7G9GANXX:2:2108:16128:97655). This is because they are readings of the same fragment, one (_R1) in the forward and another (_R2) in the reverse direction. The sequences are different, as they are readings of different parts of the fragment, and of opposing strands.
</p></details>
<br/>
<br/>



## <a id="LO3.2">LO 3.2 - Use software like FastQC to process fastq files and produce quality reports (QC)</a>

High Throughput Sequencing machines read thousands, or even millions of sequences in parallel. As you can imagine, this usually generates large fastq files, with millions of lines. Manually inspecting quality of each read is out of the question. Thus, specialized software has been developed to provide quality measures for fastq files generated by HTS machines. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular program to generate quality reports on fastq data. Running FastQC on your raw data is usually the first thing you should do once you receive a new dataset. FastQC is free and works on Windows, Mac and Linux.

**TASK**: Open a terminal. type 'fastqc' and press enter. The graphical interface of FastQC should appear. Using the FastQC interface, open the file sample_adaptors.fastq.gz and the file sample_quality_and_adaptors.fastq.gz.

**QUESTION**: What information is in a FastQC report?
<details><summary>Click Here to see the answer</summary>

A FastQC report includes, among other things:  

  * Basic statistics of the fastq file, including number of reads and sequence length

  * Per base sequence quality, displaying the boxplot distribution of the Phred Quality (Q) per base for all reads.

  * Per sequence quality scores displaying the histogram of the mean quality (Q value) of the bases of each read, for all reads

  * Per base sequence content, displaying the frequency of each nucleotide at each position of the read

  * Per sequence GC content displaying the histogram of the GC frequency of each read, for all reads

  * Sequence length distribution displaying the histogram of read lengths

  * Sequence duplication levels displaying the histograms of the number of times reads appear with exactly the same sequence

  * Overrepresented sequences (not necessarily complete reads) that appear more frequently than randomly expected

  * Adapter content indicaring the frequency of sequences of know sequencing adaptors along the length of the reads

</details>
<br/>

FastQC reports provide a series of plots that allow the user to assess the overall quality of their raw data and detect potential biases and problems. Some plots indicate distribution of base qualities along the length of reads. You can notice that, at least for illumina data, on average the quality of each base tends to decrease along the length of the read. You can also see that the reverse read (R2) is usually of worse quality than the forward read (R1). Therefore, short single-end reads usually have better average quality, and are often ready to use right out of the sequencer.

![Base Quality](images/base_quality.png) ![Tile Quality](images/tile_quality.png)

Other plots indicate biases in nucleotidic content of reads, either globally (as %GC plots), or positionally. Global bias in nucleotidic content can be useful to search for signs of contaminants. On the other hand, positional bias are useful to detect presence of artefactual sequences in your reads such as adaptors. Another insight you may obtain from this information are potential biases in the preparation of your library. For example, random hexamer priming is actually not truly random, and preferentially selects certain sequences. The currently popular transposase-based enzymatic protocol, although reasonably random, is also not completely random, and you can see this through positional bias, particularly in the beginning of reads. The presence of adaptors is a relatively common event, and therefore specific plots exist to detect the presence of the most commonly used adaptors. Finally, the presence of repetitive sequences can also suggest contaminants, PCR artifacts, or other types of bias.

![Base Bias](images/base_bias.png) ![Adaptor](images/adaptor.png)
<br/>

**QUESTION**: What are the main differences between the reports of sample_adaptors and sample_quality_and_adaptors?
<details><summary>Click Here to see the answer</summary><p>
The sample_quality_and_adaptors.fastq.gz file contains 10000 reads of 250bp, while the sample_adaptors.fastq.gz file contains 1000 reads of 76bp. The sample_quality_and_adaptors reads have a lower per base sequence quality at their end, while the reads of the sample_adaptors keep a good quality throughout. The sample_adaptors reads contain a very noticeable nucleotide positional bias particularly after position 36. sample_quality_and_adaptors also contain a bit of nucleotide positional bias, but less and only for the first 10bp. The sample_quality_and_adaptors reads display an apparently bimodal GC distribution, while the sample_adaptors reads seem closer to a single normal distribution. Finally, sample_adaptors contain a clear presence of a known Illumina adaptor after position 36 (probably the reason for the nucleotide positional bias we saw before), while sample_quality_and_adaptors contain a much smaller frequency of another Illumina adaptor (Nextera) towards the ends of the reads.
</p></details>
<br/>

**QUESTION**: What is the major difference between the two paired fastq files of the paired_example?
<details><summary>Click Here to see the answer</summary><p>
The reverse read has poorer quality bases. This is usually the case, at least for illumina. This is because the reverse reads are generated after the forward reads.
</p></details>
<br/>

**TASK**: Let's see now one simple example of the power of running tools using the command-line terminal. In a terminal window (it can be the same as you used to open fastqc, though you'll need to close the program), go to the folder fastq_examples using the ```cd``` command. Type ```fastqc *.fastq.gz``` and press enter. Inside the folder, you should now see a series of html files with FastQC reports of each of the fastq files. You can open them with the web browser by clicking on them with the mouse, or by running ```firefox *.html```.


**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Do you understand the content of a fastQ file?

  * Did you understand the difference between single-end and paired-end reads?

  * Could you run FastQC on a fastq file?

  * Can you broadly list types of information that a FastQC report contains?

  * Can you interpret information in a FastQC report to detect potential issues with data in a fastq file?

<br/>
<br/>



# <a id="LO4">Learning Outcome 4: Do simple processing operations in the raw data to improve its quality</a>

In most cases, particularly if you're sequencing short, single-end reads, the quality of your raw data is good enough to continue without any preprocessing. In fact, if you send your sequencing to an external facility, they often do these verifications and filtering for you, and you have “clean” sequences in the end. Nonetheless, it is always better to check before proceeding.

Sometimes things can go wrong, and you may need to do something about it. Some types of problems, like presence of contaminants, or some instances of positional bias will require to go back and redo the experiments. Other issues can be minimized.

## <a id="LO4.1">LO 4.1 - Remove low quality bases from your reads</a>

As you may have noticed before, reads tend to lose quality towards their end, where there is a higher probability of erroneous bases being called.

**QUESTION**: If all bases of a ficticious machine had a Q=20 (1% probability of error), what would be the probability that one 100bp read from that machine would be completely correct?
<details><summary>Click Here to see the answer</summary><p>
P(correct)=(0.99)^100 ~= 36.6%!

This serves to exemplify that many reads in current sequencing machines are likely to have at least one base incorrect.
</p></details>
<br/>
<br/>

To avoid problems in subsequent analysis, you should remove regions of poor quality in your read, usually by trimming them from the end of reads using tools such as [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Similar to FastQC, Trimmomatic is a java program, so you can use it in any operating system (such as Windows and Mac), although unlike FastQC it needs to be run only using the commandline.

Like Trimmomatic, most software for the analysis of HTS data is freely available to users. Nonetheless, they often require the use of the command line (frequently only in a Unix-like environment). User-friendly desktop software such as [CLC](https://www.qiagenbioinformatics.com/products/clc-genomics-workbench/) or [Ugene](http://ugene.net/) is available, but given the quick pace of developmpent in this area, they are constantly outdated. Moreover, even with better algorithms, HTS analysis must often be run in external servers due to the heavy computational requirements. One popular tool is [Galaxy](https://galaxyproject.org/), which allows even non-expert users to execute many different HTS analysis programs through a simple web interface. There are public instances of Galaxy where you can run your bioinformatics analysis (eg. https://usegalaxy.org, https://usegalaxy.eu). For the purpose of this course, we will run Galaxy instances locally installed in the classroom workstations. These will contain only the tools necessary to run the exercises for this course, but otherwise work very much like any other galaxy installation. After the course, you can use the public servers, or we can also set up for you a machine with Galaxy with the necessary tools.

**TASK**: Let's use Galaxy to run Trimmomatic. Open the web browser (eg. Firefox). Type [localhost:8080](localhost:8080) in the URL tab (where you put the web addresses). This means that you are accessing a galaxy instance that is running on your local machine. You should see the Galaxy interface on your web browser. The available tools are listed on the left panel, and you can search for tools by their name. Search for trimmomatic in the tool search bar. Click on the tool Trimmomatic to see the options for running the tool.

**QUESTION**: What different operations can you perform on the reads with Trimmomatic?
<details><summary>Click Here to see the answer</summary>

You can perform the following operations with Trimmomatic (either isolated, or in combination):

	* ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read

	* SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold

	* MINLEN: Drop the read if it is below a specified length

	* LEADING: Cut bases off the start of a read, if below a threshold quality

	* TRAILING: Cut bases off the end of a read, if below a threshold quality

	* CROP: Cut the read to a specified length

	* HEADCROP: Cut the specified number of bases from the start of the read

	* AVGQUAL: Drop the read if the average quality is below a specified value

	* MAXINFO: Trim reads adaptively, balancing read length and error rate to maximise the value of each read

</details>
<br/>

**TASK**:  Upload into Galaxy the fastq files from fastq_examples folder (Click on the upload icon ![upload](images/upload.jpg) on the top left of the interface). After uploading, you should now seem them on your history in the right panel. You can visualize their content by pressing the view data icon (the eye icon ![eye](images/eye.jpg)). In Galaxy, use Trimmomatic to remove low quality bases from sample_quality_and_adaptors using the default method (a 4bp window average, with a threshold of Q=20). Finally, look at the impact by running FastQC (also in Galaxy) on the trimmed reads.
<br/>

**Hint**: When uploading, Galaxy will try to guess the type of your files, but you can also explicitly specify the type of the files when uploading. For example, by default Galaxy will uncompress files when ulpoading, but if you specify that your files are fastqsanger.gz, it will keep them compressed, saving disk space.
<br/>

**QUESTION**: What was the impact of running Trimmomatic?
<details><summary>Click Here to see the answer</summary><p>
The base quality improved significantly. Nonetheless, several sequences now became smaller due to the trimming. Some became very small, and it should even be impossible to use them afterwards in the remaining of the analysis. Therefore, it is common to remove sequences that fall below a certain length (eg. 36bp).
</p></details>
<br/>

**TASK**: Although Galaxy is much simpler to use, it is hiding some details from you (and sometimes, important details are left out), so if you can, you should invest in learning the detailed command line usage of the tool, particularly if you're using it very often. In the commandline, use Trimmomatic to remove low quality bases from sample_quality_and_adaptors. Type the command ```TrimmomaticSE -phred33 sample_quality_and_adaptors.fastq.gz sample_quality_and_adaptors.trimmed.fastq.gz SLIDINGWINDOW:4:20```.  You should now have a new file with the trimmed reads. Confirm the effect of Trimmomatic by using FastQC on the file with the trimmed reads.

<br/>
<br/>


## <a id="LO4.2">LO 4.2 - Remove adaptors and other artefactual sequences from your reads</a>

Sequencing machines often require that you add specific sequences (adaptors) to your DNA so that it can be sequenced. Although sequencing facilities will generally remove these from the reads, for many different reasons, such sequences may end up in your reads, and you will need to remove them yourself. Moreover, cDNAs may contain parts of the non-genomic polyA tails that are part of mature mRNAs. Since these sequences are not part of the genome, they may prevent proper alignment and need to be removed before proceeding.

To remove these unwanted sequences, not only you have to look for the sequence in the reads, but also allow for sequencing errors, as well as the presence of incomplete sequences. Another issue of removing the adaptors is that you need to know which ones were used in your data. Since Illumina is used most of the time, these adaptors are already integrated in tools like Trimmomatic, which also take in consideration issues like reverse complement.

**TASK**: In Galaxy, Use Trimmomatic to remove adaptors from sample_adaptors.fastq.gz using Truseq3 adaptors (for this you need to select to perform an initial Illumina clip, then select the appropriate database of adaptors) and use FastQC to see the impact. Note: although Truseq3 mentions paired-end, you can also use them for single-end.

**QUESTION**: What was the impact of running Trimmomatic?
<details><summary>Click Here to see the answer</summary><p>
There was no effect in the quality of the sequences, since the original per base quality was already very good. Most reads are now only 36bp long, and the adaptor sequences are no longer present.
</p></details>
<br/>

**Hint**: In Trimmomatic you can also use your own sequences as custom adaptors. For example, in case you use uncommon adaptors, or if you want to remove polyA tails or other artefactual sequences.
<br/>

**TASK**: As you noticed, you can use Trimmommatic to do both quality and adaptor trimming. In Galaxy, use Trimmomatic to remove low quality bases from sample_quality_and_adaptors.fastq, as well as the remainings of illumina Nextera adaptors that are still left in some of the reads.

**QUESTION**: What was the impact of running Trimmomatic?
<details><summary>Click Here to see the answer</summary><p>
The base pair quality of the sequences improved and the few adaptor sequences were also removed.
</p></details>
<br/>

Paired-end data need to be handled with special care. Reads may be removed entirely if their quality is very bad (eg. if you use the MINLEN parameter in Trimmomatic). This can result in pairing information being lost, if the other member of the pair is not also removed (or placed in a special set of unpaired sequences). Software such as Trimmomatic can also take paired data as input, and handle them properly.

**TASK**: Use Trimmomatic with the 20150821.A-2_BGVR_P218 paired-end example RNA-Seq data with the default SLIDINGWINDOW parameter, as well as MINLEN of 36. Also remove adaptors (use Truseq3 paired-end adaptors). Use FastQC to evaluate the impact of the procedure. If you use trimmomatic on each individual file, you'll lose the pairing information. Therefore, you need to provide the paired data to Trimmomatic. Notice that, besides a paired fastq file, you also obtain unpaired reads that lost their pair.

**Question**: Which one has more reads - unpaired R1, or unpaired R2? Why is that?
<details><summary>Click Here to see the answer</summary><p>
The unpaired R1 has a lot more reads than unpaired R2. This is because R2 reads are usually of lower quality and are therefore more often removed.
</p></details>
<br/>
<br/>

**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Do you understand the process of quality trimming of a fastQ file?

  * Could you use Trimmomatic to improve the base quality of a fastQ file?

  * Do you understand the process of removing adaptors and other artefactual sequences from a fastQ file?

  * Could you use Trimmomatic to removing adaptors from a fastQ file?

  * Do you understand the potential issues of quality filtering and adaptor trimming in paired-end data?

<br/>
<br/>


# <a id="LO5">Learning Outcome 5: Generate alignments of processed reads against a reference genome</a>

At this stage, we've checked the quality of our raw data and performed preprocessing if necessary, so we should now be ready to use it.

## <a id="LO5.1">LO 5.1 - What is a reference genome, versioning and where to obtain genomes</a>

If the species we're using has had their genome assembled, we can now align the reads against the reference genome of our species. Genomes were (and still are) usually obtained through the efforts of large consortia, which eventually create portals making the data available for the scientific community. [ENSEMBL](http://www.ensembl.org) (in Europe) and the [UCSC genome browser](http://genome.ucsc.edu/) (in the US) emerged first as resources to display and explore the  data for the Human Genome Project, and latter agglomerated data for other model and non-model organisms, making them very convenient resources for high quality genomes.

Genome assemblies are continuously updated with new information, particularly for large eukaryotic genomes. Even the human genome, which was "completed" in 2001, is regularly being updated. More recent updates of the Human genome do not change the core sequence, but add for example alternative haplotypes for complex and highly variable regions such as the HLA. It is also very frequent to have several alternative genomes for the same species (eg. different lab strains of mice, or other model organisms).

Moreover, large genomes contain many repetitive elements, which are usually masked for secondary analysis like gene annotation. For the alignment of NGS data, it is usually recommended to use full, unmasked, sequences. It is also common to ignore alternative haplotypes, although this depends on the goals of the study.

Therefore, it is fundamental to keep track of the version of the genome used, as well as from where (and when) it was obtained. When performing analysis using resources like Galaxy, genomes are often already integrated and available for immediate use. You should always note as much information as possible about the genome you're using and, if in doubt, contact the service providers to find out missing information.



Finally, another alternative is to use transcript sequences (cDNA) directly as a reference. This is sometimes the only alternative, when full good quality genomes are not available. The presence of multiple alternative transcripts can make the alignment more difficult, but more recent approaches can actually take this information in consideration. We can also select collections of transcripts that are relevant for our analysis (eg. focusing on protein-coding mRNAs, and/or choosing a single representative transcript per gene).

**TASK**: Go to the [Ensembl](http://www.ensembl.org) website and obtain genomic fasta for Drosophila melanogaster (click on the Downloads tab on the top, and then on Download data via FTP). Also download a fasta with cDNA. See note below to understand what you need to download. Take note of the Ensembl version, as well as the version of your genome (in case later you want to integrate data that is not from Ensembl).
<br/>

**Note**: In Ensembl you will see many genomic files to download. Particularly, you'll see 3 times what appear to be the same files, but with 'rm' and 'sm'. These are masked versions of the genome (rm - repeat masked, sm - soft masked), where repetitive regions of the genome are either replaced by 'N' (in the 'rm' versions) or as lower case (in 'sm'). Repeat masking is an important process particularly for genome annotation and other applications. In RNA-Seq analysis it is usually not recommended to use repeat masked versions. You may also see the sequences of chromosomes separately, but you usually want the full genome in one single file, usually indicated as 'toplevel'. In a few highly studied species, such as Human or Mouse, there may be several alternative sequences (haplotypes) for one same region of the genome. Unless one is specifically interested in these regions, we usually do not want these alternative regions, and therefore should take only the primary assembly. For any doubt, you can read the README file in the Ensembl FTP site.
<br/>

**Note**: You need to register from where you obtained your genome, even if you know the version. For example, if you obtained your genome from UCSC, the chromosome names start with 'chr' (chr1, chr2, ...) while if you obtained from ENSEMBL, the chromosome names do not have 'chr'. Although seemingly innocent, it may make results obtained using UCSC and ENSEMBL genomes hard to integrate (even if made on the same version).  

<br/>
<br/>

## <a id="LO5.2">LO 5.2 - Alignment software: hisat; bwa; salmon</a>

To be able to align millions of short reads to a (sometimes large) reference genome, novel, more efficient, alignment methods had to be developed. The most popular are based on the [burrows-wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform), of which [bwa](http://bio-bwa.sourceforge.net/) and [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) are examples. They enable alignment of millions of reads in a few minutes, even in a common laptop.  

Methods based on the burrows-wheeler transform make assumptions to speed up the alignment process. Namely, they require the reference genome to be very similar to your sequenced DNA (less than 2-5% differences). For example, mouse data will align poorly to the human genome, although in the case of RNA-Seq this is less problematic since genes tend to be much better conserved than the rest of the genome (you would probably still bias your results to better conserved genes). Moreover, these fast alignment algorithms are not optimal, and therefore sometimes make some mistakes, although they work quite well most of the time.

![NGS Alignment Workflow](images/NGSalignmentworkflow.jpg)

Eukaryotes contain the extra complication of splicing, where your read will be spread through multiple regions of the genome (usually, different exons of the same transcript). When using small, single-end reads, this is less of a problem, since it is less likely that your reads will overlap significantly with a splice site. Nonetheless, it is a good idea to use an aligner that allows split reads. [Hisat](https://ccb.jhu.edu/software/hisat2/index.shtml) (based on bowtie) is one such splice-aware aligner (it is an update of the better known Tophat aligner). It still uses the same approach as before, but with extensions to allow splitting of reads (i.e. different parts of a read may align in different places). Recent updates of bwa (bwa mem) also allow splitting of reads, and can be used for RNA-Seq data.

![Spliced Alignment](images/SplicedAlignment.jpg)

Finally, another set of more recent approaches quickly gaining in popularity map directly against the transcriptome, without the need for a reference genome. [Salmon](https://combine-lab.github.io/salmon/) provides transcript-level estimates of gene expression without explicitly generating alignments. These methods are very fast since they rely on sequence similarity (pseudoalignments) instead of perfect alignments, and because they only search against the transcriptome. Through the use of more elaborate statistical methods they also handle the presence of different alternative splice forms that cause difficulties in the attribution of a read to a transcript. Some of these methods, such as salmon, also take explicitly in consideration bias related to differences in transcript length and nucleotide composition.

## <a id="LO5.3">LO 5.3 - Running an alignment. The SAM/BAM alignment format</a>

**TASK**: In Galaxy, run Hisat2 on the 20150821.A-2_BGVR_P218 R1 fastq file (in single-end mode) against the Drosophila genome that should be available in your Galaxy instance (you only need to set the parameter for input fastq and the genome).
<br/>

**QUESTION**:: What's the result of running Hisat2
<details><summary>Click Here to see the answer</summary><p>
The result of Hisat is an alignment file in the BAM format.
</p></details>
<br/>

The output of these aligners consist of SAM/BAM files. The [Sequence Alignment/Map (SAM) format](https://samtools.github.io/hts-specs/SAMv1.pdf) is a tabular text file format, where each line contains information for one alignment. SAM files are most often compressed as BAM (Binary SAM) files, to reduce space and allow direct access to alignments in any arbitrary region of the genome. Several tools (particularly the most recent ones) only work with BAM files. Some aligners still produce only SAM files, which may need to be converted to BAM.

![SAM Structure](images/bam_structure.png)
<br/>

**TASK**: Visualize the BAM you generated in Galaxy.
<br/>

**QUESTION**:: What is the position of the start of the first alignment in the file you generated in Galaxy (when you visualize the BAM file in Galaxy it shows you the SAM information)?
<details><summary>Click Here to see the answer</summary><p>
The first line in the SAM file (after the long header section with the lines starting with '@') corresponds to the alignment of the sequence with the identifier HWI-D00418:83:C7G9GANXX:2:1314:11646:81987. It aligns in cromosome 2L, position 69393, in the reverse orientation (that's what the 16 in the second column means), with a mapping quality of 60.
</p></details>
<br/>

The reads that do not align are also in the SAM file, but with information mentioning that they did not align the reference.

**QUESTION**:: In the description of the BAM file (result of Hisat2), you should see some text reporting a summary of the alignment. What's the overall alignment rate that is reported?
<details><summary>Click Here to see the answer</summary><p>
The overall alignment rate is 81.41%.
</p></details>
<br/>


As we mentioned before, aligners for NGS data depend on large data structures for their efficiency. These structures (like the blast databases) are built from the fasta file containing the sequence of the reference genome. This process is relatively slow and computationally intensive, although it is only necessary to do it once for each reference genome. Therefore, before aligning your reads, it is necessary to do an indexing step on the genome sequence that will be used for alignment. If using the tools on the command line, one needs to explicitly perform this step.


**TASK**: Unzip the fasta file with the Drosophila genome that you downloaded (in the terminal, move to the folder where you downloaded the file and execute ```gunzip Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz```. If you use the mouse to unzip, it may generate a strange name inside many folders, but it is the correct genome. If it is inside many folders, copy it to the folder where you downloaded the file. You may also want to rename it to the same name as the file you downloaded, but without the .gz. In the same folder where you have the unziped fasta file run the following command: ```hisat2-build Drosophila_melanogaster.BDGP6.dna.toplevel.fa Drosophila_melanogaster.BDGP6.dna.toplevel.hisat2```.

<br/>

**QUESTION**: After the command finished, can you see other files created in that same folder?
<details><summary>Click Here to see the answer</summary><p>
You can see several files, all starting with 'Drosophila_melanogaster.BDGP6.dna.toplevel.hisat2'. These contain the burrows-wheeler indexes that hisat2 will use to perform the alignment.
</p></details>
<br/>

**TASK**: After the genome indexing step, you can now perform the alignment. Make sure the file 20150821.A-2_BGVR_P218_R1.sample.fastq.gz that is inside the fastq_examples folder is in the same folder as your genome index. Now, run the following command: ```hisat2 -x Drosophila_melanogaster.BDGP6.dna.toplevel.hisat2 -U 20150821.A-2_BGVR_P218_R1.sample.fastq.gz > 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.sam```.
<br/>

This step generated a SAM file, containing the alignment information like we saw before, and which you can open using any text editor or spreadsheet.

**TASK**: We will now transform the SAM file into a BAM file. In the same terminal window where you indexed the genome and performed the alignment using Hisat, run ```samtools view -Sb 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.sam > 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.bam```.

The BAM file contains exactly the same information as the SAM file, but compressed (like a zip file). Like the SAM file, the BAM file contains information on millions of alignments. For many applications, we need to efficiently extract the alignments in certain regions of the genome (eg. the genes). For this, we can create a BAM index (NOT to be confused with the genome index), which acts like yellow pages to efficiently direct software tools to alignments in specific regions of the genome.

**TASK**: To create the BAM index, the alignments in the bam file need to be sorted by position. For this, run ```samtools sort 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.bam 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.sorted```. Finally, we can create the index ```samtools index 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.sorted.bam```. Notice now the appearance of a companion file 20150821.A-2_BGVR_P218_R1.sample.fastq.gz.sorted.bam.bai that contains the index. This file should always accompany its corresponding bam file, as programs will look for the existence of this file.
<br/>

Using services such as Galaxy, the genome indexing step is hidden from the user. When performing the alignment in Galaxy, you usually have two options: either you provide a fasta with the reference genome, or you select an available pre-built genome. When using an available genome, the indexing step was already performed, while if you provide your own fasta of the genome, an indexing step will have to be performed before the alignment step. If your genome of interest is relatively large (roughly >100Mb), it is more efficient to have it pre-built, particularly if you're reusing it often. For this, you will need to ask the persons managing the service you're using. Moreover, the transformation from SAM to BAM and the BAM indexing steps are usually also transparent in Galaxy, where the user only sees the final result.
<br/>

**TASK**: In Galaxy, run Hisat2 on the 20150821.A-2_BGVR_P218 R1 file (in single-end mode) like you did before, but using as genome the fasta for the Drosophila genome that you downloaded previously (you'll need to upload it into Galaxy). Compare the differences in the time it takes to run Hisat2.
<br/>

**QUESTION**:: Which took longer: the one with the prebuilt genome index or the one where you passed the genome fasta?
<details><summary>Click Here to see the answer</summary><p>
The alignment using the prebuilt genome index took much less time, since the genome index was already generated.
</p></details>
<br/>

Most genomes (particularly mamallian genomes) contain areas of low complexity, composed of repetitive sequences. In the case of short reads, sometimes these align to multiple regions in the genome equally well, making it impossible to know where the fragment came from. Longer reads are needed to overcome these difficulties, or in the absence of these, paired-end data can also be used. Some aligners (such as hisat or bwa) can use information from paired reads to help disambiguate some alignments. Information on paired reads is also added to the SAM/BAM file by most aligners when this data is used.
<br/>

**TASK**: In Galaxy, run Hisat2 with the 20150821.A-2_BGVR_P218 example paired-end data against the prebuilt Drosophila genome. Use the trimmed versions that you generated before with Trimmomatic.
<br/>
<br/>

In the guilgur folder, you'll have data extracted from [Guilgur et al, 2014](https://elifesciences.org/content/3/e02181). In this Drosophila melanogaster dataset, we have two conditions (WT and mut), with two replicates for each (note that nowadays, it is more common to use 3 or more replicates). To make it quick to run, we have extracted data for a very limited set of genes. This data is already of good quality, ready to align.
<br/>

**TASK**: In Galaxy, create a new empty history to avoid having too many datasets in the same history. Then, upload all R1 files from the guilgur folder (NOTE: when uploading to Galaxy, explicitly set their type to fastqsanger.gz). Next run Hisat2 on them using the Drosophila genome that is already prebuilt. After finishing, download all BAM files **and** their companion indexes (when clicking on the disk icon, you'll have two files to download, download them both).
<br/>

**Hint**: You can rename an item in you history by pressing the Edit Attributes button ![edit](images/edit.jpg). Renaming files may come in handy later.
<br/>
<br/>



**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Do you understand the concept of genome version, and where to obtain reference genomes?

  * Do you understand the need of a specialized alignment software to align millions of short RNA-Seq reads to a reference genome? Can you name at least one software that is often used for this task?

  * Could you use Hisat2 to align a RNA-Seq fastq file to a reference genome? Do you understand the basic steps involved in the use of Hisat2?

  * Do you understand what is a SAM and BAM file, the result of the alignment process?

<br/>
<br/>


# <a id="LO6">Learning Outcome 6: Assess the general quality of the alignments and detect possible problems</a>

## <a id="LO6.1">LO 6.1 - What is a reference gene annotation, versioning and where to obtain it</a>

To estimate gene expression, we need to define the genes by identifying their position in the genome. This information is stored in a hierarchical fashion (the genes, their transcripts, each transcript's exons, and so on...) in formats such as the [Generic Feature Format (GFF) files](http://gmod.org/wiki/GFF3). These consist basically of tabular text files with positions of genes (and their components) in the genome (for a specific genome version), as well as other information about the gene such as its name. Another common format used for annotations is the [BED format](http://genome.ucsc.edu/FAQ/FAQformat.html#format1).

Each gene annotation is deeply associated to one given version of the genome (because it contains positions in that genome), but the same genome version may (and usually has) several gene annotation versions. The same way one should keep in detail the version of the genome, we should also take note of the version of the gene annotation being used, and from where and when we obtained it. We also need to make sure we're using a gene version that is compatible with the genome version we used.

Gene annotations are usually complex to create, particularly for large mammalian genomes, and are a permanent work in progress (even more than the genome). Annotation of the Human genes is the work of several large groups of researchers. Other model organisms (such as the mouse) also have dedicated teams to curate their genes. Non-model organisms that are less intensively studied may suffer from having less well characterized annotations, frequently derived from other better studied organisms.

The same way Ensembl is a good source for the genome sequence, it is also a good source to obtain gene annotations. Ensembl even defined a specific derivative of the GFF format (the [GTF](http://www.ensembl.org/info/website/upload/gff.html) format) which is commonly accepted by most applications. Ensembl is updated roughly every three months. If you download the genome and gene annotation from the same ENSEMBL version, then you're sure to have compatible files, but be aware that in a few months the annotation (or even worse, the genome assembly) may change.

**QUESTION**:: Upload into Galaxy the sample GTF in the guilgur folder and inspect it (using the eye icon). What's the first gene described in the GTF?
<details><summary>Click Here to see the answer</summary><p>
It's the gene Rpn12R (common name), also known as FBgn0036465 (Flybase identifier), located in chromosome 3L, from positions 15041631 to 15042570, in the forward strand. You may notice a hierarchical structure, with genes having several transcripts, and each transcript having several exons.
</p></details>
<br/>

**TASK**: Obtain the latest Drosophila melanogaster GTF from [Ensembl](http://www.ensembl.org), similarly as how you obtained the genome fasta file.
<br/>
<br/>

## <a id="LO6.2">LO 6.2 - Visualizing alignments in IGV for single genes</a>

The SAM/BAM file contains all the information on the alignments, but it is not a practical way to visualize them for humans. To visualize the alignments along the reference genome one can use software such as [IGV](http://software.broadinstitute.org/software/igv/) or [Tablet](https://ics.hutton.ac.uk/tablet/), which work with the most common operating systems. To avoid loading all alignments simultaneously in memory, and to be able to quickly search for region-specific alignments, this software uses the BAM format (with indexes). With IGV we can have an overview of the alignments, which is particularly useful if we have specific genes of interest apriori. Usually, though, we use IGV after obtaining the final results, to confirm that certain genes (particularly important for our experiment) are really differentially expressed.

![IGV Visualization](images/IGV_visualization.jpg)
<br/>

**TASK**: Download the BAM files you generated for the guilgur dataset in Galaxy ![download](images/download.jpg). You also need to download the companion bai index files that accompany each of the BAM files. These bai index files are essential to quickly access alignments inside the BAM file. Run IGV and load the Drosophila genome fasta file as reference (in Genomes, load genome as file). Next load (in File, load file) the provided annotation file (Drosophila_melanogaster.BDGP6.85.sample.gtf - **not** the one you downloaded from ENSEMBL) inside the guilgur folder which contains information just for a subset of genes. Finally, load the BAM alignment files.
<br/>

**Hint**: To run IGV, you can open a terminal, and type 'igv'. It will then take a little bit to start.
<br/>

**QUESTION**: In IGV, look at position: 3L:15041314-15044195 (alternatively, look for the gene Rpn12R, or Fbgn0036465). What can you see?
<details><summary>Click Here to see the answer</summary><p>
You can see that the gene is expressed (there are reads aligning) in the Mutant, but not in Wild Type control samples.

![Rpnr12r](images/rpnr12r.jpg)

</p></details>
<br/>

**QUESTION**: In IGV, look at position: X:20689286-20698941 (alternatively, look for the gene run, or FBgn0003300). What can you see?
<details><summary>Click Here to see the answer</summary><p>
You can see that the gene seems slightly more expressed in the Mutant, but more noticeable you see a splicing defect in the mutant that you don't see in Wild Type control samples.

![Run](images/run.jpg)

</p></details>
<br/>

**QUESTION**: In IGV, look at position: X:5898729-5908384 (alternatively, look for the gene Act5c, or FBgn0000042). What can you see?
<details><summary>Click Here to see the answer</summary><p>
That gene is highly expressed in all samples (much greater coverage than the other genes). Looking at the coverage along the gene, it seems that samples seem to cluster by replicates instead of their genotype, since in one group of replicates the reads seem to accumulate next to the polyA.

![Act5c](images/act5c.jpg)

</p></details>
<br/>

**QUESTION**: For the gene Act5c, expression profile seems to suggests samples are clustered by replicate insted of genotype. What may this be suggesting?
<details><summary>Click Here to see the answer</summary><p>
The accumulation of reads at the ends for one group of replicates suggests a batch effect caused by RNA degradation affecting more one group of replicates than the other, irrespective of the genotype associated to the samples.
</p></details>
<br/>

**QUESTION**:  Would you be able to detect all of what you saw here using microarrays? If not, what and why?
<details><summary>Click Here to see the answer</summary><p>
You wouldn't be able to see unforeseen splicing defects or alternative transcript usage. You would also not be able to see mutations, which you can with sequencing data. It would also probably be hard to detect positional issues like we observed that suggest RNA degradation.
</p></details>
<br/>


## <a id="LO6.3">LO 6.3 - Assessing the quality of alignments</a>

After generating alignments and obtaining a SAM/BAM file, how do I know if this step went well? In fact, there are potential issues that we can only detect after we try to align against the reference genome.

One important general measure is how many (out of all reads) were properly aligned against the reference genome. In the case of bacterial sequencing one would expect >95% successful alignment, but when sequencing a mamallian genome (with many repetitive areas) it may be normal to have as low as 70-80% alignment success. RNA-Seq sequenced regions that are usually well preserved, and thus alignment rates should usually be high (more than 80-90%).

**QUESTION**:  Some aligners, such as Hisat, report the alignment rate. What were the alignment rates for the fastq files from the guilgur dataset?
<details><summary>Click Here to see the answer</summary><p>
The alignment rates were all very high (>=90%). This is because this is a selected set of reads. In the case of the example files 20150821.A-2_BGVR_P218_R1/R2.sample.fastq.gz (where reads were chosen randomly), alignment rates are around 80%, slightly higher after trimming and using paired information (which improves alignment).
</p></details>
<br/>

A low alignment rate may be caused by several reasons: the reads may not have been properly quality filtered or may contain artefactual sequence (such as adaptors and polyA tails); there may be contaminants; or an inappropriate reference genome may have been used for alignment. It can be hard to find out if there were contaminations, unless we have an idea of the possible contaminants. Finally, if we didn't use the proper genome, but there is no closer genome available, then there is not much that can be done, except perhaps trying to change parameters in the alignment software to allow for more mismatches (although this may cause biases and an increase in wrong alignments).

Another measure that can be used is the percentage of reads with duplicates (aligning exactly to the same place in the genome). Usually, duplication levels higher than 20% are not a good sign (they're a sign of low amount of sample and PCR artifacts) but again, depends on what you are sequencing and how much. In RNA-Seq it is common to have a small set of genes highly expressed, leading to the presence of duplicates. The histogram of number of duplicates per read will often look bimodal, with most reads being unique and a small subset highly duplicate (mostly from highly expressed genes). Unfortunately it is hard to distinguish PCR artifacts from highly expressed genes. When looking in IGV, PCR artifacts can be easily detected by an uneven coverage of the gene, and the presence of large stacks of duplicated reads. To be safer, one can remove duplicates, but this is not usually done in RNA-Seq, since a lot of valid information may be lost.

There are reports specific for RNA-Seq which depend on gene annotation. One report indicates how well the genes are covered by sequence, which provides a good indication of RNA integrity. One can also check how well the alignments match the known annotation. The presence of a lot of alignments outside annotated genes can mean several things: annotation is not correct (eg. if you're working with a non-model organism); there can be DNA contamination (homogeneous presence of reads in both introns and intergenic); or presence of immature RNA (more intronic, less intergenic).
<br/>
![Positional Bias](images/positional_bias.jpg)
<br/>
![Gene Coverage](images/gene_coverage.jpg)
<br/>
The same way FastQC generates reports of fastq files to assess quality of raw data, there are programs that generate global reports on the quality of BAM alignments. One popular tool for this is [qualimap](http://qualimap.bioinfo.cipf.es/).

**TASK**: In the command line, type 'qualimap'. There may be some warnings about missing packages, don't worry about those. From the Qualimap GUI, produce Rna-seq reports (File>New Analysis>RNA Seq QC) using the BAM files of the guilgur dataset. You'll need to use the example gtf file as an annotation.
<br/>

**QUESTION**: What information is in a Qualimap RNA-Seq report?
<details><summary>Click Here to see the answer</summary>

A Qualimap report includes, among other things:  

  * Number of aligned/mapped reads and other global statistics

  * Genomic origin of the reads (whether they map to exons, or introns, or intergenic regions)

  * Coverage along the genes (to identify possible positional bias)

</details>
<br/>


**QUESTION**: Can you see differences in the Qualimap reports of the different samples?
<details><summary>Click Here to see the answer</summary><p>
The Genomic origin indicates reads as being mostly exonic, although in the case of the mutant samples there is a slightly larger proportion of intronic reads. Although we know that mutant samples have splicing issues, this is not very noticeable here because the plot is dominated by highly expressed ubiquitous genes that are also of maternal origin (and thus do not have splicing issues). The greatest difference is in the plot for the coverage along the gene, but surprisingly (or not, given what we saw in IGV), the difference is between the replicates and not between genotypes. As we saw before, the libraries for the second replicate seem to have greater degradation, and thus reads get concentrated towards the polyA tail ("end" of the gene). Also note that this plot is highly dominated by the highly expressed gene Act5C.
</p></details>
<br/>


**TASK**: Run qualimap using the command line, using the command: ```qualimap rnaseq -bam mut_lib1_R1.bam -gtf Drosophila_melanogaster.BGP6.85.sample.gtf``` (you may need to adjust the file names and locations). You should see a new folder created, containing the qualimap report in html format.
<br/>

**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Do you understand the concept of genome annotation, that they are versioned, and where to obtain them?

  * Can you mention measures to evaluate the quality of your RNA-Seq alignments?

  * Could you run Qualimap to assess the quality of your RNA-Seq alignments?

  * Could you understand the output of Qualimap? Could you use it to identify potential issues in your experiment?

<br/>
<br/>


# <a id="LO7">Learning Outcome 7: Generate tables of counts using the alignment and a reference gene annotation</a>

## <a id="LO7.1">LO 7.1 - The process of generating gene counts from genome aligments</a>

To perform differential expression analysis we need to count, for each sample, how many times a different transcript/gene is read. If we align directly against the transcriptome, we just need to count the number of alignments per gene/transcript. However, if there are many alternative transcripts, aligning will become difficult. One solution may be to use just one representative transcript, or the union of all transcripts to represent the gene, although this also has issues.

What is most often done is to align against the genome, and compare the alignments (SAM/BAM) against the gene annotation (as GTF or BED). We could consider that a read counts to a gene if it overlaps with any part of the gene, but in large mammalian genomes, genes can have large introns, and it is not rare that genes overlap with each other. Moreover, the presence of DNA contamination and immature RNAs may also influence the counts. It is usually preferable that a read will count for a gene only if it overlaps to at least some part corresponding to a valid mRNA transcribed from that gene. Then, if we have strand information, we should use it to resolve other possible ambiguities.

There are stil other factors to take in consideration. What to do if a read maps equally well to multiple genome regions? This will now depends a bit on the behavior on the alignment software. Usually, these cases are marked as having a low mapping quality, so we can simply ignore them by excluding alignments with a low mapping quality. But by ignoring these cases we're losing information, and in the case of large genomes with a lot of large duplicated regions, this can be problematic. Again, if we want to use this information, we need to take into consideration what the aligner software will do. For example, bwa randomly attributes a read to one of the sites, while hisat outputs all alignmens (up to a given limit of k equally good ones).

Some counting tools will actually use the information that a read aligns to different places to estimate the likelihood that a read belongs to one or the other, depending on the local (unique) coverage. This is in fact the type of approach Salmon uses to attribute reads to transcripts. Salmon does not output an exact number of reads per transcript, but the sum of the likelihoods of reads belonging to it (eg. a read may have 60% likelihood of belonging to a transcript, and thus will count not as 1, but as 0.6).

Finally, how to avoid PCR artifacts? To be as safe as possible, we would remove duplicates to avoid PCR artifacts, and this frequently needs to be done before the counting process. Nonetheless, given that duplicates can be frequent in RNA-Seq, usually we do not remove them. Assuming that PCR artifacts occur randomly, then we should not have the same artifact in different biological replicates. In any case, for genes that are very important to us, we should always also visually check the alignments using software such as IGV.


## <a id="LO7.2">LO 7.2 - Use tools such as featureCounts to generate tables of gene counts</a>

[Featurecounts](http://bioinf.wehi.edu.au/featureCounts/) is a tool to generate gene counts from SAM/BAM alignments and GFF/GTF gene annotations. Its default behavior is to generate counts at the gene level. It assigns a read to a gene if it unambiguously overlaps at least one part of a cDNA produced by the gene (namely, exons). By default, it ignores reads mapping equally well to multiple positions.  

**TASK**: In galaxy, use featureCounts with all samples of the guilgur dataset (which are unstranded) and the sample gtf file (not the complete one you downloaded) as the annotation (leave all the other parameters as default, although take note of the options you have available).

**QUESTION**:  What are the read counts for gene Rpn12R (Fbgn0036465) in all the guilgur samples?
<details><summary>Click Here to see the answer</summary><p>

  * WT Lib1: 0

  * WT Lib2: 3

  * mut Lib1: 672

  * mut Lib2: 734

</p></details>
<br/>


**TASK**: In the commandline, run the command: ```featureCounts -a Drosophila_melanogaster.BGP6.85.sample.gtf -o mut_lib1_R1.feature.counts mut_lib1_R1.bam```. Run ```featurecounts -h``` to see all the options.
<br/>

**Hint**: featureCounts can take multiple BAM files simultaneously, generating a single table with the counts for all samples.
<br/>


## <a id="LO7.3">LO 7.3 - Use Salmon to generate counts with only the transcriptome</a>

As mentioned previously, Salmon directly matches the raw reads against a fasta with the known transcriptome, directly generating a table of "counts". Since it assigns reads to transcripts probabilistically, the result is sometimes not an integer, but a fractional number.
<br/>

**TASK**: In Galaxy, run Salmon with the guilgur data against the sample transcriptome (Drosophila_melanogaster.BGP6.88.sample.cdna.fa). To obtain gene level counts, you also need to use a table converting transcripts to genes (Drosophila_melanogaster.BDGP6.88.sample.cdna.tr_to_gene.tab). Salmon will then obtain gene counts by merging transcript counts. Notice that no SAM/BAM is generated.
<br/>

**QUESTION**:  What are the salmon read counts for gene Rpn12R (Fbgn0036465) in all the guilgur samples?
<details><summary>Click Here to see the answer</summary><p>

  * WT Lib1: 0

  * WT Lib2: 2

  * mut Lib1: 687

  * mut Lib2: 776

</p></details>
<br/>

**TASK**: Like for the aligners we used, Salmon also needs to create an index to speed up the matching process. In the terminal, run the command ```salmon index --transcripts Drosophila_melanogaster.BGP6.88.sample.cdna.fa --index  Drosophila_melanogaster.BGP6.88.sample.cdna.salmon```. Next, run the alignment using the command ```salmon quant --index  Drosophila_melanogaster.BGP6.88.sample.cdna.salmon -l A -r mut_lib1_R1.fq.gz -o mut_lib1_R1.salmon.counts```.
<br/>

**QUESTION**:  What is the result you obtain? (there should be a folder called mut_lib1_R1.salmon.counts. Inside the folder, there should be a file "quant.sf". Open that file with a text editor or spreadsheet)
<details><summary>Click Here to see the answer</summary><p>
  You obtain a table of counts, but for each transcript. The counts are fractional numbers. You also have normalized counts (per million reads), and information on the "real" transcript length and an "effective" length that can be used for normalization, which takes into account several biases.
</p></details>
<br/>


**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Do you understand the concept of obtaining a table of read counts per gene

  * Do you understand the different choices you have when counting?

  * Could you use featureCounts to generate a table of counts?

  * Could you use Salmon to produce a table of counts by using the fastq files against the transcripts?

<br/>
<br/>

# <a id="LO8">Learning Outcome 8: Generate lists of differentially expressed genes, at least for a simple pairwise comparison</a>

## <a id="LO8.1">LO 8.1 - Execute a pairwise differential expression analysis</a>

We now have a table of "raw" reads counts per gene for each sample. Since there will be different number of reads for each sample, we first need to normalize the counts. For this, we can just divide each gene count by the total number of reads for that sample (transforming it into a relative number). Since numbers are usually on the order of millions, this would make for very small numbers, so the most common measure is to divide by the total number of millions of reads (eg. divide by 5, instead of 5 million). This way, we transform "raw" counts into a relative number of "counts per million reads", or CPM.

**QUESTION**:  What are the CPM for gene Rpn12R (Fbgn0036465) in all the guilgur samples starting from the featureCounts table?
<details><summary>Click Here to see the answer</summary><p>

  * WT Lib1: 0

  * WT Lib2: 8.90

  * mut Lib1: 1389.29

  * mut Lib2: 1422.17

</p></details>
<br/>

It is also common to normalize counts according to gene length. One measure that is often used is RPKM (Reads per Kb per Million reads) in case of single-end data, or FPKM (Fragments per Kb per Million Reads) in case of paired-end data. This is simply obtained by dividing the CPM by the number of Kbs of a gene. For example, gene Rpn12R has 0.94kb (it is a single exon gene), so the RPKM in mut Lib 1 is 1185.29 (1114.17/0.94). Having a normalized RPKM value may lead us to think that we can now safely compare gene expression not only between samples, but also between genes of the same sample. Nonetheless, we should always avoid comparing expression values of different genes, because read counts depend on gene length non-linearly, and also on factors such as GC content of the gene.

But even to compare the same gene between samples, this first normalization is usually not enough. You need to be constantly aware that sequencing is a sampling experiment. If in one sample there is a gene or set of genes very highly expressed, then the amount of reads available for the rest of the genes is lower than for other samples where the highly expressed gene is not so highly expressed. So, many genes that may not be differentially expressed will appear so just because of this sampling artifact.

![Competition](images/competition.jpg)

If we assume that most genes are not differentially expressed, then the ratio of counts between samples should be the same for most genes. This assumption is the basis for many types of normalization steps. One simple process consists of obtaining a reference sample where each gene counts is the mean of that gene counts in all samples. Then, for each sample, calculate the ratio of each gene count against the gene count in the reference. For most genes, this ratio should be 1 (under our initial assumption). Therefore, the median (to avoid the issue of outliers) can be used as a sample specific normalization factor.

After normalization, the next step is to apply a statistical test for differential expression. For this we need to make a model, usually based on a statistical distribution. Given that sequencing data is based on discrete counts, most of these popular methods are based on derivations of the binomial distribution. For most experiments, the distribution that seems to fit best distribution of normalized gene expression is the negative binomial.

After normalization, we need to estimate the parameters of the distribution, as well as variance in those parameters, namely how much the (normalized) counts vary between the different samples. This variation is clearly gene dependent, since highly expressed genes vary more in terms of absolute value, and low expressed genes vary more in terms of % of gene expression (fold change). If one only looks at fold change without taking variation into account, we’re more likely to have low expressed genes as differentially expressed. Therefore, we need to accurately estimate variation per gene. But usually we have a limited number of replicates, insufficient for the accurate estimation of parameters. We thus need to find a way to produce such estimates, and we need specialized tools to do it.

There are many freely available tools to perform normalization, estimate variance and perform statistical tests for differential expression. The most commonly used, and that have been demonstrated to work well under most conditions, are the R packages [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) and [edger](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

DESeq applies the normalization with the median of the ratio of gene expression like we saw before. EdgeR applies a similar, although more sophisticated, approach (trimmed mean of M-values, or TMM in short). TMM assumes that the log fold change between any sample and a reference sample (M-value) is roughly 0 for most genes. The genes with extreme M values and extreme absolute expressions values (A) are removed (trimmed) from the calculation of the normalization factor, and a bigger weight is given to genes with less variance.

To estimate variances, a strategy that is used by edgeR and DESeq2 is to bin genes with similar expression and fit a curve, assuming that genes with similar expression have a similar variance. The parameter used to define the curve is then used as a baseline to reestimate the "true" difference of gene expression between groups of samples. Finally, the assumption that most genes are not differentially expressed is also used to rescale the variance, namely by making it constant at all bin levels.

![DESeq2 Dispersion](images/deseq2_dispersion.jpg)

![LFC reestimation](images/reestimated_lfc.jpg)

We then test each gene for differential expression, and we obtain a probability for the test. Since we test thousands of genes, some genes may get good p-values just by chance. Therefore, we need to do multiple test correction. One way is to multiply each p-value by the number of tests (a method called Bonferroni correction). This is nonetheless too strict and we usually end up not having anything differentially expressed. Other methods are less stringent while maintaining a rigorous control on the false positives. The most commonly used method is the Benjamini-Hochberg, which applies a correction proportional to the ranking of the pvalue (multiplies by number of tests and divide by rank), guaranteeing a user-specified error rate. Another way of minimizing the number of tests is to filter out the genes that have very low expression in all samples.  


## <a id="LO8.2">LO 8.2 - Visualization and interpretation of results</a>

**TASK**: In Galaxy, use DESeq2 with the featureCount tables you obtained previously for the guilgur data. Name the factor as "Genotype", with the fist factor value the "WT" and the second "Mutant". For each factor value select the appropriate results of the two replicates. Select the option to output normalized values, and leave the other parameters at default and run.

**QUESTION:** What information is in the DESeq2 result file?
<details><summary>Click Here to see the answer</summary>

	* GeneID (The identifier of the gene - in this case the Flybase identifier)

	* Base Mean (The mean normalized counts of all the samples - a measure of how much is a gene expressed)

	* log2(FC) - log2 of the Fold Change (when positive, more expressed in one group than the other, and the reverse when negative)

	* StdErr - a measure of the confidence in the true value of the estimated log2(FC)

	* Wald-Stats - A value measuring how far is the observed log2(FC) from the 0 taking the StdErr into account.

	* P-value - A value measure how likely it is to obtain the observed log2(FC) by chance.

	* P-adj - The P-value corrected for multiple testing (the value that should be used in the end)

</details>
<br/>

**QUESTION:** How many genes were detected to be differentially expressed (P-adj < 0.05)?
<details><summary>Click Here to see the answer</summary>
Only the gene FBgn0036465 (Rpn12R). The gene FBgn0003300 (run), despite having a bit more expression in the mutant, is not considered to be significantly differentially expressed.
</details>
<br/>

**QUESTION:** What information is in the DESeq2 plots?
<details><summary>Click Here to see the answer</summary>

	* PcA plot: Displays the samples through a projection that most explain the variation between samples.

	* Sample-to-Sample Distances: Like the previous plot, it displays how similar are the samples between each other.

	* Dispersion estimates - A plot displaying the approximations DESeq2 does to estimate the true log2(FC)

	* Histogram of P-values - As the name implies, it depicts the distribution of P-values for all genes.

	* MA-plot - plots the M (fold change) against the A (total expression) for all genes

	In this case, since there were only 8 genes, the plots displaying gene by gene information were basically empty.

</details>
<br/>
<br/>


Even before interpreting the results of the differential expression analysis, we should have an idea of how the samples compare to each other. For this, we can look at plots such as Principal Coordinate Analysis (PCoA) or Multi-Dimensional Scaling (MDS). The way each software implements these plots vary a bit, but in short, they provide some evidence for the relatedness of the samples. In a nutshell, samples are compared against each other, generating a matrix of distances. The, the two components explaining most of the variation between samples are calculated (each component being a linear combination of samples). Ideally, we should have the expected biological variation separated along the first component. Another important aspect to consider is how much of the variance is explained by each of the components. Again, ideally, the first component should explain as much as possible the observed variation.


**QUESTION:** How is the biological variation separated in the PCoA plot you obtained for the guilgur dataset?
<details><summary>Click Here to see the answer</summary>
The Genotype is separated clearly by the first axis in the PCoA, which explains most of the variation. Another axis seems to separate samples by the replicate. Nonetheless, this axis explains much less of the variation.

![PLOT PCoA](images/guilgur_PCA.jpg)

</details>
<br/>

Another common plot shows the hierarchical clustering of samples by explicitly displaying a heatmap with the full matrix of distances between the samples, where samples are grouped according to the distance they have to each other.

**QUESTION:** How is the biological variation separated in the heatmap plot you obtained for the guilgur dataset?
<details><summary>Click Here to see the answer</summary>
The samples sharing the genotype are clearly closer together.

![PLOT HEATMAP](images/guilgur_heatmap.jpg)

</details>
<br/>


**TASK**: In Galaxy, use DESeq2 with the salmon results you obtained previously for the guilgur data. Do the same as with the htseq-count results, but now choose as "Choice of Input Data" the option TPM values. You'll need to map transcripts to genes, and for this choose the "Gene Mapping Format" Transcript-ID and Gene-ID mapping, and select the file 'Drosophila_melanogaster.BDGP6.88.sample.cdna.tr_to_gene.tab'. Compare the results with what you obtained previously.
<br/>
<br/>

To have a dataset with a more realistic size, we will also use data from Trapnell and colleagues [(Trapnell et. al, 2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/), where the authors created an artificial Drosophila melanogaster dataset with 2 conditions and 3 replicates each, where 300 randomly selected genes (from the set of genes expressed) were perturbed in-silico. The original "raw" data and processed files can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32038).
<br/>

**TASK**: Use Hisat2 on all the R1 files from the Trapnell dataset, either in Galaxy (using the Drosophila genome that is already prebuilt) or using the command line. When uploading to Galaxy, also explicitly set their type to fastqsanger.gz to minimize the disk space used.
<br/>

**Hint**: You may want to create a new empty history to avoid having too many datasets in the same history. You can rename an item in you history by pressing the Edit Attributes button ![edit](images/edit.jpg). Renaming files may come in handy later.
<br/>

**TASK**: Run a Qualimap RNA-Seq report with one of the Trapnell BAM files (use the full Drosophila annotation). In this case you should see an ideal RNA-Seq profile, with almost all reads alignining, close to 100% exonic reads, and homogeneous coverage along the gene.
<br/>

**TASK**: Use featureCounts (in Galaxy or the command line) to generate tables of counts for the Trapnell dataset (now using the full annotation).
<br/>

**TASK**: In Galaxy, use DESeq2 to perform a pairwise comparison with the count results you obtained for the Trapnell dataset. Name the factor "Condition" and two factor variables "C1" and "C2", with 3 replicates each.

**QUESTION:** How many genes were detected to be differentially expressed (P-adj < 0.05)? Hint: you can use the filter tool in Galaxy.
<details><summary>Click Here to see the answer</summary>
There are 267 genes with adjusted p-value less than 0.05.
</details>
<br/>

DESeq2 and edgeR also show the estimates of the biological coefficient of variation (BCV), which depicts the sample variation of genes according to their expression, and illustrates the variation correction the software performed, as we discussed.

**QUESTION:** What can you see in the dispersion plot for the Trapnell dataset?
<details><summary>Click Here to see the answer</summary>
You can see a curve (in red) approximating the dispersion (y axis) along expression (x axis). As you can see, dispersion is not uniform with expression, being noticeably greater with low expression (NOTE: verify this with the Trapnell dataset). Assuming dispersion should be constant for most genes, independently of expression, then you can see dispersion being re-fitted to the curve (blue dots).

![PLOT Dispersion](images/Trapnell_dispersion.jpg)

</details>
<br/>

Finally, another type of common plot is the MA plots, which displays the log fold change between groups (M) versus the average normalized expression of genes (A). The "vulcano" plot is also commonly used, depicting logFC versus adjusted p-value. On top of these graphs it is common to signal the genes that were detected as differentially expressed under certain criteria.  

**QUESTION:** What can you see in the MA plot for the Trapnell dataset?
<details><summary>Click Here to see the answer</summary>
You can see the re-estimated log fold changes (M, in the y axis) along the range of gene expression (A, in the x axis). In red are the differentially expressed genes. Since this is an artificial dataset, these can be seen quite clearly as a group. You can also see a small group of false positives, on the positive range of fold change (and probably a few in the negative range too).

![PLOT MA](images/Trapnell_MA.jpg)

</details>
<br/>


## <a id="LO8.3">LO 8.3 - More complex settings for differential expression: paired data, batch effects</a>

So far, we just considered the simple case of pairwise comparison, where all samples are independent. But we may have cases where the samples are not independent. For example, in case of cancer, it is common (and desirable) to have tumor tissue and normal tissue for the same individual. In this case, we have paired information that needs to be taken into account in the test. There can also be other variables (eg. samples were prepared in different batches) that may confound the differential expression analysis.

**TASK**: In a first example, we have a classic paired test, with tumor and normal samples for one same patient. Let's first try a simple pairwise comparison, ignoring the pairing. In Galaxy, upload the count files for 'Tuch_*.tab' (these files are in the complex folder). Define the main Factor we're interested in (Disease), dividing it in two Factor Levels (Tumor, and Normal). Then associate samples with 'N' to the Normal factor value and with 'T' with Tumour factor value.

**QUESTION:** How can you see the effect of the variation of the patients? Hint: look at the PCA plot.
<details><summary>Click Here to see the answer</summary>
The main source of variation is related to tumour versus normal (as expected). There is nonetheless a second axis of variation related to the patient.

![PLOT Tuch Variation](images/Tuch_Patient_Variation.jpg)

</details>
<br/>


**QUESTION:** How many genes were detected to be differentially expressed (P-adj < 0.05)? Hint: you can use the filter tool in Galaxy.
<details><summary>Click Here to see the answer</summary>
There are 1482 genes with adjusted p-value less than 0.05.
</details>
<br/>

**TASK**: Now let's add the pairing information. Repeat the previous job, but now add a second factor Patient (which DESeq2 will use as 'block'). For this factor, create three factor values, one for each patient (8,33,51) and run.

**QUESTION:** How many genes were detected to be differentially expressed (P-adj < 0.05)?
<details><summary>Click Here to see the answer</summary>
Now, there are 1,943 genes with adjusted p-value less than 0.05. Because it controlled for the variation related to each patient.
</details>
<br/>
<br/>

**TASK**: In a second example, we have a classic batch effect. Here we have an experiment in Arabidopsis, where plants were subjected to a treatment with a chemical (hrcc) or mock. In Galaxy, upload the count files for 'cumbie_*.tab'. Define the main Factor we're interested in (Treatment), dividing it in two Factor Levels (hrcc, and mock). Then associate samples corresponding to each case.  

**QUESTION:** Can you see the batch effect? Hint: look at the PCA plot and sample to sample distances.
<details><summary>Click Here to see the answer</summary>
Now, the main source of variation seems to be associated to replicate number 2 (hrrc2 and mock2). The other axis with the biological component is explaining less of the variation between samples (although still 34%).

![PLOT Cumbie PCA](images/Cumbie_PCA.jpg)

![PLOT Cumbie Sample Distance](images/Cumbie_sample_distances.jpg)

</details>
<br/>

What can you do in this case? You can either control for the batch effect giving the batch as an extra factor, or you can remove that batch altogether. What is the best way? Unfortunately there is no clear answer. Having more replicates is usually better, but given that the batch seems to be responsible for the larger fraction of the variation, it might be better to remove it. In any case, if you keep the second batch (without controling for it), your results are valid, although you'll be losing genes because of the variation.

**TASK**: Try both alternatives: Try the same as before, but removing the second batch. Next, use the second batch but include another factor (batch), with three values, one for each batch, and correspond the samples to each.  

**QUESTION:** How many genes were detected to be differentially expressed (P-adj < 0.05) in each case?
<details><summary>Click Here to see the answer</summary>
In the original case with all samples, there are 752 genes with adjusted p-value less than 0.05. There were only 251 when only two samples were used. When controling for the batch, we could obtain 1220 genes. So in this case, it seems to be better to control for the batch.
</details>
<br/>
<br/>

**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Do you understand the need to normalize your table of counts, different types of normalization and underlying assumptions?

  * Do you understand some of the issues underlying the estimation of the log fold change and its variance in order to estimate differentially expressed genes?

  * Could you use DESeq2 in Galaxy to produce a table of differentially expressed genes?

  * Could you interpret the results coming from DESeq2?

  * Did you understand how to identify issues like paired data and batch effects and how to handle them using DESeq2?


<br/>
<br/>



## <a id="LO8.4">LO 8.4 - Gain control over your analysis using R and Rstudio </a>

The tools available in Galaxy are limited in terms of the ability to express more complex experimental designs. Moreover, we are limited to the plots and outputs Galaxy gives us. To have full control on our analysis, we need to go to R and explore all the flexibility that it allows.

The following exercises are designed to illustrate the use of two of the most commonly used packages for differential expression analysis in R: DESeq2 and edgeR. In the first exercise we will analyse again the Trapnell datasets, as we did before, using DESeq2 (like we did before) and edgeR in its "classic" pairwise exact test mode. We will then compare the results from both approaches.

**TASK**: Open RStudio on your computer. Using the R console, perform the steps described in this document: [Exercise 1](http://htmlpreview.github.io/?https://github.com/dsobral/ADER/blob/master/material/complex/tutorial1.html)

In the second exercise we will perform a differential expression analysis of the Tuch et al. dataset. This example aims to demonstrate the steps necessary to use edgeR's GLM functions, that are applicable both to a simple pairwise comparison, as well as more complex designs.

**TASK**: Open RStudio on your computer. Using the R console, perform the steps described in this document: [Exercise 2](http://htmlpreview.github.io/?https://github.com/dsobral/ADER/blob/master/material/complex/tutorial2.html)

**TASK**: In Rstudio, open edgeR_example2_Cumbie.R and run the commands inside.

The final example we will explore contains several factors, and one of the factors have 3 different possible values. This introduces many possibilities of experimental questions to test. We just need to decide which ones are relevant biological questions. To define which questions make sense, it is usually best to start with the use of unbiased exploratory tools such as PCoA and Hierarchical clustering.

**TASK**: In Rstudio, run edgeR_example4_GSE60450.R.


**Hint**: The [edgeR manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) is a very good source of information.
<br/>
<br/>

**NOTE**: Assess how well you achieved the learning outcome. For this, see how well you responded to the different questions during the activities and also make the following questions to yourself.

  * Could you use Rstudio to run a differential expression analysis in R?

  * Could you ... ?

  * Did you ...?

<br/>
<br/>


# <a id="LO9">Learning Outcome 9 - Understand specificies of differential gene expression in the case of single-cell RNAseq</a>

## <a id="LO9.3">LO 9.1 - Generate a count matrix for a single-cell RNAseq dataset (Chromium - 10x Genomics)</a>

## <a id="LO9.3">LO 9.2 - Generate groups of cells by clustering gene expression</a>

## <a id="LO9.3">LO 9.3 - Obtain marker genes for the different groups of cells</a>


# <a id="LO10">Learning Outcome 10 - Perform simple functional enrichment analysis and understand the concepts involved</a>

## <a id="LO10.1">LO 10.1 - How to extract meaning from a list of genes</a>

A list of genes of “interest” produced by an ‘omics experiment such as RNA-seq is mostly meaningless: gene identifiers are opaque, and we’re generally interested in understanding phenomena at the cellular and/or organismal level, rather than the gene level. To do this, we must abstract from the genes to their functions, or whatever other aspect we’re interested in studying (e.g., chromosome locations, transcription regulation networks).
In order to abstract to the functional level, we need functional descriptions of the genes in a consistent manner, i.e., we need all functional aspects to be described in the same manner for all genes that have those aspects. In other words, we need genes to be <b>annotated</b> using a <b>functional classification scheme</b>. Moreover, this scheme should typically be organized hierarchically: the fine-grained gene functions generally occur in only one or a few genes, so further abstraction is required in order to enable integration of our gene set.

There are several suitable functional classification schemes available for genes, covering different aspects and/or levels of gene function. For instance, the Enzyme Commission (EC) classification covers individual enzymatic function, whereas KEGG covers metabolic pathways. The Gene Ontology (GO) is the broadest, most detailed, and most widely used classification scheme for gene function, and consequently also the most commonly used for the analysis of differential expression results.

GO is divided into three major functional aspects, called GO types: molecular function, which covers individual gene functions; biological process, which covers how gene functions integrate into cellular and/or organismal processes; and cellular component, which covers where gene functions take place. Each of these aspects is organized as a directed acyclic graph, which is essentially a relaxed hierarchy with multi-parenting. In addition to subclass (‘is a’) relation, GO includes other relations such as ‘part of’, ‘occurs in’, and ‘regulates’. While the three aspects of GO are ‘is a’ orthogonal, molecular functions can be ‘part of’ biological processes, and both can ‘occur’ in cellular components.

Genes can be (directly) annotated to multiple GO terms, even within the same aspect. Furthermore, according to the true path rule, a gene annotated to a GO term is implicitly annotated to all ancestors of that term. For instance, a gene annotated with ‘cytochrome c oxidase activity’ is inherently annotated with ‘catalytic activity’, ‘electron carrier activity’, ‘transporter activity’, and all other GO terms in the path between them.

GO annotations of genes are available on an individual basis in most genetic databases, as well as in dedicate GO browsers such as [AmiGO](http://amigo.geneontology.org) and [QuickGO](https://www.ebi.ac.uk/QuickGO). They can also be downloaded on a genome-wide scale from [GO’s annotation repository](http://www.geneontology.org/page/download-annotations) or [BioMart](http://www.ensembl.org/biomart).
Viewing the annotations of your gene set on an individual gene basis is unfeasible and insufficient: there are too many genes to analyze manually and integrate, and even if you could, this doesn’t tell you how significant the patterns you find are.

**Task**: Go to the [Gene Ontology download page](http://www.geneontology.org/page/downloads). Under GO annotations, download the GO annotation file for <i>Drosophila melanogaster</i>. Under ontology, download the Gene Ontology itself in OBO format (which is less verbose than OWL).

**Task**: Go to [BioMart](http://www.ensembl.org/biomart) and download the GO annotations for <i>Mus musculus</i>. Select "ENSEMBL Genes 92" database, then "Mouse Genes" dataset, then under Attributes, select "Gene stable ID", and from the External section, "GO term accession".


## <a id="LO10.2">LO 10.2 - Understand the concept of functional enrichment analysis, and the statistics involved</a>

Enrichment analysis is the application of statistical tests (usually the one-tailed Fisher’s exact test) to ascertain whether a sample set of entities is enriched in relation to the overall population with respect to particular features of interest. By enriched, we mean that the sample frequency of the feature is greater than would be expected by chance given the population frequency.

The one-tailed Fisher's exact test, a.k.a hypergeometric test for over-representation, is based on the hypergeometric distribution, which is used to determine probabilities in sampling events without replacement. More concretely, the hypergeometric distribution measures the probability of obtaining <i>k</i> successes in a sample of <i>n</i> elements, given a population with <i>K</i> successes in <i>N</i> total elements. For example, we could use it to determine the probability of drawing 4 Aces in a Poker hand of 5 cards (4 successes in a sample of 5, given 4 successes in a population of 52) which is a mere 0.002%. For Fisher's test, we want to measure the probability of getting at least <i>k</i> successes, meaning we have to sum hypergeometric probabilities from <i>x</i>=<i>k</i> up to <i>x</i>=min(<i>n</i>,<i>K</i>).  

Functional enrichment analysis is the application of enrichment analysis to 'omics gene lists, which can be considered samples of the genome (or the genes covered by the experiment), with the features of interest being gene functional categories (such as GO terms). It enables us to assess how meaningful statistically are the functional patterns we observe when going from the gene level to the functional level.

Some care is needed when defining the sample and the population sets of genes:
- The sample can be all differentially expressed genes, only the under-expressed, or only the over-expressed, depending on the  biological question being addressed. It may make sense to perform enrichment analysis with all three options, as each gives you a different insight into your dataset. The genes we consider differentially expressed can simply be those with significant p-values in the differential expression test, but we can also be more strict and enforce at least 2-fold expression differences.
- The population should only include genes observed in the experiment, for which presence or absence from the sample could be accurately determined. In RNA-seq experiments, it should not be the whole genome, but only the set of genes with meaningful expression observed in either of the conditions being compared. How we define "meaningful" expression is also subject to interpretation: we can simply exclude genes with no counts, or be stricter and exclude those with very small counts as well.
- <b>The sample must be a subset of the population</b>, so we should first apply whatever criteria we deem fit to select the population, then select the sample from within that selection.

Some care is also needed when counting the frequencies and sizes of the sample and population:
- The sample and population frequencies (<i>k</i> and <i>K</i>) should be the number of genes in the respective set annotated with the GO term being tested or any of its descendants (according to the true path rule).
- The sample and population sizes (<i>n</i> and <i>N</i>) should be the total number of genes in the respective set that have at least one GO annotation under the GO type of the term being tested. That is to say, we should only count genes that have determined functions, as otherwise we cannot be certain whether they should be counted as successes for the test.

Most of our statistical tests − including Fisher’s exact test − rely on controlling type I errors. When we accept an event/observation as significant because it has a p-value of 0.001, we are accepting that statistically, one time in a thousand, we’ll be wrong − the event/observation in question will be the product of chance alone. This is a problem when we perform multiple related tests, as the chance of getting a statistical “extreme” in at least one of them will be greater the more tests we perform. Because GO enrichment analysis relies on performing hundreds (or sometimes thousands) of Fisher’s tests, we must correct the statistics for multiple testing.

There are two families of multiple test corrections: the family-wise error rate (FWER) and the false discovery rate (FDR). In the former, we control the probability of making at least one false discovery, which is a conservative but “safe” approach that produces corrected p-values. In the latter, we control the ratio of false discoveries, which is a more powerful but less “safe” approach. It produces q-values, which indicate the ratio of false discoveries you are accepting if you reject the null hypothesis.

There are several GO enrichment analysis tools available, for instance:
- Webtools: [GOrilla](http://cbl-gorilla.cs.technion.ac.il/), [GO’s own tool](http://www.geneontology.org/page/go-enrichment-analysis)
- Galaxy/Command Line tools: [GOEnrichment](https://github.com/DanFaria/GOEnrichment) (IGC Galaxy); [Ontologizer](http://ontologizer.de/) (Galaxy test toolbox)
- R tools: gsea, GOstats, topGO

For some of these tools, you have to provide the version of GO and the GO annotations you want to use, in addition to the population and sample sets of genes. For other, you only need to provide the latter, and identify the organism in question. The first set of tools should be preferred, as they give you greater control over what you're testing with. In either case, the sample and population frequencies will be computed automatically by the tool. The tool should exclude from both sets the genes that don’t have GO annotations (of the GO type being tested) and perform the tests independently for each GO type, but not all tools do this.

**Task**: Picking up the differential expression results from the [Trapnell et al](https://raw.githubusercontent.com/dsobral/ADER/master/material/complex/trapnell_GSE32038_gene_exp.diff) example with 300 random differentially expressed <i>Drosophila melanogaster</i> genes, define the sample set and population set of genes for performing functional enrichment analysis. You can use your own EdgeR or DESeq results for this dataset instead of the published differential expression results. You can do the filtering and selection of genes either in a spreadsheet or in Galaxy.

**Task**: Using the sample and population files you generated in the previous task, as well as the GO file and <i>Drosophila melanogaster</i> GO annotation file you downloaded earlier, perform functional enrichment analsysis using the GOEnrichment tool in Galaxy, with "Summarize Output" set to off and otherwise default options. Are there significantly enriched terms at 0.01 significance without multiple test corrections? And with the correction?

**NOTE**: Assess how well you achieved the learning outcome by asking yourself the following questions:

  * Do you understand the concept of enrichment analysis and the underlying statistical test?

  * Could you define a population set and a sample set of genes from an RNA-seq experiment, after differential expression tests have been applied?

  * Why do we need to correct for multiple tests?

  * What is the difference between a p-value, a corrected p-value, and a q-value?


## <a id="LO10.3">LO 10.3 - Interpreting the results of functional enrichment analysis</a>

What we can get out of functional enrichment analysis results hinges heavily on what we put into them, i.e., on the biological context of our experiment, and the biological question(s) we are seeking to address. The clearer the question, the more straightforward it should be to interpret the results.
In general, functional enrichment analysis can be used for:
- Validation (e.g., of a protocol for extracting membrane proteins)
- Characterization (e.g., of the effects of a stress in an organism)
- Elucidation (e.g., of the functions impacted by the knock-out of a transcription factor)

With respect to the results, it is essential to keep in mind that **statistically significant does not mean biologically meaningful**.
On the one hand, we can have functional enrichment of functional aspects that seem too broad to derive any relevant conclusion, or that appear to be unrelated to the experiment in question. You should look at these with a critical eye − there may be some underlying meaning that is not readily apparent, even in the case of very generic terms such as “binding” or “single organism process”. On the other hand, aspects that are too specific may not be very interesting. In the extreme case of a biological process associated with a single gene in a given organism, if that gene appears in the study set, it is likely to be statistically enriched (if the study set is relatively small in comparison with the population set), but that doesn’t give us any insight into the study set as a whole. In general, we’re interested in GO terms that are sufficiently generic to integrate a significant part of our dataset, but sufficiently specific to give us some conclusive insights.
Despite multiple test corrections, sporadic outliers may occur. After all, we’re making a statistical test (of enrichment) on top of another (of differential expression) which in turn is preceded by a statistical normalization. Even though we’re comfortable with the assumptions and p-values in each individual step, the likelihood of error propagates across the steps, and even fine differences in each step can affect the final results.
You should also keep in mind that enrichment analysis is qualitative, rather than quantitative: you are treating genes as either “on” or “off” (be “on” differentially expressed, overexpressed, or underexpressed) and consequently only assessing which functional aspects are statistically affected, rather than by how much they are affected.

The p-values (or q-values) are mostly meaningless for interpreting the results—their role is just to identify enriched terms. The sample frequency and the semantics of the term (its definition, its place in the ontology's structure) are the aspects you should draw upon for biological interpretation. While the former can be checked on results tables, the latter requires graph views of the results. Because of GO’s hierarchical structure, we may get related enriched terms with different levels of specificity, and we should consider them together as a cluster when drawing conclusions. These clusters may not be readily apparent from a results table, but are easy to detect in a graph view of the results.

The size and complexity of GO often result in huge results tables and graphs (particularly for <i>biological process</i>), which be quite overwhelming to analyze.
We can reduce the number of tests perfomed by excluding unnecessary tests:
- Singletons, i.e., terms that occur in a single sample gene, may be enriched but are not very useful for understanding our dataset in an integrated manner, and thus can generally not be tested.
- Testing a superclass when its sample frequency is the same as its subclass's is redundant and unnecessary. The superclass can only be enriched if the subclass also is, and we gain neither specificity nor integration by testing the superclass.
Reducing the number of tests performed has the added benefit of reducing the multiple-test correction factor, which in the case of redundant tests would be artificially inflated (since the tests are fully dependent).

We can also reduce the number of tests performed by using GO slims, which are "trimmed" version of GO with only broad terms, instead of using the full GO. Different GO slims are available for different taxa. Using GO slims will greatly simplify the results. However, it will also lead to substantial loss in specificity, and the enriched terms may be more generic than we would like. Furthermore, using GO slims requires first converting the GO annotations from full GO to GO slim (which you can do using GOSlimmer in Galaxy).

Instead of using GO slims, we may simplify/summarize the enrichment analysis results a posteriori, using one of three strategies:
- The family-based clustering algorithm integrated into [GOEnrichment](https://github.com/DanFaria/GOEnrichment) which reduces complexity while keeping branch information, but loses some specificity
-The semantic similarity-based [REVIGO](http://revigo.irb.hr/) tool, which not only loses specificity but may merge branches
- An <a>ad hoc</a> filter, as long as it is not based on the enrichment analysis p-values, but rather on fixed features of GO (such as vertical or horizontal cuts of the ontology). Here the rationale is that your initial enrichment analysis is exploratory, aimed at discovering what features there are in your dataset, and you are refining it to focus on the feature that matter to you.

**Task**: Pick up the differential expression results file for [mouse brain vs. heart](https://raw.githubusercontent.com/dsobral/ADER/master/material/complex/mouse_brain_vs_heart.txt). Create a population file and two sample files, one for overexpressed genes and another for underexpressed ones.

**Task**: Run GOEnrichment in Galaxy using both sample files, as you did previously, using the Mouse annotation file you got from BioMart earlier. Run the program with the "Summarize Output" parameter set to off (and otherwise default parameters), then analyze the biological process results tables and graph files.

**Task**: Repeat the previous runs, but this time with "Summarize Output" set to on. Compare the results tables and graph files.

**Task**: Download the generic [GO Slim](http://www.geneontology.org/ontology/subsets/goslim_generic.obo). Use the GOSlimmer tool in Galaxy to convert the mouse annotations from GO to GO slim. Then repeat the GOEnrichment runs using GO Slim and the corresponding annotations. How do the results compare in terms of simplicity and specificity?
