---
layout: tutorial_page
permalink: /CAN_2021_Module6_Lab
title: BiCG 2019 Module 6 Lab
header1: Workshop Pages for Students
header2: Cancer Analysis 2021
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/CAN_2021
description: Cancer Analysis 2021 module 6Lab  - Gene Expression
author: Florence Cavalli
---

# Lab Module 6 - Gene Expression Profiling

## Background
The PCA3 gene plays a role in Prostate Cancer detection due to its localized expression in prostate tissues and its over-expression in tumour tissues. This gene expression profile makes it a useful marker that can complement the most frequently used biomarker for prostate cancer, PSA. There are cancer assays available that test the presence of PCA3 in urine. 


In this assignment, we will be using a subset of the [GSE22260 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22260), which consists of 30 RNA-seq tumour/normal pairs,

Experimental information and other things to keep in mind:
- The libraries are polyA selected. 
- The libraries are prepared as paired end. 
- The samples are sequenced on a Illumina Genome Analyzer II (this data is now quite old). 
- Each read is 36 bp long 
- The average insert size is 150 bp with standard deviation of 38bp. 
- We will only look at chromosome 9 in this exercise. 
- The dataset is located here: [GSE22260 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22260)
- 20 tumour and 10 normal samples are available 
- For this exercise we will pick 3 matched pairs (C02, C03, C06 for tumour and N02, N03, N06 for normal). We can do more if we have time. 

## Objectives
We overall want to assess the prostate cancer specific expression of the PCA3 gene. 
For this we will:
- Obtain data and references 
- Align the RNA-seq fastq files
- Quantify the transcripts expression
- Perform differential expression analysis



## PART 1: Obtaining Data and References
Goals:
- Obtain the files necessary for data processing 
- Familiarize yourself with reference and annotation file format 
- Familiarize yourself with sequence FASTQ format 

## Setup

First login into the server, and then enter the workspace directory:
In order to keep our files in one location, we're going to create a new directory for this module and enter it:

```
cd  ~/workspace
mkdir -p ~/workspace/Module6/
mkdir -p ~/workspace/Module6/Module6_Lab/
export RNA_LAB=~/workspace/Module6/Module6_Lab
```
## Obtain reference, annotation and data files and place them in the Module6_Lab directory
Note: when initiating an environment variable, we do not need the $; however, everytime we call the variable, it needs to be preceded by a $.

```
echo $RNA_LAB
cd $RNA_LAB
cp ~/CourseData/CAN_data/Module6/data.zip .
unzip data.zip
```

Check the files and set up the following environment variables

```
ls
ls fasta/*
ls refs/*
export RNA_DATA_DIR=$RNA_LAB/fasta
export RNA_REFS_DIR=$RNA_LAB/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/Homo_sapiens.GRCh38.dna.chromosome.9
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.86.chr9.gtf
export RNA_ALIGN_DIR=$RNA_LAB/hisat2
```

**Q1)** How many items are there under the “refs” directory (counting all files in all sub-directories)? 

**A1)** The answer is 12. Review these files so that you are familiar with them.

```
cd $RNA_LAB/refs/
tree
find *
find * | wc -l
```
What if this reference file was not provided for you? How would you obtain/create a reference genome fasta file for chromosome 9 only. How about the GTF transcripts file from Ensembl? How would you create one that contained only transcripts on chromosome 9?

**Q2)** How many exons does the gene PCA3 have?

**A2)** The answer is 4. Review the GTF file so that you are familiar with it. What downstream steps will we need this file for? What is it used for?

```
cd $RNA_LAB/refs
grep -w "PCA3" Homo_sapiens.GRCh38.86.chr9.gtf
grep -w "PCA3" Homo_sapiens.GRCh38.86.chr9.gtf | egrep "exon"
grep -w "PCA3" Homo_sapiens.GRCh38.86.chr9.gtf | egrep "exon" |wc -l
```

**Q3)** How many cancer/normal samples do you see under the data directory?

**A3)** The answer is 12. 6 normal and 6 tumor samples.
```
cd $RNA_LAB/fasta/
ls -l
ls -1 | wc -l
```

NOTE: The fasta files you have copied above contain sequences for chr9 only. We have pre-processed those fasta files to obtain chr9 and also matched read1/read2 sequences for each of the samples. You do not need to redo this.

**Q4)** What sample has the highest number of reads?

**A4)** The answer is that 'carcinoma_C06' has the most reads (288428/2 = 144214 reads).
An easy way to figure out the number of reads is to make use of the command ‘wc’. This command counts the number of lines in a file. Keep in mind that one sequence can be represented by multiple lines. Therefore, you need to first grep the read tag ">" and count those.

```
cd $RNA_LAB/fasta/
head carcinoma_C02_read1.fasta
```
```
> HWUSI-EAS230-R:2:88:885:1584#0/1
GCTCTTCGGTTCTTTCCTTCTTCAAGTGGTATGCTC  
> HWUSI-EAS230-R:2:15:691:382#0/1
GGATTTTGACAAATCCTTATCTCCGGCCACCCCATA  
> HWUSI-EAS230-R:2:52:111:1052#0/1
GCTGGAAAGCCACCAAGATGCTGACATTGAAGACTT  
> HWUSI-EAS230-R:2:28:1640:236#0/1
CTCTGAGACGTCACCAAGTGCGGCGCCGGCAGCCTG  
```

Running this command only gives you 2 x read number (*note: replace `YourFastaFile.fasta` with the actual name of a fasta file*):
```
cd $RNA_LAB/fasta/
wc -l YourFastaFile.fasta
wc -l *
```
Or  for the number of reads directly (*note: replace `YourFastaFile.fasta` with the actual name of a fasta file*):
```
egrep ">" YourFastaFile.fasta |wc -l
```

## PART 2: Data alignment
Goals:
- Familiarize yourself with HISAT2 alignment options 
- Perform alignments 
- Obtain alignment summary 

**Q5)** Create HISAT2 alignment commands for all of the six samples and run alignments
```
echo $RNA_ALIGN_DIR
mkdir -p $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR

hisat2 -p 8 --rg-id=carcinoma_C02 --rg SM:carcinoma --rg LB:carcinoma_C02 -x $RNA_REF_INDEX --dta -f -1 $RNA_DATA_DIR/carcinoma_C02_read1.fasta -2 $RNA_DATA_DIR/carcinoma_C02_read2.fasta -S ./carcinoma_C02.sam
hisat2 -p 8 --rg-id=carcinoma_C03 --rg SM:carcinoma --rg LB:carcinoma_C03 -x $RNA_REF_INDEX --dta -f -1 $RNA_DATA_DIR/carcinoma_C03_read1.fasta -2 $RNA_DATA_DIR/carcinoma_C03_read2.fasta -S ./carcinoma_C03.sam
hisat2 -p 8 --rg-id=carcinoma_C06 --rg SM:carcinoma --rg LB:carcinoma_C06 -x $RNA_REF_INDEX --dta -f -1 $RNA_DATA_DIR/carcinoma_C06_read1.fasta -2 $RNA_DATA_DIR/carcinoma_C06_read2.fasta -S ./carcinoma_C06.sam

hisat2 -p 8 --rg-id=normal_N02 --rg SM:normal --rg LB:normal_N02 -x $RNA_REF_INDEX --dta -f -1 $RNA_DATA_DIR/normal_N02_read1.fasta -2 $RNA_DATA_DIR/normal_N02_read2.fasta -S ./normal_N02.sam
hisat2 -p 8 --rg-id=normal_N03 --rg SM:normal --rg LB:normal_N03 -x $RNA_REF_INDEX --dta -f -1 $RNA_DATA_DIR/normal_N03_read1.fasta -2 $RNA_DATA_DIR/normal_N03_read2.fasta -S ./normal_N03.sam
hisat2 -p 8 --rg-id=normal_N06 --rg SM:normal --rg LB:normal_N06 -x $RNA_REF_INDEX --dta -f -1 $RNA_DATA_DIR/normal_N06_read1.fasta -2 $RNA_DATA_DIR/normal_N06_read2.fasta -S ./normal_N06.sam
```

Convert sam alignments to bam. How much space did you save by performing this conversion?

```
samtools sort -@ 8 -o carcinoma_C02.bam carcinoma_C02.sam
samtools sort -@ 8 -o carcinoma_C03.bam carcinoma_C03.sam
samtools sort -@ 8 -o carcinoma_C06.bam carcinoma_C06.sam
samtools sort -@ 8 -o normal_N02.bam normal_N02.sam
samtools sort -@ 8 -o normal_N03.bam normal_N03.sam
samtools sort -@ 8 -o normal_N06.bam normal_N06.sam
ls -lah 
```

Merge the bams for visualization purposes

```
cd $RNA_LAB/hisat2
java -Xmx2g -jar /usr/local/picard/picard.jar MergeSamFiles OUTPUT=carcinoma.bam INPUT=carcinoma_C02.bam INPUT=carcinoma_C03.bam INPUT=carcinoma_C06.bam
java -Xmx2g -jar /usr/local/picard/picard.jar MergeSamFiles OUTPUT=normal.bam INPUT=normal_N02.bam INPUT=normal_N03.bam INPUT=normal_N06.bam
```

**Q6)** How would you obtain summary statistics for each aligned file?

**A6)** There are many RNA-seq QC tools available that can provide you with detailed information about the quality of the aligned sample (e.g. FastQC and RSeQC). However, for a simple summary of aligned reads counts you can use samtools flagstat. You can also look for the logs generated by TopHat. These logs provide a summary of the aligned reads.

```
cd $RNA_LAB/hisat2
samtools flagstat carcinoma_C02.bam > carcinoma_C02.flagstat.txt
samtools flagstat carcinoma_C03.bam > carcinoma_C03.flagstat.txt
samtools flagstat carcinoma_C06.bam > carcinoma_C06.flagstat.txt

samtools flagstat normal_N02.bam > normal_N02.flagstat.txt
samtools flagstat normal_N03.bam > normal_N03.flagstat.txt
samtools flagstat normal_N06.bam > normal_N06.flagstat.txt

more carcinoma_C02.flagstat.txt
grep "mapped (" *.flagstat.txt
```

Note: ***more** is a unix command to view the content of a file on screen. Press "enter" to see more lines. Press "q" to exit*

## PART 3: Expression Estimation
Goals:
- Familiarize yourself with Stringtie options 
- Run Stringtie to obtain expression values 
- Obtain expression values for the gene PCA3 

Create an expression results directory, run Stringtie on all samples, and store the results in appropriately named subdirectories in this results dir
```
cd $RNA_LAB/
mkdir -p expression/stringtie/ref_only/
cd expression/stringtie/ref_only/

stringtie -p 8 -G $RNA_REF_GTF -e -B -o carcinoma_C02/transcripts.gtf $RNA_ALIGN_DIR/carcinoma_C02.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o carcinoma_C03/transcripts.gtf $RNA_ALIGN_DIR/carcinoma_C03.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o carcinoma_C06/transcripts.gtf $RNA_ALIGN_DIR/carcinoma_C06.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o normal_N02/transcripts.gtf $RNA_ALIGN_DIR/normal_N02.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o normal_N03/transcripts.gtf $RNA_ALIGN_DIR/normal_N03.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o normal_N06/transcripts.gtf $RNA_ALIGN_DIR/normal_N06.bam

ls
tree *
```

**Q7)** How do you get the expression of the gene PCA3 across the normal and carcinoma samples?

**A7)** To look for the expression value of a specific gene, you can use the command ‘grep’ followed by the gene name and the path to the expression file

```
cd $RNA_LAB/expression/stringtie/ref_only
```
Have a look at a transcripts.gtf file
```
more carcinoma_C02/transcripts.gtf
```
(Look for Ensembl Id (ENSGxxxxxxxxxxx)  for PCA3 in /workspace/Module6/Module6_Lab/refs/Homo_sapiens.GRCh38.86.chr9.gtf)
```
grep ENSG00000225937 ./*/transcripts.gtf | cut -f1,9 | grep FPKM
```


## PART 4: Differential Expression Analysis
Goals:
- Perform differential analysis between tumor and normal samples 
- Check if PCA3 is differentially expressed 

```
mkdir -p $RNA_LAB/de/ballgown/ref_only/
cd $RNA_LAB/de/ballgown/ref_only/
```

Perform carcinoma vs. normal comparison, using all samples, for known (reference only mode) transcripts:
First create a file that lists our 6 expression files, then view that file, then start an R session where we will examine these results:
```
printf "\"ids\",\"type\",\"path\"\n\"carcinoma_C02\",\"carcinoma\",\"$RNA_LAB/expression/stringtie/ref_only/carcinoma_C02\"\n\"carcinoma_C03\",\"carcinoma\",\"$RNA_LAB/expression/stringtie/ref_only/carcinoma_C03\"\n\"carcinoma_C06\",\"carcinoma\",\"$RNA_LAB/expression/stringtie/ref_only/carcinoma_C06\"\n\"normal_N02\",\"normal\",\"$RNA_LAB/expression/stringtie/ref_only/normal_N02\"\n\"normal_N03\",\"normal\",\"$RNA_LAB/expression/stringtie/ref_only/normal_N03\"\n\"normal_N06\",\"normal\",\"$RNA_LAB/expression/stringtie/ref_only/normal_N06\"\n" > carcinoma_vs_normal.csv

more carcinoma_vs_normal.csv 
```

```
R
```

The file [Module6_Lab_ballgown.R](https://raw.githubusercontent.com/bioinformatics-ca/CAN_2021/main/Module6/Module6_Lab_ballgown.R) will help you to run differential expression analysis with the Ballgown package. Copy and paste the commands from this file into your terminal once R has started (you should see ">" as a prompt).

Have a look at the Ballgrown package and manual document as well
https://www.bioconductor.org/packages/release/bioc/html/ballgown.html

**Q8)** Are there any significant differentially expressed genes? What about PCA3? 

**A8)** Due to the small sample size, the PCA3 signal is not significant at the adjusted p-value level. You can try re-running the above exercise on your own by using all of the samples in the original data set. Does including more samples change the results?

**Q9)** What plots can you generate to help you visualize this gene expression profile

**A9)** The CummerBund package provides a wide variety of plots that can be used to visualize a gene’s expression profile or genes that are differentially expressed. Some of these plots include heatmaps, boxplots, and volcano plots. Alternatively, you can use custom plots using ggplot2 command or base R plotting commands such as those provided in the supplementary tutorials. Start with something very simple such as a scatter plot of tumor vs. normal FPKM values.

**see [Module6_Lab_plots.R](https://raw.githubusercontent.com/bioinformatics-ca/CAN_2021/main/Module6/Module6_Lab_plots.R) for plotting options**.

When you finished running the commands from Module6_Lab_ballgown.R, you should have exited R (no more ">" prompt). Type R into the terminal again to restart R:
```
R
```

If you ever get stuck in R and want to return to the terminal, type `q()` to exit R.

Your plots will be at http://##.oicrcbw.ca/Module6/Module6_Lab/de/ballgown/ref_only/Module6_Lab_Supplementary_R_output.pdf

Remember to replace "##" with your instance number.

Compare your plots to the expected plots: [Module6_Lab_Supplementary_R_output.pdf](https://github.com/bioinformatics-ca/CAN_2021/raw/main/Module6/Module6_Lab_Supplementary_R_output.pdf)

