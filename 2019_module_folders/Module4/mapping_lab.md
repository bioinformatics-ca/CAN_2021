---
layout: tutorial_page
permalink: /BiCG_2019_Module4_lab
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: BiCG Module 4 Lab - Mapping
author: Jared Simpson
modified: May 29th, 2019
---

This tutorial was created by Jared Simpson.

# Introduction

In this tutorial we will use [bwa](https://github.com/lh3/bwa) to map reads from a breast cancer cell line to the human reference genome.

# Mapping Tutorial


## Preparation

First, let's create a directory to hold our results:

```
cd ~/workspace
mkdir Module4-mapping
cd Module4-mapping
```

## Mapping using bwa-mem

bwa mem is the leading algorithm for mapping short reads to a reference genome. You can read about how it works [here](http://arxiv.org/abs/1303.3997). We'll now run bwa mem to map HCC1395 reads to the human reference genome. 

In the following command we provide bwa with the location of the reference genome - in this exercise we use the human reference genome prepared by the [1000 Genomes project](http://www.1000genomes.org/category/reference/) - and a FASTQ file containing the paired end reads for the tumour sample. In this case the paired end reads are in *interleaved* format where the two halves of a pair are in consecutive FASTQ records. The output file is in the SAM format. When mapping a whole-genome sequencing run it will take many hours to run - in this tutorial we only use a subset of the reads so this step doesn't take very long.

```
bwa mem -t 4 -p ~/CourseData/CG_data/Module4/human_g1k_v37.fasta ~/CourseData/CG_data/Module4/reads.tumour.fastq > tumour.sam
```

In this command:  

`mem` specifies the bwa algorithm to run.  
`-t 4` specifies the number of threads to use for the alignment; in this case we use 4 threads.  
`-p` indicates the reference genome (in fasta format) that we want to align to.  
We also need to specify the fastq file we want to align.  
`>` redirects the output from bwa mem to tumour.sam instead of printing it to the screen.  


## Exploring the alignments

SAM is a plain-text format that can be viewed from the command line. You can use the head command to look at the alignments:

```
head -100 tumour.sam
```

You will see the SAM header containing metadata followed by a few alignments. You can refer to the slides from the lecture to determine the meaning of each field.

In this SAM file, the reads are ordered by their position in the original FASTQ file. Most programs want to work with the alignments ordered by their position on the reference genome. We'll use [samtools](https://github.com/samtools/samtools) to sort the alignment file. To do this, we need to first convert the SAM (text) to the BAM (binary) format. We use the `samtools view -Sb` command to do this, and pipe the output directly into samtools sort.

```
samtools view -Sb tumour.sam | samtools sort -o tumour.sorted.bam
```

In this command, `-o` specifies the output file. The `|` sends the output from the first half of the command to the second half.  

The samtools view command can also be used to convert BAM to SAM:

```
samtools view tumour.sorted.bam | head -100
```

samtools also provides functions to view the alignments for a particular region of the reference genome. To do this we first need to build an index of the BAM file, which allows samtools to quickly extract the alignments for a region without reading the entire BAM file:

```
samtools index tumour.sorted.bam
```

Now that the BAM is indexed, we can view the alignments for any region of the genome:

```
samtools view tumour.sorted.bam 9:14,196,000-14,197,000
```

As a tab-delimited file, the SAM format is easy to manipulate with common unix tools like grep, awk, cut, sort. For example, this command uses cut to extract just the reference coordinates from the alignments:

```
samtools view tumour.sorted.bam 9:14,196,000-14,197,000 | cut -f3-4
```


The `samtools flagstat` command gives us a summary of the alignments:

```
samtools flagstat tumour.sorted.bam
```

We can use `samtools idxstats` to count the number of reads mapped to each chromosome:

```
samtools idxstats tumour.sorted.bam
```

This command will just display the number of reads mapped to chromosome 9:

```
samtools idxstats tumour.sorted.bam | awk '$1 == "9"'
```

## Examining alignments

You can view all of the aligned reads for a particular reference base using mpileup.
The following command will show the read bases and their quality scores at a heterozygous SNP.
The "." and "," symbols indicate bases that match the reference. 
There are 38 reads that show a "G" base at this position. 
The individual's genotype at this position is likely A/G.

```
samtools mpileup -f ~/CourseData/CG_data/Module4/human_g1k_v37.fasta -r 9:14,196,087-14,196,087 tumour.sorted.bam
```

In this command, `-f` indicates our reference (in fasta format), `-r` specifies the region we want to look at.  

Load the data into IGV by performing the following:

```
   Open IGV
   Choose 'Load from URL' from the file menu
   Type: http://##.oicrcbw.ca/Module4-mapping/tumour.sorted.bam where # is your student ID
   Navigate to 9:14,196,087
```

Notice that the alignments have high mapping quality.

This is the end of the mapping tutorial. In the remaining time you can also map reads for the matched-normal sample for this cell line. The reads are in `~/CourseData/CG_data/Module4/reads.normal.fastq`.
