---
layout: tutorial_page
permalink: /BiCG_2019_Module5_lab
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: BiCG Module 5 Lab
author: Jared Simpson
modified: March 5th, 2019
---

# Genome Assembly for Short and Long Reads

by [Jared Simpson](https://simpsonlab.github.io)

## Introduction

In this lab we will perform de novo genome assembly of a bacterial genome using three different sequencing technologies. You will be guided through the genome assembly. At the end of the lab you will know:

1. How to run a short read assembler on Illumina data
2. How to run a long read assembler on Pacific Biosciences or Oxford Nanopore data
3. How to improve the accuracy of a long read assembly using short reads
4. How to assess the quality of an assembly

## Data Sets

In this lab we will use a bacterial data set to demonstrate genome assembly. This data set consists of sequencing reads for Escherichia coli. E. coli is a great data set for getting started on genome assembly as it is a small genome (4.6 Mbp) with relatively few repeats, and has a high quality reference. We have provided Illumina, PacBio and Oxford Nanopore reads for this genome. The assemblies you will run later using `spades` and `canu` use a read set that is restricted to a one megabase region of the E. coli genome. This is to reduce the amount of compute time the assemblies take, so that they complete during this lab. Running a whole genome assembly uses the exact same set of commands and the results you will obtain are comparable to assembling the entire genome.

## Data Preparation

First, lets create and move to a directory that we'll use to work on our assemblies:

```
mkdir -p ~/workspace/Module5
cd ~/workspace/Module5
```

For convenience, we'll make symbolic links to the data sets that we'll work with. Run this command from the terminal, which will find all of the sequencing data sets we provided (using the `ls` command) and symlink those files into your current working directory.


```
ls ~/CourseData/CG_data/Module5/ecoli* | xargs -i ln -s {}
```

If you run `ls` you should now be able to see three files of sequencing data.

## E. coli Genome Assembly with Short Reads

Now we'll assemble the E. coli 50x Illumina data using the [spades](http://bioinf.spbau.ru/spades) assembler. Parameterizing a short read assembly can be tricky and tuning the parameters (for example the size of the *k*-mer used) is often quite time consuming. Thankfully, spades will automatically select values for its parameters, making it particularly easy to use. You can start spades with this command (it will take a few minutes to run):

```
spades.py -o ecoli-illumina-50-spades/ -t 4 --12 ecoli.illumina.50x.fastq
```

After the assembly completes, let's move the results to a new directory that we'll use to keep track of all of our assemblies:

```
mkdir -p assemblies
cp ecoli-illumina-50-spades/contigs.fasta assemblies/ecoli.illumina.50x.spades-contigs.fasta
```

We can now start assessing the quality of our assembly. We typically measure the quality of an assembly using three factors:

- Contiguity: Long contigs are better than short contigs as long contigs give more information about the structure of the genome (for example, the order of genes)
- Completeness: Most of the genome should be assembled into contigs with few regions missing from the assembly (remember for this exercise we are only assembling one megabase of the genome, not the entire genome)
- Accuracy: The assembly should have few large-scale *misassemblies* and *consensus errors* (mismatches or insertions/deletions)

We'll use `abyss-fac.pl` to calculate how contiguous our spades assembly is. Typically there will be a lot of short "leftover" contigs consisting of repetitive or low-complexity sequence, or reads with a very high error rate that could not be assembled. We don't want to include these in our statistics so we'll only use contigs that are at least 500bp in length (protip: piping tabular data into `column -t` will format the output so the columns nicely line up):

```
abyss-fac.pl -t 500 assemblies/ecoli.illumina.50x.spades-contigs.fasta | column -t
```

The N50 statistic is the most commonly used measure of assembly contiguity. An N50 of x means that 50% of the assembly is represented in contigs x bp or longer. What is the N50 of the spades assembly? How many contigs were produced?

## E. coli Genome Assembly with Long Reads

Now, we'll use long sequencing reads to assemble the E. coli genome. Long sequencing reads are better at resolving repeats and typically give much more contiguous assemblies. Long reads have a much higher error rate than short reads though, so we need to use a different assembly strategy. In this tutorial, we'll use [canu](https://github.com/marbl/canu) to assemble the 100X PacBio dataset. The canu assembly of the pacbio data should take about 15 minutes on your computer (maybe take a break while it is running).  Run this command to generate the assembly:

```
canu -fast -p ecoli-pacbio-canu -d ecoli-pacbio-auto genomeSize=1.0m -pacbio-raw ecoli.pacbio.100x.fastq
```

When it completes, copy the Pacbio assembly to our results directory:

```
# Copy the pacbio assembly from canu's directory
cp ecoli-pacbio-auto/ecoli-pacbio-canu.contigs.fasta assemblies/ecoli.pacbio.100x.canu-contigs.fasta
```

Our data set also includes an Oxford Nanopore data set. We can assemble the genome using canu, this time providing command line arguments specific to nanopore data.

```
canu -fast -p ecoli-nanopore-canu -d ecoli-nanopore-auto genomeSize=1.0m -nanopore-raw ecoli.nanopore.100x.fastq
```

Now let's copy the assembly:

```
cp ecoli-nanopore-auto/ecoli-nanopore-canu.contigs.fasta assemblies/ecoli.nanopore.100x.canu-contigs.fasta
```

## Assessing the Quality of your Assemblies using a Reference

The accuracy of the genome assembly is determined by how many misassemblies (large-scale rearrangements) and consensus errors (mismatches, insertions or deletions) the assembler makes. Calculating the accuracy of an assembly typically requires the use of a reference genome. We will use the [QUAST](http://quast.bioinf.spbau.ru/) software package to assess the accuracy of the assemblies.

Run QUAST on your three E. coli assemblies by running this command:

```
quast.py -R ~/CourseData/CG_data/Module5/references/ecoli_k12.fasta assemblies/*.fasta
```

Using the web browser for your instance, open the QUAST PDF report (Module5/quast_results/latest/report.pdf) and try to determine which of the assemblies was a) the most complete b) the most contiguous and c) the most accurate.

## Assembly Polishing

Both the nanopore and pacbio assemblies have errors in their consensus sequence as indicated by the "mismatches per 100kb" and "indels per 100kb" lines in the QUAST output. To help improve the accuracy of the assembly, we can use a post-assembly consensus improvement step called "polishing". There are many assembly polishing programs available for both pacbio data (racon, arrow) and nanopore data (nanopolish, racon). To demonstrate polishing we will use a program called `medaka` that is particularly fast and easy to run. While we're only polishing the nanopore assembly today as a demonstration, please note that the pacbio assembly could also be improved by polishing it with arrow.

We're now going to use `medaka` to improve our assembly. Medaka uses a neural network which is trained to calculate a better consensus sequence for nanopore assemblies.

```
medaka_consensus -i ecoli.nanopore.100x.fastq -d assemblies/ecoli.nanopore.100x.canu-contigs.fasta -o ecoli_medaka_polished -t 4 -m r941_flip235
```

Now we can copy the medaka assembly to our output directory:

```
cp ecoli_medaka_polished/consensus.fasta assemblies/ecoli.nanopore.100x.canu-contigs-polished.fasta
```

Now, re-run the QUAST step from above:

```
quast.py -R ~/CourseData/CG_data/Module5/references/ecoli_k12.fasta assemblies/*.fasta
```

The report will be updated in Module5/quast_results/latest/report.pdf (all versions will also be stored in their own time-stamped directories in Module5/quast_results). Did the quality of your nanopore assembly improve?
