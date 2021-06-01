---
layout: tutorial_page
permalink: /BiCG_2019_IA_Day2
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: Integrative Assignment Day 2
author: Hamza Farooq
modified: May 28th, 2019
---

# Integrated Assignment - Day 2

In this assignment, we will do the preprocessing steps for **copy number alteration (CNA)** analysis. In tomorrow's lab for [Module 6](https://bioinformaticsdotca.github.io/BiCG_2019_Module6_Lab) we will do the full analysis using genome-wide exome data. For today we are only focusing on chromosome 20 for the sake of time, so we can show how the preprocessing steps are used. We will get better results from our analysis tomorrow.

We will align to the human reference genome a subset of reads from a tumour sample (the HCC1395 breast cancer cell line) and a matched normal sample (a lymphoblastoid cell line derived from the same patient).

Task list:  
1) Build a bwa and bowtie2 index, as well as index our fasta  
2) Align our reads to the reference genome using bwa  
3) Convert the sam file to a sorted bam  
4) Generate a reference genome mappability file \*  
5) Generate a reference genome GC content file  
6) Calculate tumor and normal depth  
7) Identify heterozygous germline positions in the normal

\* Because this step takes about an hour, we're going to use pre-generated results for the sake of time.

We will wait to run Titan, the final step in our CNA analysis, until tomorrow so we can spend more time discussing how it works.  

---

First let's set up our working folder and create an environment variable to help navigate our paths. We're also going to make a folder to hold our logs and errors for our processes in a separate folder called "jobs":

```
cd ~/workspace
mkdir IA_tuesday
IA_HOME=/home/ubuntu/workspace/IA_tuesday
cd $IA_HOME
mkdir jobs
JOB_OUT=$IA_HOME/jobs
```

`$IA_HOME` is an environment variable. Think of it as a shortcut or nickname for our working directory, so we don't have to type out `/home/ubuntu/workspace/IA_tuesday` every time.

Now that we've setup our working directory, we can link our reference data and our tumour/normal data:

```
ln -s ~/CourseData/CG_data/IA_tuesday/data
ln -s ~/CourseData/CG_data/IA_tuesday/bams
mkdir ref
ln -s ~/CourseData/CG_data/IA_tuesday/refs/chr20_adj/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa ref/
```
This allows us to refer to the files in \~/CourseData/CG_data as though they were in our working directory $IA_HOME.

Let's start with building our reference indices for bwa and bowtie2 (Make sure to consult each respective tool's page in the future to view all the options available):
```
bwa index $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa > $JOB_OUT/bwa_index.log 2>$JOB_OUT/bwa_index.err
bowtie2-build -t 8 $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa > $JOB_OUT/bowtie2_index.log 2>$JOB_OUT/bowtie2_index.err
samtools faidx $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa
```

Currently we don't have write permissions in our `bams` file, only read permissions. This can be checked with `ls -l $IA_HOME/bams/*`. You should see something like this:
```
-r--r--r-- 1 ubuntu ubuntu  343049561 May 27 19:20 bams/HCC1395_norm.chr6.bam
-r--r--r-- 1 ubuntu ubuntu  552640940 May 27 19:22 bams/HCC1395_tum.chr6.bam
-r--r--r-- 1 ubuntu ubuntu  346959643 May 27 19:20 bams/norm.bam
-r--r--r-- 1 ubuntu ubuntu 1714996931 May 27 19:22 bams/norm.sam
-r--r--r-- 1 ubuntu ubuntu  346579564 May 27 19:21 bams/sorted_norm.bam
-r--r--r-- 1 ubuntu ubuntu    1882904 May 27 19:21 bams/sorted_norm.bam.bai
-r--r--r-- 1 ubuntu ubuntu  557279533 May 27 19:22 bams/sorted_tum.bam
-r--r--r-- 1 ubuntu ubuntu    1975568 May 27 19:22 bams/sorted_tum.bam.bai
-r--r--r-- 1 ubuntu ubuntu  558005885 May 27 19:20 bams/tum.bam
-r--r--r-- 1 ubuntu ubuntu 2767858102 May 27 19:21 bams/tum.sam
```

The `r` means read: we currently have permission to read the files in `bams`. We want to be able to edits these files, or write (`w`) them as well. We can do this with the following command:
```
chmod +w $IA_HOME/bams/*
```

Our directory should now look like this:
```
-rw-rw-r-- 1 ubuntu ubuntu  343049561 May 27 19:20 bams/HCC1395_norm.chr6.bam
-rw-rw-r-- 1 ubuntu ubuntu  552640940 May 27 19:22 bams/HCC1395_tum.chr6.bam
-rw-rw-r-- 1 ubuntu ubuntu  346959643 May 27 19:20 bams/norm.bam
-rw-rw-r-- 1 ubuntu ubuntu 1714996931 May 27 19:22 bams/norm.sam
-rw-rw-r-- 1 ubuntu ubuntu  346579564 May 27 19:21 bams/sorted_norm.bam
-rw-rw-r-- 1 ubuntu ubuntu    1882904 May 27 19:21 bams/sorted_norm.bam.bai
-rw-rw-r-- 1 ubuntu ubuntu  557279533 May 27 19:22 bams/sorted_tum.bam
-rw-rw-r-- 1 ubuntu ubuntu    1975568 May 27 19:22 bams/sorted_tum.bam.bai
-rw-rw-r-- 1 ubuntu ubuntu  558005885 May 27 19:20 bams/tum.bam
-rw-rw-r-- 1 ubuntu ubuntu 2767858102 May 27 19:21 bams/tum.sam
```

Once this is complete, we're going to use bwa to align our fastq files to the reference genome. This will take just under 3 minutes:
```
bwa mem -t 4 $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa $IA_HOME/data/HCC1395_norm.chr20.comb.fq | samtools view -bS | samtools sort > $IA_HOME/bams/HCC1395_norm.chr20.sorted.bam
bwa mem -t 4 $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa $IA_HOME/data/HCC1395_tum.chr20.comb.fq | samtools view -bS | samtools sort > $IA_HOME/bams/HCC1395_tum.chr20.sorted.bam
```

When our alignment is completed, we're going to index our bam files for easier parsing:
```
samtools index $IA_HOME/bams/HCC1395_norm.chr20.sorted.bam 2>$JOB_OUT/index_norm.err
samtools index $IA_HOME/bams/HCC1395_tum.chr20.sorted.bam 2>$JOB_OUT/index_tum.err
```

Now that the setup has been done, we can start with actually analyzing our samples. Let's set up some more environment variables to where a few key scripts are stored:
```
SCRIPTS_DIR=~/CourseData/CG_data/Module6/scripts
INSTALL_DIR=~/CourseData/CG_data/Module6/install
HMMCOPY_DIR=~/CourseData/CG_data/Module6/install/hmmcopy/HMMcopy
```

Mappability refers to how well reads map to a genome, and is roughly correlated with uniqueness of sequences in the genome (unique regions will have high mappability; repetitive regions will have lower mappability). Generating a reference genome mappability file typically takes about an hour for our analysis, even though it is focused on just chromosome 20, so we will use a pre-generated file available [here](https://github.com/bioinformaticsdotca/BiCG_2018/blob/master/IntegrativeAssignment2/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.map.ws_1000.wig). Download it to your AWS instance with the following command:
```
cd ref; curl -o Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.map.ws_1000.wig https://github.com/bioinformaticsdotca/BiCG_2018/blob/master/IntegrativeAssignment2/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.map.ws_1000.wig
```
If you would like to learn how to run it yourself, see [this file](https://github.com/bioinformaticsdotca/BiCG_2019/blob/master/IA_Day2/RefGenomeMappability.md).

In order to correct the binned read counts for GC bias (which affects read coverage), we need to generate files based on the GC coverage of the reference genome:
```
$HMMCOPY_DIR/bin/gcCounter \
    $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa \
    > $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.gc.wig
```

Now that our coverage and mappability files have been created, we can calculate our tumor and normal read depth for chromosome 20. Our analysis uses binned read counts to infer copy number. We used the HMMcopy C package to extract the number of reads in each bin for 1000 bp bins across the genome, from both the tumour and normal .bam files. The output is in the form of .wig files for the tumour and normal. These files specify the number of reads that fall into each genomic bin.
```
cd $IA_HOME;
mkdir -p $IA_HOME/hmmcopy
$HMMCOPY_DIR/bin/readCounter \
    -c 20 \
    $IA_HOME/bams/HCC1395_tum.chr20.sorted.bam > $IA_HOME/hmmcopy/HCC1395_exome_tumour.wig
$HMMCOPY_DIR/bin/readCounter \
    -c 20 \
    $IA_HOME/bams/HCC1395_norm.chr20.sorted.bam > $IA_HOME/hmmcopy/HCC1395_exome_normal.wig
```

Next we will identify heterozygous germline positions in the normal and extract allele counts in the tumour. We will use these positions to extract b-allele frequencies (BAF) from the tumour genome. Since this requires running mutationseq which calls python, we're going to add python to our PATH so that our system can use the program. We'll also use environment variables for ease of use:
```
PYTHON_DIR=$INSTALL_DIR/miniconda/miniconda2
MUSEQ_DIR=$INSTALL_DIR/mutationseq/mutationseq/museq
REF_DIR=$IA_HOME/ref
export PATH=$PYTHON_DIR/lib/python2.7/site-packages/:$PATH
export LD_LIBRARY_PATH=$PYTHON_DIR/lib:$LD_LIBRARY_PATH
```
```
mkdir -p $IA_HOME/mutationseq
$PYTHON_DIR/bin/python $MUSEQ_DIR/preprocess.py \
	reference:$REF_DIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa \
	model:$MUSEQ_DIR/models_anaconda/model_single_v4.0.2_anaconda_sk_0.13.1.npz \
	tumour:$IA_HOME/bams/HCC1395_tum.chr20.sorted.bam \
	normal:$IA_HOME/bams/HCC1395_norm.chr20.sorted.bam \
	--verbose --single --coverage 4 --threshold 0.85 --buffer_size 2G \
	--out $IA_HOME/mutationseq/HCC1395_mutationseq.vcf \
	--log_file $JOB_OUT/run_mutationseq.log \
	--config $MUSEQ_DIR/metadata.config \
	2> $JOB_OUT/run_mutationseq.err
```

The file HCC1395_mutationseq.vcf contains the raw output of MutationSeq as a [variant call format](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/) file.
Once mutationseq is completed, we'll use our custom script to transform the vcf file to a counts file:
```
$PYTHON_DIR/bin/python $SCRIPTS_DIR/transform_vcf_to_counts.py \
	--infile $IA_HOME/mutationseq//HCC1395_mutationseq.vcf \
	--outfile $IA_HOME/mutationseq/HCC1395_mutationseq_postprocess.txt
```

This is the type of file that is used when running Titan to perform CNA analysis, which we will do tomorrow.