---
layout: tutorial_page
permalink: /BiCG_2019_Module6_Lab-DataProcessing
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
---

## Introduction

As mentioned in the lab, the data retrieval and some of the preprocessing steps have been completed for you in order to save time. This page provides instructions on how the retrieval and preprocessing was done, in case you wish to reproduce the results in your own time. 

## Environment

In this section, we will set some environment variables to help facilitate the execution of commands. These variables will store the location of some important paths we will need.

>Note: these environment variables will only persist in the current session. If you log out and log back into the server, you will have to set the variables again.

Set your analysis directory:

```
CNA_WORKSPACE=~/workspace/Module6
```

The preprocessing helper scripts can be downloaded from the wiki. Set the directory of where these scripts are stored:

```
CNA_SCRIPTS_DIR=$CNA_WORKSPACE/scripts
```

Set your install directory:

```
INSTALL_DIR=$CNA_WORKSPACE/install
```

Set the directory where HMMcopy is installed:

```
HMMCOPY_DIR=$INSTALL_DIR/hmmcopy/HMMcopy
```

## Array Data

In this section, we cover the retrieval and preprocessing of data used in the lab section 
"Analysis of CNAs using Arrays".

### Download the HCC1395 Affymetrix SNP 6.0 Microarray Data

CEL files containing the probe intensity data for cell line HCC1395 can be downloaded from the NCBI Geo website as follows:

```
mkdir -p $CNA_WORKSPACE/data/HCC1395/snp6
cd $CNA_WORKSPACE/data/HCC1395/snp6
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM888nnn/GSM888107/suppl/GSM888107.CEL.gz
gunzip GSM888107.CEL.gz
```

Details of subsequent data processing are in the lab module.

## Sequencing Data

In this section, we cover the retrieval and preprocessing of data used in the lab section 
"Analysis of CNAs using Sequencing Data"

### Download the HCC1395 Illumina Sequencing Data

Illumina sequencing data for the HCC1395 tumour-derived cell line and normal-derived lymphoblasoid cell line can be downloaded [here](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data).

Since the link contains unaligned bam files, we have downloaded and aligned the exome data for you. Link the aligned exome data to your data directory:

~~~bash
ln -s /home/ubuntu/CourseData/CG_data/data/HCC1395/exome $CNA_WORKSPACE/data/HCC1395
~~~

### Calculate Tumour and Normal Read Depth

Copy number prediction is based in large part on read depth.  The genome is sub-divided into bins (in this case equal-sized bins of 1000 bp) and the number of reads per bin is counted. 

We extract this information from the tumour and normal bam files using the HMMcopy C package, and store it in `.wig` format.  We use the `-c` option to specify the chromosomes of interest:

```
mkdir -p $CNA_WORKSPACE/analysis/exome/hmmcopy
cd $CNA_WORKSPACE/analysis/exome/hmmcopy
$HMMCOPY_DIR/bin/readCounter \
	-c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y \
	$CNA_WORKSPACE/data/HCC1395/exome/HCC1395_exome_tumour.bam > $CNA_WORKSPACE/analysis/exome/hmmcopy/HCC1395_exome_tumour.wig
$HMMCOPY_DIR/bin/readCounter \
	-c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y \
	$CNA_WORKSPACE/data/HCC1395/exome/HCC1395_exome_normal.bam > $CNA_WORKSPACE/analysis/exome/hmmcopy/HCC1395_exome_normal.wig 
```

### Identify Heterozygous Germline Positions in the Normal and Extract Allele Counts in the Tumour

We use MutationSeq to identify heterozygous germline positions in the normal sample, and then extract allele counts at those positions from the tumour sample. 

Before we begin, we need to set the paths to our Python and MutationSeq installations, as well as to the index reference genome that our samples were aligned to:

~~~bash
PYTHON_DIR=$INSTALL_DIR/miniconda/miniconda2
MUSEQ_DIR=$INSTALL_DIR/mutationseq/mutationseq/museq
REF_DIR=$CNA_WORKSPACE/ref_data
~~~

We also need to export some paths for Python to run:

~~~bash
export PYTHONPATH=$PYTHON_DIR/lib/python2.7/site-packages/:$PYTHOPATH
export LD_LIBRARY_PATH=$PYTHON_DIR/lib:$LD_LIBRARY_PATH
~~~

The command below will run MutationSeq in single-sample mode on the normal to identify heterozygous germline positions, filter for high-confidence positions, and then extract the allele counts for those positions from the tumour sample:

~~~bash
mkdir -p $CNA_WORKSPACE/analysis/exome/mutationseq
cd $CNA_WORKSPACE/analysis/exome/mutationseq
$PYTHON_DIR/bin/python $MUSEQ_DIR/preprocess.py \
	reference:$REF_DIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
	model:$MUSEQ_DIR/models_anaconda/model_single_v4.0.2_anaconda_sk_0.13.1.npz \
	tumour:$CNA_WORKSPACE/data/HCC1395/exome/HCC1395_exome_tumour.bam \
	normal:$CNA_WORKSPACE/data/HCC1395/exome/HCC1395_exome_normal.bam \
	--verbose --single --coverage 4 --threshold 0.85 --buffer_size 2G \
	--out $CNA_WORKSPACE/analysis/exome/mutationseq/HCC1395_mutationseq.vcf \
	--log_file run_mutationseq.log \
	--config $MUSEQ_DIR/metadata.config \
	2> run_mutationseq.err &
~~~

You can monitor the run with `jobs` or by examining the log file with `less run_mutationseq.log`. Press `q` to exit the `less` program. This will take some time.

Once the MutationSeq run is complete, we use a custom script to convert the MutationSeq counts into the right format for input to Titan:

~~~bash
$PYTHON_DIR/bin/python $CNA_SCRIPTS_DIR/transform_vcf_to_counts.py \
	--infile $CNA_WORKSPACE/analysis/exome/mutationseq/HCC1395_mutationseq.vcf \
	--outfile $CNA_WORKSPACE/analysis/exome/mutationseq/HCC1395_mutationseq_postprocess.txt
~~~

You can also add an argument `--positions_file /path/to/file` with the path to a file listing known SNPs from the dbSNP database in the format `chr:position` with one SNP per line. This will filter the MutationSeq output to only include known positions from the dbSNP database. This is good practice but we will not do it in this case.

Finally, you will need to filter this file so that it contains only positions on chromosomes that are part of the Titan analysis (e.g. the autosomes and X). We called this filtered file `*_mutationseq_postprocess_filtered.txt`.

### Generate a Reference Genome Mappability File

We use the HMMCopy C package to generate the genome reference mappability file. This a two step process where we first need to generate a BigWig `.bw` file followed by converting it into a `.wig` file. Note that this will take some time.

```
$HMMCOPY_DIR/util/mappability/generateMap.pl \
	$CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
	-o $CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw
```

This may fail if you have not built a bowtie index of your reference yet. To do this, make sure first that bowtie is in your $PATH variable, then run the command above with the `-b` option:

```
$HMMCOPY_DIR/util/mappability/generateMap.pl \
	-b $CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
	-o $CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw
```

This will build the bowtie index. Now you can re-run the first command. In the end, you should get a BigWig file which you will need to convert into a wig file. This is done using the following command:

```
$HMMCOPY_DIR/bin/mapCounter \
	-w 1000 \
	$CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw \
	> $CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.ws_1000.wig
```

### Generate a Reference Genome GC Content File

We use the HMMCopy C package to generate the genome reference GC content file. This is done in a single step:

```
$HMMCOPY_DIR/bin/gcCounter \
	$CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
	> $CNA_WORKSPACE/ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.gc.wig
```

You are now ready to run Titan.
