---
layout: tutorial_page
permalink: /BiCG_2019_Module6_Lab-Installation
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
---

## Introduction

This page describes where to get and how to install the software used for this lab module. The table of contents for this page is:

* [Installing Array Data Preprocessing Software](#installing-array-data-preprocessing-software)
* [Installing OncoSNP](#installing-oncosnp)
* [Installing Sequencing Data Preprocessing Software](#installing-sequencing-data-preprocessing-software)
* [Installing HMMcopy](#installing-hmmcopy)
* [Installing Titan](#installing-titan)

To start, create an environment variable specifying a path to your install directory:
~~~bash
INSTALL_DIR=/path/to/install_dir
~~~

## Installing Array Data Preprocessing Software

We will use the procedure described on the [PennCNV site](http://www.openbioinformatics.org/penncnv/penncnv_tutorial_affy_gw6.html). This page also contains links to most of the software we need to download.

Note: You will need to register with Affymetrix to download the library files, specifically the `.cdf` file, from the [Affimetrix site](http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6). The [Aroma Project](http://www.aroma-project.org) has a copy available which we will download since it can be done from the cloud.

Download PennCNV, which we will use to extract BAF and LRR measurements:
~~~bash
mkdir -p $INSTALL_DIR/penncnv
cd $INSTALL_DIR/penncnv
wget https://github.com/WGLab/PennCNV/archive/v1.0.3.tar.gz
tar -xzvf v1.0.3.tar.gz
~~~

Download the PennCNV package for Affymetrix Genome-Wide Human SNP Array 6.0:
~~~bash
cd $INSTALL_DIR/penncnv
wget http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz
tar -xzvf gw6.tar.gz
~~~

Add the hg19 version of the `.pfb` file, since the package above contains only the hg18 version:
~~~bash
cd $INSTALL_DIR/penncnv/gw6/lib
wget http://www.openbioinformatics.org/penncnv/download/affygw6.hg19.pfb.gz
gunzip affygw6.hg19.pfb.gz
~~~

Download Affymetrix Power Tools (APT), which we will use to normalize the array data:
~~~bash
mkdir -p $INSTALL_DIR/apt
cd $INSTALL_DIR/apt
wget https://downloads.thermofisher.com/Affymetrix_Softwares/APT_1.19.0_Linux_64_bit_x86_binaries.zip
unzip APT_1.19.0_Linux_64_bit_x86_binaries.zip
cd $INSTALL_DIR/apt/apt-1.19.0-x86_64-intel-linux/bin
chmod +x apt*
~~~

Download the SNP6 cel definition file:
~~~bash
cd $INSTALL_DIR/apt
wget http://www.aroma-project.org/data/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf.gz
gunzip GenomeWideSNP_6.cdf.gz
~~~

## Installing OncoSNP

We will use [OncoSNP](https://sites.google.com/site/oncosnp/) to infer CNVs from the processed array data. The OncoSNP files can be downloaded [here](https://sites.google.com/site/oncosnp/user-guide/downloads). Note that you will need to register with the author and he will supply a password to unlock the downloaded files.

Here we will download OncoSNP version 2.0:
~~~bash
mkdir -p $INSTALL_DIR/oncosnp
cd $INSTALL_DIR/oncosnp
wget http://www.well.ox.ac.uk/~cyau/oncosnp/oncosnp_v2.0.run
bash oncosnp_v2.0.run          # enter your oncosnp password
~~~

Then follow the on-screen instructions. 

You'll also need to install the MATLAB MCR to run OncoSNP:

~~~bash
mkdir -p $INSTALL_DIR/matlab_mcr
cd $INSTALL_DIR/matlab_mcr
wget http://www.mathworks.co.uk/supportfiles/downloads/R2013b/deployment_files/R2013b/installers/glnxa64/MCR_R2013b_glnxa64_installer.zip
unzip MCR_R2013b_glnxa64_installer.zip
echo destinationFolder=$INSTALL_DIR/matlab_mcr > $INSTALL_DIR/matlab_mcr/input_file.txt
echo agreeToLicense=yes >> $INSTALL_DIR/matlab_mcr/input_file.txt
echo mode=silent >> $INSTALL_DIR/matlab_mcr/input_file.txt
./install -mode silent -inputFile $INSTALL_DIR/matlab_mcr/input_file.txt

~~~

Make sure that the LD_LIBRARY_PATH contains the path to the Matlab MCR, with the correct version vXX (replace XX in the command below with the installed version number):

~~~bash
export LD_LIBRARY_PATH="$MCR_DIR/vXX:$LD_LIBRARY_PATH"
~~~

Finally, you will need to download GC content files for the relevant build of the genome. For this tutorial we will use the hg19 (build GRCh37) files:

~~~bash
mkdir -p $INSTALL_DIR/gc_content
cd $INSTALL_DIR/gc_content
wget http://www.well.ox.ac.uk/~cyau/gc/b37.zip
unzip b37.zip
~~~

## Installing Sequencing Data Preprocessing Software

We will use [MutationSeq](https://bitbucket.org/shahlabbcca/mutationseq) to identify heterozygous germline variants for input to Titan. To install MutationSeq run:
~~~bash
mkdir -p $INSTALL_DIR/mutationseq
cd $INSTALL_DIR/mutationseq
wget ftp://ftp.bcgsc.ca/public/shahlab/MutationSeq/museq_4.3.8_with_anaconda_models.tar.gz
tar -zxvf museq_4.3.8_with_anaconda_models.tar.gz
~~~

Mutationseq requires Python 2.7 and several key packages to run (numpy, scipy, matplotlib, scikit-learn, and intervaltree). We will install the [Minoconda distribution](https://conda.io/miniconda.html):

~~~bash
mkdir -p $INSTALL_DIR/miniconda
cd $INSTALL_DIR/miniconda
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
~~~

Follow the on-screen instructions, and specify the installation directory ($INSTALL_DIR/miniconda/miniconda2). 

Now use conda to install the MutationSeq dependencies:
~~~bash
$INSTALL_DIR/miniconda/miniconda2/bin/conda install -c bioconda numpy=1.7.1 scipy=0.12.0 scikit-learn=0.13.1 matplotlib=1.2.1 intervaltree
~~~

We also need to install the Boost C libraries (note the version is important, later versions may cause errors):
~~~bash
mkdir -p $INSTALL_DIR/boost
cd $INSTALL_DIR/boost
wget http://sourceforge.net/projects/boost/files/boost/1.51.0/boost_1_51_0.tar.gz
tar -zxvf boost_1_51_0.tar.gz
~~~

Finally, we compile a dependency `pybam.so` using the makefile in the MutationSeq directory:
~~~bash
cd $INSTALL_DIR/mutationseq/mutationseq/museq
make PYTHON=$INSTALL_DIR/miniconda/miniconda2/bin/python BOOSTPATH=$INSTALL_DIR/boost/boost_1_51_0 -B
~~~

MutationSeq is installed correctly if you run the following and get the version number printed to the terminal (e.g. 4.3.8):
~~~bash
$INSTALL_DIR/miniconda/miniconda2/bin/python $INSTALL_DIR/mutationseq/mutationseq/museq/classify.py --version
~~~

## Installing HMMcopy

HMMcopy has two components, a Bioconductor R package and a C package. We will use the HMMcopy C package to extract binned read counts from our tumour and normal samples, and use this as input to Titan. Note that the HMMcopy R package can be used to call CNAs based on read depth, but it does not make use of b-allele frequencies, and does not model ploidy, normal cell contamination, and tumour heterogeneity. For these reasons, we will not use it for our bulk tumour analysis.

Information about the Bioconductor package can be found [here](http://bioconductor.org/packages/release/bioc/html/HMMcopy.html).

To install the HMMcopy R package, first start R:

~~~bash
R
~~~

In the R environment, enter the following commands:

~~~
source("http://bioconductor.org/biocLite.R")
biocLite("HMMcopy")
~~~

Follow along and accept any suggested updates.

The HMMcopy C package can be downloaded from the BCCRC:

~~~bash
mkdir -p $INSTALL_DIR/hmmcopy
cd $INSTALL_DIR/hmmcopy
wget http://compbio-bccrc.sites.olt.ubc.ca/files/2013/12/HMMcopy.zip
unzip HMMcopy.zip
cd HMMcopy
cmake .
make
~~~

## Installing Titan

Titan is an R package available through Bioconductor. To install it, first start R:

~~~bash
R
~~~

In the R environment, enter the following commands:

~~~
source("http://bioconductor.org/biocLite.R")
biocLite("TitanCNA")
~~~
