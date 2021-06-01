---
layout: tutorial_page
permalink: /BiCG_2019_Module5_lab_supplemental
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: BiCG Module 5 Lab Supplemental
author: Jared Simpson
modified: June 27th, 2018
---

## Generating a preqc report for the E. coli data set

```
# First build an FM-index of the two E. coli Illumina data sets:
sga index -t 4 -a ropebwt ecoli.illumina.15x.fastq
sga index -t 4 -a ropebwt ecoli.illumina.50x.fastq

# Next, we can run `sga preqc` to run the calculations and generate the PDF report.
sga preqc -t 4 ecoli.illumina.15x.fastq > ecoli.illumina.15x.preqc
sga preqc -t 4 ecoli.illumina.50x.fastq > ecoli.illumina.50x.preqc
python sga-preqc-report.py ecoli.illumina.15x.preqc ecoli.illumina.50x.preqc
```
