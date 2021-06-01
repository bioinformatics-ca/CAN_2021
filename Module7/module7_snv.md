---
layout: tutorial_page
permalink: /BiCG_2019_Module7_Lab
title: BiCG 2019 Module 7 Lab
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: BiCG 2019 Module 7 Lab - Somatic Mutations and Annotations
author: Sorana Morrissy
modified: March 15th, 2018
---

# Lab Module 7 - Somatic Mutations and Annotation

This lab is designed to provide an introduction into Somatic Nucleotide Variation detection using two common programs: [Strelka](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts271) and [MuTect](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html). This lab will also go over simple annotation of the files using [Annovar](http://annovar.openbioinformatics.org/en/latest/) and manipulation of vcf files using [SnpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)

## Setup

First login into the server, and then enter the workspace directory:

```
cd ~/workspace
```

In order to keep our files in one location, we're going to create a new directory for this module and enter it:

```
mkdir Module7_snv
cd Module7_snv
```

Now to ease the use of commands that will follow, we're going to set some environment variables. These variables will stand in for the directories where our softwares of interest are location.

**Environment variables are only temporary. This means that once you log out of the server, these variables will be erased and must be reset if you'd like to use them again.**

Now we'll set out environment variables for MuTect and SnpEff (which will be used for processing the outputs from the callers):

```
JAVA_6_DIR=/home/ubuntu/CourseData/CG_data/Module7/install/java/jre1.6.0_45/bin
MUTECT_DIR=/home/ubuntu/CourseData/CG_data/Module7/install/mutect
SNPEFF_DIR=/usr/local/snpEff
ANNOVAR_DIR=/home/ubuntu/CourseData/CG_data/Module7/install/annovar
STRELKA_DIR=/home/ubuntu/CourseData/CG_data/Module7/install/strelka/bin/
```

## Linking the Sequencing and Referencing Data

For this lab module, we'll be using exome data on the HCC1395 breast cancer cell line. The tumour and normal bams have already been processed and placed on the server. So we'll create a soft link to the directory that it's stored. We'll also create a soft link to where the reference genome is stored, as well as a folder we'll use later on in the lab:

```
mkdir HCC1395 ref_data snv_analysis
ln -s /home/ubuntu/CourseData/CG_data/Module7/HCC1395/* HCC1395
ln -s /home/ubuntu/CourseData/CG_data/Module7/ref_data/* ref_data
ln -s /home/ubuntu/CourseData/CG_data/Module7/snv_analysis/* snv_analysis
```

For this lab we're going to limit our analysis to just the 7MB and 8MB region of chromosome 17 to ensure processing occurs quickly. The files we'll be focusing on can be viewed using the following command:

```
ls HCC1395/HCC1395_exome*.17*
```

You should see the files:  
* HCC1395_exome_normal.17.7MB-8MB.bam  
* HCC1395_exome_normal.17.7MB-8MB.bam.bai  
* HCC1395_exome_tumour.17.7MB-8MB.bam  
* HCC1395_exome_tumour.17.7MB-8MB.bam.bai  

We can take a look at the contents of our bam files by using samtools. The header contains information about the reference genome, commands used to generate the file, and on read groups (information about sample/flow cell/lane etc.), as well as the main contents of the file contain the one aligned read (read name, chromosome, position etc.) per line:

```
samtools view -h HCC1395/HCC1395_exome_normal.17.7MB-8MB.bam | less -s
```

# Predicting SNVs

## MuTect

In order to run MuTect, we simply run the java application of the program in one go. For our MuTect run, we're going to use two files to help filter out noise and focus on damaging mutations: common SNV's present in the general population as provided by dbSNP, and potentially damaging SNV's as documented in COSMIC. These files are stored in the same path as our MUTECT_DIR variable.

Looking at the dbsnp file:

```
less $MUTECT_DIR/dbsnp_132_b37.leftAligned.vcf
```

And now the cosmic file:

```
less $MUTECT_DIR/b37_cosmic_v54_120711.vcf
```

In order to run MuTect, we'll be running the following command:

```
mkdir results; 
mkdir results/mutect;
$JAVA_6_DIR/java -Xmx4g -jar $MUTECT_DIR/muTect-1.1.4.jar \
--analysis_type MuTect \
--reference_sequence ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.fa \
--cosmic $MUTECT_DIR/b37_cosmic_v54_120711.vcf \
--dbsnp $MUTECT_DIR/dbsnp_132_b37.leftAligned.vcf  \
--intervals 17:7000000-8000000 \
--input_file:normal HCC1395/HCC1395_exome_normal.17.7MB-8MB.bam \
--input_file:tumor HCC1395/HCC1395_exome_tumour.17.7MB-8MB.bam \
--vcf results/mutect/HCC1395.17.7MB-8MB_summary.vcf \
--out results/mutect/HCC1395.17.7MB-8MB_stats.out
```
*You can ignore the final warning:*  
`WARN  03:59:35,152 RestStorageService - Error Response: PUT '/GATK_Run_Reports/2PPjkmQ1fN3n7utCzrguPxiwxnT7XxwV.report.xml.gz' -- ResponseCode: 403, ResponseStatus: Forbidden...`

This is a failed attempt to report usage statistics + errors addressed in this [forum post](https://gatkforums.broadinstitute.org/gatk/discussion/1721/warnings-and-errors-during-score-calibration).

The `-Xmx4g` option allocates 4 GB of RAM for MuTect to run it's analysis. MuTect will output a vcf file with the calls.

The calls file contains information about the each SNV including information such as the reference allele, alternate allele and whether the SNV should be kept based on MuTect's quality filtering.

```
less -s results/mutect/HCC1395.17.7MB-8MB_summary.vcf
```

Now to subset for just the SNV's that have passed quality filtering, we use the following commands:

```
grep -v "REJECT" results/mutect/HCC1395.17.7MB-8MB_summary.vcf > results/mutect/mutect_passed.vcf
grep -v "REJECT" results/mutect/HCC1395.17.7MB-8MB_stats.out > results/mutect/mutect_passed.stats.out
```

To take a look into the end of files, we can use the command _tail_:

```
tail results/mutect/mutect_passed.vcf
```


## Strelka

Next, we'll call mutations using Strelka. Strelka provides aligner specific config files for bwa, eland, and isaac, each with preconfigured default parameters that work well with the aligner used. Since the bam files we're working with were processed using bwa, we'll make a local copy of that file:

```
mkdir config
cp /home/ubuntu/CourseData/CG_data/Module7/install/strelka/etc/strelka_config_bwa_default.ini config/strelka_config_bwa.ini
```

Our data is exome and so the coverage of the file is different, we need to change the `isSkipDepthFilters` parameter in the `strelka_config_bwa.ini` file. The default setting of `isSkipDepthFilters = 0` must simply be changed to `isSkipDepthFilters = 1`, and we'll accomplish this using the nano text editor:

```
nano config/strelka_config_bwa.ini
```
* Scroll to isSkipDepthFilters.
* Replace 0 with 1.
* *Ctrl+X* to exit.
* *Y* to save changes.
* *Enter* to maintain the same filename and overwrite.


OOOPS! That didn't work. Why?  
"Error writing config/strelka_config_bwa.ini: Permission denied"  

To see the file permissions:  
```
ls -lh config/strelka_config_bwa.ini
```
We only see the letter "r" meaning we can only read.  
We need write permission before we can edit a file.  

Give yourself write permissions using the command chmod.  
Then, check the permissions on this file again.  
```
chmod u+w config/strelka_config_bwa.ini
ls -lh config/strelka_config_bwa.ini
```
Now, edit this file with nano as specified above.  
  
    


The reason why we do this is described on the [Strelka FAQ page](https://sites.google.com/site/strelkasomaticvariantcaller/home/faq):

*The depth filter is designed to filter out all variants which are called above a multiple of the mean chromosome depth, the default configuration is set to filter variants with a depth greater than 3x the chromosomal mean. If you are using exome/targeted sequencing data, the depth filter should be turned off...*  
 
*However in whole exome sequencing data, the mean chromosomal depth will be extremely small, resulting in nearly all variants being (improperly) filtered out.*  

If you were doing this for whole genome sequencing data, then you should leave this parameter set to 0 as the depth of coverage won't be as high. 

A Strelka analysis is performed in 2 steps.  In the first step we provide Strelka with all the information it requires to run the analysis, including the tumour and normal bam filenames, the config, and the reference genome.  Strelka will create an output directory with the setup required to run the analysis.

```
$STRELKA_DIR/configureStrelkaWorkflow.pl \
    --tumor HCC1395/HCC1395_exome_tumour.17.7MB-8MB.bam \
    --normal HCC1395/HCC1395_exome_normal.17.7MB-8MB.bam \
    --ref ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.fa \
    --config config/strelka_config_bwa.ini \
    --output-dir results/strelka/
```

The output directory will contain a _makefile_ that can be used with the tool _make_.

To run the Strelka analysis, use make and specify the directory constructed by `configureStrelkaWorkflow.pl` with make's `'-C'` option.

```
make -C results/strelka/ -j 2
```

The `-j 2` parameter specifies that we want to use 2 cores to run Strelka. Change this number to increase or decrease the parallelization of the job. The more cores the faster the job will be, but the higher the load on the machine that is running Strelka. 

*If you have access to a grid engine cluster, you can replace the command `make` with `qmake` to launch Strelka on the cluster.*  

Strelka has the benefit of calling SNVs and small indels.  Additionally, Strelka calculates variant quality and filters the data in two tiers.  The filenames starting with `passed` contain high quality candidates, and filenames starting with `all` contain high quality and marginal quality candidates.

The Strelka results are in VCF format, with additional fields explained on the [strelka website](https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output).

```
less -s results/strelka/results/passed.somatic.snvs.vcf
less -s results/strelka/results/passed.somatic.indels.vcf
```

## Annotating our mutations

Now that we have our list of mutations, we can annotate them with multiple layers of information, such as the gene name, the function of the location (intronic, exonic), whether the mutation belongs to dbSNP or the 1000 genome project, and multiple other layers.

Let's place these results in a new folder within our results folder:

```
mkdir results/annotated/
```

To run Annovar on the MuTect output, we use the following command:

```
$ANNOVAR_DIR/table_annovar.pl results/mutect/mutect_passed.vcf $ANNOVAR_DIR/humandb/ -buildver hg19 -out results/annotated/mutect -remove -protocol refGene,cytoBand,genomicSuperDups,1000g2015aug_all,avsnp147,dbnsfp30a -operation g,r,r,f,f,f -nastring . --vcfinput
```

We can view the resulting file by going _less_ on it:

```
less results/annotated/mutect.hg19_multianno.vcf
```

Unfortunately our strelka vcf doesn't have the mandatory info field for genotype _GT_ in the format that annovar requires, and so some manipulation of our strelka results are needed. Annovar requires a minimum of the chromosome, start position, end position, reference allele, and alternate allele.

The following command will remove the header information, format the file as needed by annovar while retaining other information, and storing it into our new file _passed.somatic.snvs.txt_.

```
grep -v "##" results/strelka/results/passed.somatic.snvs.vcf | awk '{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > results/strelka/results/passed.somatic.snvs.txt
```

Now we can annotate the Strelka output as we did before:

```
$ANNOVAR_DIR/table_annovar.pl results/strelka/results/passed.somatic.snvs.txt $ANNOVAR_DIR/humandb/ -buildver hg19 -out results/annotated/strelka -remove -protocol refGene,cytoBand,genomicSuperDups,1000g2015aug_all,avsnp147,dbnsfp30a -operation g,r,r,f,f,f -nastring . -csvout
```

The output from this is stored in comma separated format as results/annotated/strelka.hg19_multianno.csv which can be downloaded and viewed locally in excel.


## Parsing specific fields from our vcf files

The VCF format is sometimes not useful for visualization and data exploration purposes which often requires the data to be in tabular format. We can convert from VCF format to tabular format using the extractField() function from SnpSift/SnpEff. Since each mutation caller has a different set of output values in the VCF file, the command needs be adjusted for each mutation caller depending on the fields present in the header.

From the [SnpSift documentation](http://snpeff.sourceforge.net/SnpSift.html#filter):  
All VCF fields can be used as variables names, as long as they are declared in the VCF header OR they are "standard" VCF fields.

For our example:  
Standard fields: `CHROM POS REF ALT` - Defined by VCF 4.1 specification  
Info field variables: `Func.refGene Gene.refGene cytoBand ExonicFunc.refGene avsnp147` -  Defined in header lines ##INFO=<ID=...  
Genotype fields: `GEN[0].FA GEN[1].FA` - Defined in header lines ##FORMAT=<ID=FA... where:  
* GEN[0] = first sample = sample_HCC1395_tumour
* GEN[1] = second sample = sample_HCC1395_normal

Therefore, to convert the MuTect VCF file into a tabular format containing only the chromosome, position, reference allele, alternate allele, function of gene, gene name, cytoband, exonic function, avsnp147 annotation, allele frequency of normal, and allele frequency of control:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "." results/annotated/mutect.hg19_multianno.vcf CHROM POS REF ALT Func.refGene Gene.refGene cytoBand ExonicFunc.refGene avsnp147 GEN[0].FA GEN[1].FA > results/mutect/mutect_passed_tabbed.txt
```

We can view the new tab delimited file in our browser at http://XX.oicrcbw.ca/Module7_snv/results/mutect (click mutect_passed_tabbed.txt)

To convert the Strelka VCF file into a tabular format into a similar format containing the chromosome, position, reference allele, alternate allele, depth of normal, depth of tumour, # of A's, C's, G's, and T's in normal and tumour (two values each representing "stringent,permissive" filter tiers):

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "."  results/strelka/results/passed.somatic.snvs.vcf CHROM POS REF ALT GEN[0].DP GEN[1].DP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU > results/strelka/strelka_snvs_tabbed.txt 
```

The -e parameter specifies how to represent empty fields. In this case, the "." character is placed for any empty fields. This facilitates loading and completeness of data. For more details on the extractField() function see the [SnpSift documentation](http://snpeff.sourceforge.net/SnpSift.html#Extract).

We can view the first 10 lines of our tabbed file by running the following commands (one at a time):

```
head results/mutect/mutect_passed_tabbed.txt
head results/strelka/strelka_snvs_tabbed.txt
```
## Data Exploration

### Integrative Genomics Viewer (IGV)

A common step after prediction of SNVs is to visualize these mutations in IGV. Let's load these bam into IGV. Open IGV, then:

1. Change the genome to hg19 (if it isn't already)  
2. File -> Load from URL ...  
    * http://xx.oicrcbw.ca/Module7_snv/HCC1395/HCC1395_exome_tumour.17.7MB-8MB.bam
    * http://xx.oicrcbw.ca/Module7_snv/HCC1395/HCC1395_exome_normal.17.7MB-8MB.bam


Where the xx is your student number. Once the tumour and normal bam have been loaded into IGV, we can investigate a few predicted positions in IGV:

* 17:7491818  
* 17:7578406  
* 17:7593160  
* 17:7011723  
* 17:7482929  

Manually inspecting these predicted SNVs in IGV is a good way to verify the predictions and also identify potential false positives:

*When possible, you should always inspect SNVs*

### Exploration in R

While IGV is good for visualizing individual mutations, looking at more global characteristics would require loading the data into an analysis language like R.  

We will use exome-wide SNV predictions for MuTect for these analyses; specifically, we're only going to look at the _stats.out_ output from MuTect that has been run on the whole exome file. The processed tabular text files along with the `snv_analysis.Rmd` RMarkdown file that contains the R code is available in our `snv_analysis` folder.

Now let's launch our RStudio instance and continue our analysis there.  

In a web browser, navigate to:  
http://xx.oicrcbw.ca:8080/  
Replacing xx with your student number.

Within RStudio, click File > Open File > workspace/Module7_snv/snv_analysis/snv_analysis.Rmd


(Alternatively, open a new RMarkdown file and paste the contents of the file [here](https://raw.githubusercontent.com/bioinformaticsdotca/BiCG_2018/master/module7_snv/snv_analysis.txt).)    
***
