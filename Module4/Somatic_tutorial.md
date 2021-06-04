---
layout: tutorial_page
permalink: /CAN_2021_module3_lab
title: CAN 2021 Module 4 Lab
header1: Workshop Pages for Students
header2: Cancer Analyis 2021 Module 4
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/CAN_2021
description: CAN 2021 Module 3 lab
author: Aaron Gillmor
 
---


================================

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

================================


In this workshop, we will work with common tools to process and analyze cancer sequencing data. Using the command line we will analyze DNA sequencing from the whole genome. The focus of this workshop will be on single nucelotide variants, insertion, deletions (commonly referred to as SNV/indels) as well as calling copy-number variations (CNVs). We will also annotate SNV & CNV files, such that you will known the functional conseqeunce of these variants.


Data Source
We will be working on the CageKid samples from Module3, specifically --> patient C0098. 
Whole genome sequencing and analysis can take multiple days to run, as such we have downsampled the files so that we can proceed more quickly. For the SNV analsyis we have selected the region on chromsome 9 between 130215000 & 130636000. 

The tools we are going to using for assessing SNV's is mutect2 and varscan2 (data and provide command lines that allow detecting Single Nucleotide Variants (SNV).


**Running Mutect2**

```
cd ~/workspace
mkdir Module4_somaticvariants
cd Module4_somaticvariants
```


```
/usr/local/GATK/gatk Mutect2 -R /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -I /home/ubuntu/CourseData/CAN_data/Module4/alignments/normal/normal.sorted.realigned.bam -I /home/ubuntu/CourseData/CAN_data/Module4/alignments/tumor/tumor.sorted.realigned.bam -normal normal -tumor tumor -O pairedVariants_mutect2.vcf -L 9:130215000-130636000
```

#This will generate an intial vcf but does not contain any filters which tell us important information about the variants.

```
gatk FilterMutectCalls -R /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta --filtering-stats pairedVariants_mutect2.vcf.stats -V pairedVariants_mutect2.vcf -O pairedVariants_mutect2_filtered.vcf
```

Let's look at our file now

```
less -S pairedVariants_mutect2_filtered.vcf
```

#Although we have all of the information to proceed forward we still have multiple variants at a single loci. And since we will be comparing to another tool it is best to consider these multiallelic variants which can be accomplished by using bcftools norm functionality.
#Next I want to split multiallelic (example Ref == A & Alt == AT,ATT)

```
bcftools norm -m-both -f /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -Oz -o pairedVariants_mutect2_filtered_normalized.vcf pairedVariants_mutect2_filtered.vcf
```

#Running varscan2 (http://varscan.sourceforge.net/)

First we create pileup files from both the normal sample and the blood sampled. This is done with a useful tool called samtools:
```
samtools mpileup -B -q 1 -f /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -r 9:130215000-130636000 /home/ubuntu/CourseData/CAN_data/Module4/alignments/normal/normal.sorted.realigned.bam > normal.mpileup &&
```
```
samtools mpileup -B -q 1 -f /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -r 9:130215000-130636000 /home/ubuntu/CourseData/CAN_data/Module4/alignments/normal/normal.sorted.realigned.bam > normal.mpileup &&
```
#Next we use varscan to execute the somatic function which does the initial variant calling from our mpileups

```
java -Xmx2G -jar /usr/local/VarScan.v2.3.9.jar somatic normal.mpileup tumor.mpileup --output-vcf 1 --strand-filter 1 --somatic-p-value 0.05 --output-snp varscan2.snp.vcf --output-indel varscan2.indel.vcf
```

#Running we can run the intersection and union using bcftools. 
```
varscan processSomatic varscan2.indel.vcf --min-tumor-freq 0.05
varscan processSomatic varscan2.snp.vcf --min-tumor-freq 0.05
```

#zip and tabix of high confidence variants
```
bgzip varscan2.snp.Somatic.hc.vcf
tabix varscan2.snp.Somatic.hc.vcf.gz

bgzip varscan2.indel.Somatic.hc.vcf
tabix varscan2.indel.Somatic.hc.vcf.gz
```

#Combine high-confidence snp and indels
```
bcftools concat -Ov -a -r 9:130215000-130636000 varscan2.snp.Somatic.hc.vcf.gz varscan2.indel.Somatic.hc.vcf.gz -o -o pairedVariants_varscan2_filtered.vcf
```
```
bgzip  pairedVariants_varscan2_filtered.vcf
tabix  pairedVariants_varscan2_filtered.vcf.gz
```

#Now we have finished variant calling using two different variants callers.

##combining with multiple callers
```
bcftools isec -n+1 -f PASS -p intersectionDirectory pairedVariants_mutect2_filtered_normalized.vcf.gz pairedVariants_varscan2_filtered.vcf.gz -o pairedVariants_mutect2_varscan2.vcf -Ov
```
This creates a directory named intersectionDirectory that contains all variants that pass in at-least one vcf. It's at this point where you can choose which variants to include.

```
cd intersectionDirectory/
```
- All variants -- The union will represent the most variants but has the higher false positives count
- Common variants -- The intersection will represerent the least variants but has a lower false positive 

In this directory we see the 0000.vcf and 0001.vcf. These corresponds to the mutect2 vcf and the varscan2 vcf, respectively. (confirmed in the README.txt)

If we review the sites.txt files:

less -S sites.txt | head -n 4

    9       130223126       C       T       01
    9       130264057       T       A       11
    9       130296899       A       T       11
    9       130308743       C       A       11

This reads that the first variants on chromosme 9, position 130223126 changes from a C to a T: and that this variant was called in the varscan.vcf (01). 11 corresponds to the being called in both and 10 is being called in only mutect2. 

less -S sites.txt | cut -f  5 | sort | uniq -c 
      
      4 01
      6 10
     11 11

This show's us that varscan has 4 unique variants, mutect2 has 6 unique variants and there are 11 shared variants.

#Finally we combine our mutect2 and varscan2 into a single vcf.
```
bcftools merge pairedVariants_mutect2_filtered_normalized.vcf.gz pairedVariants_varscan2_filtered.vcf.gz -o pairedVariants_mutect2_varscan2.vcf.gz -Oz
```
#Now that we have our two vcf's we can look at their contents. A full description of the [vcf_format](https://www.internationalgenome.org/wiki/Analysis/vcf4.0/):

## Every vcf has a metadata sections two ##
```
zless -S pairedVariants_varscan2_filtered.vcf | egrep ' ##' | head -n 10

```
      ##fileformat=VCFv4.2
      ##FILTER=<ID=PASS,Description="All filters passed">
      ##FILTER=<ID=FAIL,Description="Fail the site if all alleles fail but for different reasons.">
      ##FILTER=<ID=base_qual,Description="alt median base quality">
      ##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">
      ##FILTER=<ID=contamination,Description="contamination">
      ##FILTER=<ID=duplicate,Description="evidence for alt allele is overrepresented by apparent duplicates">
      ##FILTER=<ID=fragment,Description="abs(ref - alt) median fragment length">
      ##FILTER=<ID=germline,Description="Evidence indicates this site is germline, not somatic">
      ##FILTER=<ID=haplotype,Description="Variant near filtered variant on same haplotype.">


This section contains details about the vcf including 
- explanations of the filters, format and ID's. if you want to know what something means in a vcf: this is where you look! <question >
- contigs used and their length
- Processes taken to generate the vcf eg bcftools_mergeCommand=merge -o pairedVariants_mutect2_varscan2.vcf.gz -Oz 0000.vcf.gz 0001.vcf.gz

Every vcf has a column (with one #) and;
```
less -S pairedVariants_varscan2_filtered.vcf | egrep '#CHROM' | head -n 1
```
```
CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  normal  tumor   NORMAL  TUMOR 
```
 
The header has 8 mandatory columns. CHROM, POS, ID, REF, ALT, QUAL, FILTER & INFO
These are followed by a FORMAT column header and sample ID's. These sample id's will follow from their vcf's 

# Every vcf content section which has the variants.

```
less -S pairedVariants_varscan2_filtered.vcf | egrep -v '#' | head -n 1
 
```
 
```
9       130223126       .       C       T       .       PASS    SOMATIC;SS=2;SSC=18;GPV=1;SPV=0.014587;DP=55    GT:GQ:DP:RD:AD:FREQ:DP4 ./.:.:.:.:.:.:. ./.:.:.:.:.:.:. 0/0:.:26:26:0:0%:14,12,0,0      0/1:.:29:22:6:21.43%:15,7,1,5
 ```
 
 This shows a variant on chromosome 9 position 130223126 that goes from the reference C to a T. This variant was called uniquely by varscan2. This variant has passed all filters. The INFO and FORMAT field has additional information that can be explained in more detail using the '##' information. 


#Annotating our combined files
```
mkdir results
```

#annovar 
```
/usr/local/annovar/table_annovar.pl pairedVariants_mutect2_varscan2.vcf /usr/local/annovar/humandb/ -buildver hg19 -out results/annotated_mutect2_varscan2 -remove -protocol refGene -operation g -nastring . --vcfinput
```
This will produce two annotation files: annotated_mutect2_varscan2.hg19_multianno.vcf & annotated_mutect2_varscan2.hg19_multianno.txt. These files contain the gene annotations and infer the functional consequence of each variant. 

```
 less -S annotated_mutect2_varscan2.hg19_multianno.txt | cut -f 6 | sort | uniq -c
```
 This shows the genomic region change where each variant is occuring. Typically exonic variants are more likely to produce a phenotypic chance. 

        5 exonic
        1 intergenic
        15 intronic

 
```
 less -S annotated_mutect2_varscan2.hg19_multianno.txt | cut -f 9 | sort | uniq -c
```

        16 .
        2 frameshift insertion
        3 nonsynonymous SNV

 This shows the functional consequence of each variant, we have 5 exonic variants that results in two frameshifts and 3 nonsynonymous SNV's.
 
 
Next we can examine some variants in IGV. 
 
 start IGV 
 
 load your bam files by clicking file --> load from file 
 
 input the variant region of intereset. Input chr9:130,634,091-130 which corresponds to a exonic mutation called in AK1.
 ![image](https://user-images.githubusercontent.com/15352153/120568604-72690b00-c3d1-11eb-8164-329868b51e7c.png)

 Here we can see evidence of the C --> T in the tumor bam but not in the normal sample. This variant was called correctly. 
 
 
![image](https://user-images.githubusercontent.com/15352153/120569761-d42a7480-c3d3-11eb-8897-12cbfa600a65.png)

  
    
    

 

 
