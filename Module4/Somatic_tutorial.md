---
layout: tutorial_page
permalink: /CAN_2021_module3_lab
title: CAN 2021 Module 4 Lab
header1: Workshop Pages for Students
header2: Cancer Analyis 2021 Module 4
image: 
home: https://bioinformaticsdotca.github.io/CAN_2021
description: CAN 2021 Module 4 lab
author: Aaron Gillmor
 
---

================================

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

================================


In this workshop, we will work with common tools to process and analyze cancer sequencing data. Using the command line we will analyze DNA sequencing from the whole genome. The focus of this workshop will be on calling single nucelotide variants, insertion, deletions (commonly referred to as SNV/indels) as well as calling copy-number variations (CNVs). We will also annotate SNV & CNV files, such that you will known the functional conseqeunce of these variants.


Data Source:
we will be working on the [CageKid](https://www.cnrgh.fr/cagekid/) samples from Module3. 
Whole genome sequencing and analysis can take multiple days to run, as such we have downsampled the files so that we can proceed more quickly. For the SNV analsyis we have selected a region on chromsome 9 between 130215000 & 130636000. 

The tools and their general function we are going to using for calling SNV's and CNV's are:
* mutect2
* varscan2
* samtools
* bcftools
* bgzip and tabix
* annovar
* controlfreec
* R & R-Studio

Files for variant calling:
1) Bam --> Sequence alignments from Module 3 or /home/ubuntu/CourseData/CAN_data/Module4/alignments/normal/normal.sorted.realigned.bam
2) Reference genome --> /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta
3) Germline-reference file --> /home/ubuntu/CourseData/CAN_data/Module4/accessory_files/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz

**Running Mutect2**

```
cd ~/workspace
mkdir -p Module4_somaticvariants
cd Module4_somaticvariants
```

```
/usr/local/GATK/gatk Mutect2 -R /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -I /home/ubuntu/workspace/CBW_CAN_2021/Module4/alignment/normal/normal.sorted.dup.recal.bam -I /home/ubuntu/workspace/CBW_CAN_2021/Module4/alignment/tumor/tumor.sorted.dup.recal.bam -normal normal -tumor tumor -O pairedVariants_mutect2.vcf -L 9:130215000-130636000 --germline-resource /home/ubuntu/CourseData/CAN_data/Module4/accessory_files/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz
```


This will generate an intial vcf but does not contain any filters which tell us important information about the variants.

```
less -S pairedVariants_mutect2.vcf
```

```
gatk FilterMutectCalls -R /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta --filtering-stats pairedVariants_mutect2.vcf.stats -V pairedVariants_mutect2.vcf -O pairedVariants_mutect2_filtered.vcf
```

Let's look at our file now

```
less -S pairedVariants_mutect2_filtered.vcf
```

Although we have all of the information to proceed forward we have multiple variants at a single loci. And since we will be comparing to another tool it is best to consider these multiallelic variants which can be accomplished by using bcftools norm functionality. This is to say: I want to split multiallelic (Ref == A & Alt == AT,ATT) into two seperate variants A --> AT and A --> ATT. This is because some variants callers (there are many) will report the multiallelic (mutect2) and some will report the singleallele variant (varscan2).

```
bcftools norm -m-both -f /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -Oz -o pairedVariants_mutect2_filtered_normalized.vcf pairedVariants_mutect2_filtered.vcf
```

Running varscan2 (http://varscan.sourceforge.net/)

First we create pileup files from both the normal sample and the control sample. This is done with a tool called samtools:
```
samtools mpileup -B -q 1 -f /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -r 9:130215000-130636000 /home/ubuntu/CourseData/CAN_data/Module4/alignments/normal/normal.sorted.realigned.bam > normal.mpileup &&
```
```
samtools mpileup -B -q 1 -f /home/ubuntu/CourseData/CAN_data/Module4/references/human_g1k_v37.fasta -r 9:130215000-130636000 /home/ubuntu/CourseData/CAN_data/Module4/alignments/normal/normal.sorted.realigned.bam > normal.mpileup &&
```
Next we use varscan2 to execute the somatic function which does the initial variant calling from our mpileups

```
java -Xmx2G -jar /usr/local/VarScan.v2.3.9.jar somatic normal.mpileup tumor.mpileup --output-vcf 1 --strand-filter 1 --somatic-p-value 0.05 --output-snp varscan2.snp.vcf --output-indel varscan2.indel.vcf
```

#Running we can run the intersection and union using bcftools. 
```
varscan processSomatic varscan2.indel.vcf --min-tumor-freq 0.05
varscan processSomatic varscan2.snp.vcf --min-tumor-freq 0.05
```

zip and tabix of high confidence variants
```
bgzip varscan2.snp.Somatic.hc.vcf
tabix varscan2.snp.Somatic.hc.vcf.gz

bgzip varscan2.indel.Somatic.hc.vcf
tabix varscan2.indel.Somatic.hc.vcf.gz
```

Combine high-confidence snp and indels
```
bcftools concat -Ov -a -r 9:130215000-130636000 varscan2.snp.Somatic.hc.vcf.gz varscan2.indel.Somatic.hc.vcf.gz -o -o pairedVariants_varscan2_filtered.vcf
```
```
bgzip  pairedVariants_varscan2_filtered.vcf
tabix  pairedVariants_varscan2_filtered.vcf.gz
```

Now we have finished variant calling using two different variants callers. The next step is to combine the output of the multiple callers
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
      
      3 01
      6 10
     11 11

This show's us that varscan has 4 unique variants, mutect2 has 6 unique variants and there are 11 shared variants.

Finally we combine our mutect2 and varscan2 vcfs into a single vcf. This will generate the union of variants. 
```
cd ../
```

```
bcftools merge pairedVariants_mutect2_filtered_normalized.vcf.gz pairedVariants_varscan2_filtered.vcf.gz -o pairedVariants_mutect2_varscan2.vcf.gz -Oz
```
Now that we have our two vcf's combined we can look at their contents. A full description of the [vcf_format](https://www.internationalgenome.org/wiki/Analysis/vcf4.0/):
The generate fields will depend on the caller's used. common ones are:
* Reference allele (REF) and Alternate allele (ALT)
* Variant quality (QUAL)
* Genotype per-sample (GT) note this is not in Strelka which is a popular variant caller


Every vcf has a metadata sections two ##
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


This section contains details about the vcf including:
- Explanations of the FILTERS, INFO and FORMAT.
- Contigs used and their length
- Processes taken to generate the vcf eg bcftools_mergeCommand=merge -o pairedVariants_mutect2_varscan2.vcf.gz -Oz 0000.vcf.gz 0001.vcf.gz

Every vcf has a column (with one #) and
 
```
less -S pairedVariants_varscan2_filtered.vcf | egrep '#CHROM' | head -n 1
 
```
```
CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  normal  tumor   NORMAL  TUMOR 
```
 
The header has 8 mandatory columns. CHROM, POS, ID, REF, ALT, QUAL, FILTER & INFO
These are followed by a FORMAT column header and sample ID's. These sample id's will follow from their vcf's 

Every vcf content section which has the variants.

```
less -S pairedVariants_varscan2_filtered.vcf | egrep -v '#' | head -n 1
 
```
 
```
9       130223126       .       C       T       .       PASS    SOMATIC;SS=2;SSC=18;GPV=1;SPV=0.014587;DP=55    GT:GQ:DP:RD:AD:FREQ:DP4 ./.:.:.:.:.:.:. ./.:.:.:.:.:.:. 0/0:.:26:26:0:0%:14,12,0,0      0/1:.:29:22:6:21.43%:15,7,1,5
 ```
 
This shows a variant on chromosome 9 position 130223126 that goes from the reference C to a T. This variant was called uniquely by varscan2. This variant has passed all filters. The INFO and FORMAT field has additional information that can be explained in more detail using the '##' information. Please not that missing data is denoted as a '.'.


Annotating our combined files
```
mkdir -p results
```

annovar is a tool to annotate our vcf files with gene names & functional information 
```
/usr/local/annovar/table_annovar.pl pairedVariants_mutect2_varscan2.vcf.gz /usr/local/annovar/humandb/ -buildver hg19 -out results/annotated_mutect2_varscan2 -remove -protocol refGene -operation g -nastring . --vcfinput
```
This will produce two annotation files: annotated_mutect2_varscan2.hg19_multianno.vcf & annotated_mutect2_varscan2.hg19_multianno.txt. These files contain the gene annotations and infer the functional consequence of each variant. 

```
 less -S annotated_mutect2_varscan2.hg19_multianno.txt | cut -f 6 | sort | uniq -c
```
 This shows the genomic region change where each variant is occuring. Typically exonic variants are more likely to produce a phenotypic change. 
 
        5 exonic
        15 intronic
        
```
 less -S annotated_mutect2_varscan2.hg19_multianno.txt | cut -f 9 | sort | uniq -c
```

        15 .
        2 frameshift insertion
        3 nonsynonymous SNV

This shows the functional consequence of each exonic variant, we have 5 exonic variants that results in two frameshifts and 3 nonsynonymous SNV's. 

 
Next we can examine some variants in IGV. 
 
 start IGV 
 ensure hg19 is selected
 
 load your bam files by clicking file --> load from file 
 If you have not downloaded the bam or are missing it. It is provided here (https://drive.google.com/drive/folders/1f_1pDbpNaNUT6oC-pEPXInemyZwMV4Hc)
 
Input chr9:130,634,091.
![image](https://user-images.githubusercontent.com/15352153/120568604-72690b00-c3d1-11eb-8164-329868b51e7c.png)

![image](https://user-images.githubusercontent.com/15352153/120906146-438ea700-c614-11eb-9fb1-6dc9fb7c8d75.png)
Here we can see evidence of the C --> T in the tumor bam but not in the normal sample. This variant is called correctly and is called by both mutect2 and varscan2.

```
9       130634110       130634110       C       T       exonic  AK1     nonsynonymous_SNV       0/0:48,0:0.02:48:22,0:26,0:21,27,0,0:.:.:.:.                0/1:50,30:0.376:80:25,10:25,19:26,24,14,16:.:.:.:.       0/0:0:.:47:.:.:.:.:47:0%:21,26,0,0      0/1:29:.:79:.:.:.:.:50:36.71%:26,24,15,14
```

Now input chr9:130634993. 

![image](https://user-images.githubusercontent.com/15352153/120906346-f14e8580-c615-11eb-9e22-e7f8cce02ef0.png)
Here we can see evidence of two variants in the AKT1 exonic region. These variants are called only using mutect2.

```
9       130634993       130634993       C       T       exonic  AK1     nonsynonymous SNV       0/0:46,0:0.021:46:24,0:21,0:18,28,0,0   0/1:48,8:0.155:56:24,5:24,3:20,28,5,3   ./.:.:.:.:.:.:. ./.:.:.:.:.:.:.
9       130635011       130635011       C       G       exonic  AK1     nonsynonymous SNV       0/0:39,0:0.024:39:23,0:16,0:17,22,0,0   0/1:46,9:0.175:55:21,5:25,4:16,30,5,4   ./.:.:.:.:.:.:. ./.:.:.:.:.:.:.
```

Now input chr9:130316821. This will be an intronic variant in the NIBAN2 gene. 

![image](https://user-images.githubusercontent.com/15352153/120906325-b8aeac00-c615-11eb-8163-6a74c044f85a.png)
This variant is real and only called by varscan2.

```
9       130316821       130316821       C       T       intronic        NIBAN2  .       ./.:.:.:.:.:.:. ./.:.:.:.:.:.:.                                     0/0:.:23:20:0:0%:13,7,0,0        0/1:.:24:19:5:20.83%:14,5,1,4
```

Short break time 
   
**Copy number variations**

In this workshop, we will present the main steps for calling copy number variations (CNVs). 
Normally we would perform a copy number variation analysis on alignment files (BAMs). Due to time and resource constraints we will not be able to do this analysis on the full sequence or from shortened bam files. Instead we will be working with some preprocessed data using gc correction files and mpileups. Our data sample is again the cagekid sample c0098.

To examine the copy number states within the cage kid tumors we are going to be using a tool called [controlfreec] (http://boevalab.inf.ethz.ch/FREEC/)

With this tool we can call 
* Copy number gains
* Copy number losses
* Copy-neutral loss of heterozygosity
* This tool can also be used on exome seqencing and can even be used without a control (non-tumor) sample.

#controlfreec requires additional files
* A reference genome and a reference index (fasta and fai):  
* List of single nucelotide polymorphisms 


Running controlfreec 

```
cd ~/workspace/Module4_somaticvariants
mkdir copynumber
cd copynumber
```

The first step is to generate a configuration file which contain the parameters and file locations.
```
bash /home/ubuntu/CourseData/CAN_data/Module4/scripts/generate_controlfreec_configurationfile.sh > CBW_config.txt
```
Let's now look at the configuration file. 

```
less CBW_config.txt
```
Here we see four categories general, sample, control and BAF.

We have commented out the BAM files, however if you you were running this for the first time you should start from the BAM files. 
Instead we use pileups and gc content. 
 * Pileup contain base information that will help us assess which alleles.
 * cpn contain gc content for normalization.
 
- [General] indicates parameters for how we want the algorithm to run; parameters like minCNAlength, intercept are indicating WGS is being used
- [sample] indicates the input data of our tumor either in bam format or pileup/cpn
- [control] indicates our input data of our normal sample, again in different formats
- [BAF] indicates the necesarry files for calculatibng the b-allele frequency which will help us indicated LOH

#Now to run controlfreec
```
/usr/local/bin/freec -conf CBW_config.txt > output.log
``` 
Note here we will get somewarning that segmentation is completed, this due to subsetting the files to allow for a quicker runtime. 

Let's list the files we generated here. 
```
ls *
```
     CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_BAF.txt : Containing the allele frequency 
     CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_ratio.txt : Contain read-depth ratios
     CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_CNVs : Containing our copy numbers calls th
     CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_info.txt : Containing sample information, purity and ploidy

While we don't have annovar we can simply annotate our segment file with a gene transfer format file or [GTF] (http://m.ensembl.org/info/website/upload/gff.html). 

```
bedtools intersect -wb -b <(less CBW_regions_c0098_Tumor.sorted.markduplicates.bam_CNVs | awk 'NF==7' | awk '{print "chr"$1"\t"$2"\t"$3"\t"$0}' | less -S) -a <(less ~/references/genes_and_GTFS/Homo_sapiens.GRCh38.Ensemble100.FullGeneAnnotations.txt | awk '$4=="ensembl_havana"') | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$15}' >  AnnotatedCBW_regions_c0098_Tumor.sorted.markduplicates.bam_CNVs.tsv
```
```
less -S AnnotatedCBW_regions_c0098_Tumor.sorted.markduplicates.bam_CNVs.tsv 
```
We can examine a gene of interest, such as a loss of VHL which is common in kidney cancers. 

```
less AnnotatedCBW_regions_c0098_Tumor.sorted.markduplicates.bam_CNVs.tsv | awk '$4=="VHL"' | less -S
```

If you don't have a specific gene or gene set you can use [COSMIC](https://cancer.sanger.ac.uk/census): Which contains genes that are known and speculated to be associated with cancer. 

Now lets go visualize these results using R --> You will need to open r studio for this.
This script is availabe with controlfreec but due to the subsettting we will have to plot the ratio and BAF ourselves.

* If you are missing any of the data it can be found at github  
* https://github.com/bioinformatics-ca/CAN_2021/raw/main/Module4/Data/CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_BAF.txt
* https://github.com/bioinformatics-ca/CAN_2021/raw/main/Module4/Data/CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_ratio.txt

 Let's read in our data
 ```
ratio_dataTable <- read.table(file = "CBW_regions_c0098_Tumor.sorted.markduplicates.bam_ratio.txt", header=TRUE)
#ratio_dataTable <- read.table(file = "https://github.com/bioinformatics-ca/CAN_2021/raw/main/Module4/Data/CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_ratio.txt",header=TRUE)

ratio <- data.frame(ratio_dataTable)
BAF_dataTable <-read.table(file="CBW_regions_c0098_Tumor.sorted.markduplicates.bam_BAF.txt", header=TRUE)
#BAF_dataTable <- read.table(file = "https://github.com/bioinformatics-ca/CAN_2021/raw/main/Module4/Data/CBW_regions_c0098_Tumor.sorted.markduplicates.bam_sample.cpn_BAF.txt"", header=TRUE)
BAF<-data.frame(BAF_dataTable)
ploidy <- 2                   
```
                   
Now we can plot all of our data in a way that allows us to see gain, losses and neutral regions         
```{r}         
for (chrom in c(3,5,11))
{ 
tt <- which(ratio$Chromosome==chrom)
if (length(tt)>0)
{
        plot(ratio$Start[tt],log2(ratio$Ratio[tt]),xlab = paste ("position, chr",i),ylab = "normalized copy number profile (log2)",pch = ".",col = colors()[88],cex=4)
        tt <- which(ratio$Chromosome==chrom  & ratio$CopyNumber>ploidy )
        points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = ".",col = colors()[136],cex=4)

        tt <- which(ratio$Chromosome==chrom  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
        points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = ".",col = colors()[461],cex=4)
        tt <- which(ratio$Chromosome==chrom)
        }
        tt <- which(ratio$Chromosome==chrom)
}
```      
   
```{r}
 

for (i in c(3,5,11))
{
tt <- which(BAF$Chromosome==i)
if (length(tt)>0)
  {
                lBAF <-BAF[tt,]
                plot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste ("position, chr",i),ylab = "BAF",pch = ".",col = colors()[1])

                tt <- which(lBAF$A==0.5)
                points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[92])
                tt <- which(lBAF$A!=0.5 & lBAF$A>=0)
                points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[62])
                tt <- 1
                pres <- 1

                if (length(lBAF$A)>4) {
                        for (j in c(2:(length(lBAF$A)-pres-1))) {
                                if (lBAF$A[j]==lBAF$A[j+pres]) {
                                        tt[length(tt)+1] <- j
                                }
                        }
                        points(lBAF$Position[tt],lBAF$A[tt],pch = ".",col = colors()[24],cex=3)
                        points(lBAF$Position[tt],lBAF$B[tt],pch = ".",col = colors()[24],cex=3)
                }

                tt <- 1
                pres <- 1
                if (length(lBAF$FittedA)>4) {
                        for (j in c(2:(length(lBAF$FittedA)-pres-1))) {
                                if (lBAF$FittedA[j]==lBAF$FittedA[j+pres]) {
                                        tt[length(tt)+1] <- j
                                }
                        }
  points(lBAF$Position[tt],lBAF$FittedA[tt],pch = ".",col = colors()[463],cex=3)
  points(lBAF$Position[tt],lBAF$FittedB[tt],pch = ".",col = colors()[463],cex=3)
}}}

 
``` 
                                                            
![image](https://user-images.githubusercontent.com/15352153/120875901-6dd85a00-c56b-11eb-9bf9-523629c0550a.png)
![image](https://user-images.githubusercontent.com/15352153/120876052-e0493a00-c56b-11eb-9381-54e82ce7d7b7.png)
                                                            
Based off of the ratio data we see that chromosome 3 has a diploid region that is then going into a gain:          
However we see here that this is a copy-neutral loss of heterozygosity. This is a form of copynumber variation where the cell losses one allelic region but gains the other. So instead of AB it will be AA or BB. 
                                                            
 Now we see a chromsome 5 region which represents a gain and it is supported by the BAF
![image](https://user-images.githubusercontent.com/15352153/120875932-89dbfb80-c56b-11eb-8509-ad6e5a435717.png)
![image](https://user-images.githubusercontent.com/15352153/120876231-d411ac80-c56c-11eb-8dff-8732ed98dedd.png)

 And we also see another gain on chromsome 11, which is again supported by the BAF.                                                           
![image](https://user-images.githubusercontent.com/15352153/120876295-323e8f80-c56d-11eb-88f7-5af71c879b49.png)
![image](https://user-images.githubusercontent.com/15352153/120876298-38cd0700-c56d-11eb-9fe9-05bfc48f7b74.png)
                                                   
                                
    
    

 

 
