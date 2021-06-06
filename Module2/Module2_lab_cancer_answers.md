---
layout: tutorial_page
permalink: /CAN_2021_Module2_lab_cancer
title: CAN 2021 Module 2 Lab answers
header1: Workshop Pages for Students
header2: Cancer Analysis 2021
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/CAN_2021
description: CAN 2021 Module 2 Lab answers
author: Heather Gibling, Sorana Morrissy, and Florence Cavalli
modified: June 6, 2021
---

# CAN Module 2 - Viewing cancer alignments in IGV 

## Suggested Answers
Below are suggested answers to the questions from lab 2. You might have thought of additional reasons and answers!


# Visualization Part 2: Inspecting small variants in the normal sample

## Heterozygous and Homozygous SNPs
1. Which variant is heterozygous and which is homozygous?
    * Left is heterozygous (~half the reads are marked as mismatches at this location) and right is homoszygous (all are marked as mismatches)
3. What are the variant allele frequencies for each SNP? To find out you can click/hover on the colored bars in the coverage track
    * 59% (left) and 100% (right)
5. Do these look like true SNPs? What evidence is there for this?
    * Yes, the allele frequencies are close to what we would expect to see for heterozygous and homozygous SNPS, and all of the mismatched bases are of high quality, and there is no strand bias. Additionally, they both line up with known SNPs in the dbSNP track.
7. Are these sequencing errors, SNPs, or SNVs?
    * Sequencing errors and SNVs. The low quality mismatches are most likely sequencing errors. The higher quality mismatches that only occur in one or two reads are most likely somatic SNVs. They are not SNPs because they are not known germline mutations.
9. What does “Shade base by quality” do and how might this be helpful?
    * Adjusts the mismatch color to be dark when the base is of high quality and light when low quality. This helps us differentiate between sequencing errors and SNVs.
11. Can a normal sample have somatic SNVs?
    * Yes! There is always a chance there will be a mutation or error during genome replication

## Homozygous deletion 
1. How large is this deletion?
    * 24bp as indicated in the black horizontal line indicating a break in the read (relative to the reference genome)
3. Is it homozygous or heterozygous?
    * Homozygous-no reads span this region and there is a clear sharp drop in the coverage track.

## GC coverage
1. What do you notice about the coverage track and GC Percentage track?
    * The coverage track slowly dips down to 0 and then back up, and the GC Percentage track has a clear drop to 0%
3. What does coloring alignments by “read strand” tell you?
    * No forward reads are to the left of this dip and no reverse reads are to the right.
5. Do you think this is a deletion? Compare this region to the previous region.
    * This is a loss in coverage due to no GC content, not a deletion. Regions of the genome with a low GC% are notoriously difficult to sequence. The lack of reads spanning this region indicates these fragments were not able to be sequenced at all.

# Visualization Part 3: Inspecting small somatic variants in the tumor sample

## Somatic SNVs
1. How many SNVs are in this region for each sample?
    * Normal: 3, Tumor: 4
3. What is the variant allele frequency for the extra SNV in the tumor sample? How did it get this high?
    * 37%: The cells with this mutation were had an advantage during growth of the tumor sample and the mutation became more frequent

## Somatic SNP with change in heterozygosity
1. What are the variant allele frequencies for each sample?
    * Normal: 59%, Tumor: 80%
3. Why are these frequencies different?
    * Tumor cells with the alternate allele had an advantage during growth of the tumor sample and the variant became more frequent

## Somatic indel next to SNP with change in heterozygosity
1. What type of variant is in your centre line?
    * Small deletion
3. What do you notice about the variant in the normal sample to the right?
    * The variant allele is less frequent in the tumor sample and the reads with the deletion do not have the variant allele
5. What might be an explanation for what happened?
    * This deletion is in linkage disequilibrium with the reference allele (they are very close together). The deletion originated in tumor cells without the variant SNP allele, and as these cells had an advantage during growth of the tumor sample, the SNP variant became less frequent

# Visualization Part 4: Inspecting structural variants in NA12878

## Inversion
1. What do you notice about the blue and teal reads at the top?
    * The blue read pairs are both reverse facing, whereas the teal read pairs are both forward facing
3. How does “View as pairs” help understand that this is an inversion?
    * It shows us that not all read pairs in this region are forward-reverse as expected
5. Is this inversion heterozygous or homozygous?
    * Heterozygous, as there are lots of properly oriented read pairs as well

## Duplication
1. What do you notice about the green reads at the top?
    * The read pairs are in the opposite orientation as expected (reverse-forward instead of forward-reverse)

## Large deletion
1. What do the red read pairs indicate?
    * Read pairs with a larger insert size than expected. Since the senquence in between is missing in the sample, the read pairs came from an appropriately sized fragment, but when mapped to the reference genome they are much further apart than expected if the reference sequence was present
3. What other track can we look at to see that this is a deletion?
    * There is a dip in the coverage track indicating fewer reads mapped here


## Resources
* [IGV user manual](http://software.broadinstitute.org/software/igv/UserGuide)
