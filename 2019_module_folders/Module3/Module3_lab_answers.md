# BiCG Module 3 - IGV lab

## Suggested Answers

Below are suggested answers to the questions from lab 3. You might have thought of additional reasons and answers!

---

# Visualization Part 2: Inspecting SNPs, SNVs, and SVs

## Neighbouring Somatic SNV and Germline SNP 

**Questions:**
* What does "Shade base by quality" do? How might this be helpful?
    * Distinguishes high-quality bases from lower quality bases. If the T alternate allele was predominantly light/transparent, we wouldn't be very confident that it is a probable SNV. The only light-colored T in this example is at the very end of a read, where base quality is generally lower.
* How does "Color by read strand" help?
    * Lets us see that the T alternate allele is present on both forward and reverse sequencing reads. If it was present on only forward reads or only reverse reads, that could indicate a sequencing artifact as opposed to a probable SNV.

## Homopolymer Repeat with Indel 

**Question:**
* Is the "T" likely a valid SNV? What evidence suggests it is or isn't?
    * Probably not. 6/36 reads show a T at this positition (numbers obtained by clicking/hovering on the coverage track block at this position), and 4 of those 6 are light colored (low quality). Additionally, only forward (red) reads have the T alternate allele, indicating a strand bias sequencing error.

## Coverage by GC 

**Question:**
* Does the coverage correspond to the GC content?
    * Yes! Read coverage is low where GC content of the reference genome is very low (and also very high, but in these data the low coverage effect is easier to see). Seqeuencing technologies are not perfect and have biases, meaning not all regions of the genome are covered equally. You can read more about the effect of GC content on short-read sequencing coverage [here](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-5-r51).

## Low Mapping Quality 

**Question:**
* Why do LINE elements affect mapping quality?
    * LINEs (long interspersed nuclear elements) are retrotransposons found throughout the genome. Because these sequences are repeated in many different locations, when a read aligns to part of this sequence, it can often map equally well to multiple locations in the genome. This is not helpful, so aligners generally give low mapping quality scores reads that do not map uniquely to one location (i.e. 0; the reads will also be white instead of grey in IGV).

## Homozygous Deletion 

**Question:**
* What other track provides evidence of a deletion at this location? 
    * The coverage track drops down to zero at the same location as the gap in sequencing reads.

---

# Visualization Part 4: Visualizing Long Reads

## Cleaning up Sequencing Error Noise

**Question:**
* What does the abundance of dashed lines and purple "I"s tell us about the types of errors produced by nanopore sequencing?
    * Nanopore sequencing is prone to insertions and deletions. "I"s are insertions and the dashed lines are deletions, relative to the reference genome.

## Viewing Variants

**Question:**
* Some of the coloring options we used for viewing the previous sample (HCC1143) are not available for this NA12878 bam, such as *View as pairs* and *Color alignments by insert size and pair orientation*. Why is this?
    * Nanopore sequencing doesn't do paired read sequencing! Instead, long fragments of DNA are pulled through pores and are sequenced until the end of the fragment (or until the pore wears out). Reads are considered to be independent from each other in this regard. This is different from short-read sequencing platforms like Illumina that sequence both ends of a DNA fragment but are usually unable to sequence the middle segment.