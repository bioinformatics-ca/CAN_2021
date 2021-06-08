Initializing module4 SNV, Indels and CNVs

Question: What is the definition of AD and AF in this mutect2-VCF?
Answer: AD is the allelic depth for the ref and alt alleles (respectively) and AF is the allele-frequency

Question: what does the "germline" in the FILTER field mean?
Answer: Evidence indicates this site is germline, not somatic

Question: what does -m-both accomplish here
Answer: This will split multiallelic variants. The "-" splits multiallelic variants and the both will work on both snps and indels.

Question: What is the definition of AD and RD in this varscan-VCF & how would you calculate the variant allele frequency?
Answer: AD is the depth of variant alleles and RD is the depth of References alleles.

Question: What happened in the LOH variant below? (chr9:130311404)
Answer: We increased the VAF meaning we lost counts of the reference.  

Question: What happened in the LOH variant below? (chr9:130346618)
Answer: We decreased the VAF meaning we 

Question: In what cancer type would having germline variants be useful?
Answer: Pediatric cancers have a low mutation burden and would not have many somatic variants. The cancer was therefore initated by something other than a high mutation burden. Germline variants are also useful in predisposition studies, for example if you have germline missense mutations in TP53, your odds of cancer are higher than those who don't have the mutation

Question: Compare the variant allele frequency between these three variants in AK1. Which one's came first?
Answer: The first variant likely came first, this is indicated by it's higher VAF. 37% > ~16%

Question: This variant is only called by varscan2, what red flags exist for this variant?
Answer: Here we can look at the reference sequence an we can see a repeat region CTTTTT repeat. We also see a strand bias, there are far more reads in favor of one strand than the other. 

Question: What does each dot mean in the top CNV plots and what does each dot mean in the bottom BAF plot?
Answer: In the CNV plot each dot represents a window which we set at 50000bp & In the BAF plot each dot represents a SNP location and the relative b-allele-frequency 

Question: What is the absolute copy-number of the red region in the top plot?
Answer: log(2)x=0.6 --> solving for x=1.59 * 2 --> absolute copy number of 3.03

Question: Comparing the green region to the BAF plots, what kind of copy number variation is taking place?
Answer: Here we look to have a diploid region as is indicated by the green color, however we have a see evidence for a loss, as the pink lines (BAF) are near 1.0 and 0 meaning we have a loss. This is a scenario known as copy neutral LOH. One way in which this can occur: AB ---lose B---> A ---gain A---> AA.

Question: What CNV events are supported by the TFRI GBM data?
Answer: We have loss on 6, whole gain of chromosome 7, isochromosome 9, CopyNeutral LOH on 10, isochromosome 18, lage gain on 19

Question: There are hundred of segments as seen in both CNV and BAF plots, what is one way you can modify controlfreec to call fewer segments?
Answer: One way would be to increase the window and step size. 
