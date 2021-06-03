---
layout: tutorial_page
permalink: /BiCG_2019_Module11_lab_answers
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: BiCG 2019 Module 11 Lab Answers
author: Robin Haw
modified: March 5, 2019
---

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Module 14 Practical Lab: Reactome

By Robin Haw

CBW Lab Module 11 Answers

Example 1.
1.	The overall sub-network consists of 287 nodes and 713 edges. The largest component of the subnetwork consists of 250 nodes and 614 edges, with the remainder of nodes and edges distributed amongst 8 other small subnetworks and interactions.
2.	Couple of ways to answer this. The driver mutations are probably the frequently mutated gene in the samples. The node size is proportional to the number of samples where the gene is mutated. Method 1- Look for the largest nodes in the diagram. Method 2 – Click Node Table and sort by “sampleNumber”. The largest node is TP53, ie. mutations in the TP53 gene are highly prevalent, occurring in 100 samples. Other driver mutations include EGFR (95) and PTEN (93). Additional mutations of interest include NF1, PIK3R1, PIK3CA, PIK3R1, RYR2, RB1.
3.	Search for “TP53 PEG3” in search bar in top right of Cytoscape tool. Annotated Functional Interaction based upon data from the TRED database. This targeted interaction describes an interaction between TP53 (regulator) and PEG3 (target). An immunoprecipitation experiment demonstrates the interaction, and the supporting evidence has been published in the paper with a PubMed ID: 11679586.
4.	Search for “TAF1 TAF7L” in search bar in top right of Cytoscape tool. Predicted Functional Interaction based upon data (2/9 sources are true) from a mouse interaction database and GO (GO BP sharing). FI score: 0.53
5.	21 modules, with 8 modules of 10 ≥ genes.
6.	18 modules, depending on the results of the enrichment analysis. Some pathways gene sets at the cutoff threshold may come or go but those highly significant gene sets are always there.
7.	0: TP53 Signaling 1: RTK signalling, 2: ECM and Integrin signalling.


Example 2.
1.	The overall sub-network consists of 272 nodes and 569 edges. The largest component of the subnetwork consists of 240 nodes and 510 edges, with the remainder of nodes and edges distributed amongst 14 other small subnetworks and interactions.
2.	The largest node is TP53, ie. mutations in the TP53 gene are highly prevalent, occurring in at least 96% of HGS-OvCa samples.
3.	After clustering, there are 23 modules with 11 modules of 10 ≥ genes.
4.	8 modules, depending on the results of the enrichment analysis. Some pathways gene sets at the cutoff threshold may come or go but those highly significant gene sets are always there.
5.	0: TP53 Signaling 1: ECM and Integrin signalling, 2: EGFR Signaling and Inositol phosphate metabolism, 7: Calcium signalling-Adreneric Signaling-Cardiac Muscle Contraction,
6.	Yes, ECM and Integrin signalling.
7.	Nuclear components - Nucleoplasm, nuclear membrane, nuclear pore, chromatin, etc.
8.	Modules 0, 2 and 7 will be highlighted. Navigate through hierarchy. Neoplasm > Neoplasm_by_Site > Breast Neoplasm > Maligant_Breast_Neoplasm > Breast Carcinoma > Stage_IV_Breast_Cancer.  Go back to the Network Module Browser. Genes in the modules that have ‘Stage IV Breast Cancer’ annotations will be yellow-highlighted: BRCA1, NRG1, TP53, INSR, EGFR.
9.	EGFR.
10.	3 module: 1, 3 and 7.

 
![img1](https://github.com/bioinformaticsdotca/CSHL_2019/blob/master/Module14/Reactome1-1.png?raw=true)  

11.	The ReactomeFIViz app splits samples into two groups: samples having genes mutated in a module (green line), and samples having no genes mutated in the module (red line). The plugin uses the log-rank test to compare the two survival curves, and estimates p-values. In Modules 4 (KM: p=0.00489), patient with genes mutated (green line) have a better prognosis than patients with no gene mutations (red line). Module 7 is most statistically significant modules from the CoxPH and KM analysis.
  
![img2](https://github.com/bioinformaticsdotca/CSHL_2019/blob/master/Module14/Reactome2-1.png?raw=true)   

11.	In Module 7, the Calcium signaling, Chemical Synapse/Neurotransmission and Muscle Contraction annotations reflect a shared set of genes. These genes represent voltage-gated ion channels, which are a group of transmembrane ion channels that activated by changes in electrical potential difference. Even though ion channels are especially critical in neurons and muscle tissue, they are common in many types of cells, controlling the influx and outflux of ions. There are a number of genetic disorders, which disrupt normal functioning of ion channels. Calcium homeostasis is essential for cell migration, and tumor metastasis in particular. It may be that mutations in Module 7 genes disrupt calcium homeostasis, thereby impairing the tumour’s ability to metastasize, and extending patient’s overall survival.


