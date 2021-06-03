---
layout: tutorial_page
permalink: /BiCG_2019_Module11_lab
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
description: BiCG Module 11 Lab
author: Robin Haw
modified: June 26th, 2018
---

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Module 11 Practical Lab: Reactome

By Robin Haw

**Aim:** This exercise will provide you with an opportunity to perform pathway and network analysis using the Reactome Functional Interaction (FI) and the ReactomeFIViz app. 

**Goal:** Analyze gene lists and somatic mutation data to identify biology that contributes to GBM and ovarian cancer.

**Example 1: Network-based analysis of GBM gene-sample data** 
1. Open up Cytoscape.   
a.	Go to Apps>Reactome FI and Select “Gene Set/Mutational Analysis”.    
b.	Choose “2018 (Latest)” Version.   
c.	Upload/Browse [GBM_genesample.txt](https://raw.githubusercontent.com/bioinformaticsdotca/HT-Biology_2017/master/GBM_genesample.txt) file.   
d.	Select “Gene/sample number pair” and Choose sample cutoff value of 4.   
e.	Select “Fetch FI annotations”.   
f.	Click OK.  

2.	Describe the size and composition of the GBM sub-network?  
3.	What are the most frequently mutated genes?
4.	Describe the TP53-PEG3 interaction, and the source information to support this interaction?  
5.	Describe the data sources for the TAF1-TAF7L FI?  
6.	After clustering, how many modules are there?   
7.	How many pathway gene sets are there in Module 3 when the FDR Filter is set to 0.005 and Module Size Filter to 10?   
a.	Hint: Analyze Module Functions>Pathway Enrichment. Select appropriate filters at each step.  
8.	What are the most significant pathway gene sets in Module 0, 1, 2?  
a.	Hint: You don’t need to list them all!   

**Example 2: Network-based analysis of OvCa somatic mutation**   
1.	Open up Cytoscape.   
a.	Go to Apps>Reactome FI and Select “Gene Set/Mutational Analysis”.    
b.	Choose “2018 (Latest)” Version.   
c.	Upload/Browse [OVCA_TCGA_MAF.txt](https://raw.githubusercontent.com/bioinformatics-ca/bioinformatics-ca.github.io/master/2016_workshops/cancer/OVCA_TCGA_MAF.txt) file.   
d.	Select “NCI MAF” (Mutation Annotation File) and Choose sample cutoff value of 4.   
e.	Do not select “Fetch FI annotations”.   
f.	Click OK.  

2.	Describe the size and composition of the OvCa network?  
3.	What are the most frequently mutated genes?  
4.	After clustering, how many modules are there?   
5.	How many pathway gene sets are there in Module 3 when the FDR Filter is set to 0.005 and Module Size Filter to 10?  
a.	Hint: Analyze Module Functions>Pathway Enrichment. Select appropriate filters at each step.  
6.	What are the most significant pathway gene sets in modules, 0, 1, 2 and 7?   
7.	Do the GO Biological Process annotations correlate with the significant pathway annotations for Module 1?   
a.	Hint: Analyze Module Functions>GO Biological Process. Select appropriate filters at each step.  
8.	What are the most significant GO Cell Component gene sets in Module 0 when the FDR Filter is set to 0.005 and Module Size Filter to 10? [Optional]  
a.	Hint: Analyze Module Functions>GO Cell Component. Select appropriate filters at each step.  
9.	Are any of the modules annotated with the NCI Disease term: “Stage_IV_Breast_Cancer” [malignant cancer]?  
a.	Hint: Load Cancer Gene Index>Neoplasm>Neoplasm_by_Site>Breast Neoplasm>…….  
10.	What are the targets of Docetaxel?
a.	Hint: Overlay Cancer Drugs>Fetch Cancer Drugs. Maybe apply filters? 
11.	How many modules are statistically significant in the CoxPH analysis?   
a.	Hint: Analyze Module Functions>Survival Analysis>Upload/Browse [OVCA_TCGA_Clinical.txt](https://raw.githubusercontent.com/bioinformatics-ca/bioinformatics-ca.github.io/master/2016_workshops/cancer/OVCA_TCGA_Clinical.txt). Click OK.  
12.	What does the Kaplan-Meyer plot show for the most clinically significant modules?  
a.	Hint: Click the most statistically significant module link [blue line] from the CoxPH results panel. Click OK. Click #_plot.pdf to display Kaplan-Meyer plot. Repeat this for the other significant module links. KM plot: samples having genes mutated in a module (green line), and samples having no genes mutated in the module (red line).    
13. Taking into consideration what you have about module 7, what is your hypothesis?

[Answer key](https://bioinformaticsdotca.github.io/BiCG_2019_Module11_lab_answers)
