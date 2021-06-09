---
layout: tutorial_page
permalink: /CAN_2021_module9_lab
title: CAN 2021 Module 9 Lab
header1: Workshop Pages for Students
header2: Cancer Analyis 2021 Module 9
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/CAN_2021
description: CAN 2021 Moule 9 lab
author: Robin Haw
 
---

# CBW Lab Module 9 
***By Robin Haw***  

================================

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

================================


## Data Source
please download this [zip file](https://drive.google.com/file/d/14gT-Qgwcvp8TL1hPxUePg1eQIIyFQFJr/view?usp=sharing)


## Aim
This exercise will provide you with an opportunity to perform pathway and network analysis using the Reactome Functional Interaction (FI) network and the ReactomeFIViz Cytoscape app. 

**Goal**: Analyze gene lists and somatic mutation data to understand the biology that contributes to multiple cancers.

### Example 1: Pathway-based analysis of GBM genelist (PMCID: PMC2671642)
* Open up Cytoscape. 
* Go to Apps >Reactome FI>Reactome Pathways.
* Unfurl the “Signal Transduction” events in the event hierarchy, by clicking the triangle to the left of the event name, in the “Reactome” tab on the left. 
* Click on your favourite pathway. 
* Right-click on the highlighted pathway name in the display drop-down menu, select “Analyze Pathway Enrichment” 
* Upload/Browse/"Copy-and-paste" the “GBM_genelist.txt” data into the configuration panel, and click “OK”. 

**Question**
What are the most significant biological pathways when the FDR Filter is set to 0.05?
<details>
  <summary>
**Hint ** (click here)
  </summary>
  
* Right-click on a pathway in the Table Panel, and click “View in Diagram”. Purple-coloured nodes reflect hits in the dataset. Right-click on highlighted nodes to invoke additional features
* "Open Reactome Reacfoam" from the pathway tree. You may also download the Reacfoam view by clicking the download button at the top-right corner. For windows 10 users, to open the Reacfoam view, you need to allow "public" access to Cytoscape by checking "public" in the settings for "Allow an app through Windows Firewall" in the "System and Security" control settings

</details> 
</details> 

### Example 2: Gene Set Enrichment Analysis of AML dataset (PMCID: PMC6849941)

* Open up Cytoscape. 
* Go to Apps >Reactome FI>Reactome Pathways.
* Unfurl the “Signal Transduction” events in the event hierarchy, by clicking the triangle to the left of the event name, in the “Reactome” tab on the left. 
* Click on your favourite pathway. 
* Right-click on the highlighted pathway name in the display drop-down menu, select “Perform GSEA Analysis” 
* Upload/Browse/"Copy-and-paste" the “AML_trametinib_genelist_rnk.txt” data into the configuration dialog, and click “OK”. 
* (Optional) Enter the gene score file and choose the minimum and maximum size of pathways along with the permutation number.

**Question**
What does the GSEA results tell you about the inflammatory profile associated with trametinib sensitivity?

### Example 3: Network-based analysis of OvCa somatic mutation (PMCID: PMC3163504)
* Open up Cytoscape. 
* Go to Apps>Reactome FI and Select “Gene Set/Mutational Analysis”.  
* Choose “2020 (Latest)” Version. 
* Upload/Browse “OVCA_maf.txt” file. 
* Select “NCI MAF” (Mutation Annotation File) and Choose sample cutoff value of 4. 
* Do not select “Fetch FI annotations”. 
* Click OK.

**Questions**
* 1. Describe the size and composition of the OvCa network?
* 2. What are the most frequently mutated genes?
* 3. Describe the TP53-PEG3 interaction, and the source information to support this interaction?
* 4. Describe the data sources for the EGFR- PTPRG FI?
* 5. After clustering, how many modules are there? 
* 6. How many pathway gene sets are there in Module 1 when the FDR Filter is set to 0.005 and Module Size Filter to 10?
<details>
  <summary>
**Hint** (click here)
  </summary>
 Analyze Module Functions>Pathway Enrichment. Select appropriate filters at each step.
</details> 

* 7. What are the most significant pathway gene sets in Module 0, 1, 2, and 3? 
* 8. Do the GO Biological Process annotations correlate with the significant pathway annotations for Module 0? 
<details>
  <summary>
**Hint ** (click here)
  </summary>
Analyze Module Functions>GO Biological Process. Select appropriate filters at each step.
</details> 

 * 9. What are the most significant GO Cell Component gene sets in Module 2 when the FDR Filter is set to 0.005 and Module Size Filter to 10? [Optional]
 <details>
  <summary>
**Hint ** (click here)
  </summary>
Analyze Module Functions>GO Cell Component. Select appropriate filters at each step.
</details> 

* 10. Are any of the modules annotated with the NCI Disease term: “Stage_IV_Breast_Cancer” [malignant cancer]?
 <details>
  <summary>
**Hint ** (click here)
  </summary>
Load Cancer Gene Index>Neoplasm>Neoplasm_by_Site>Breast Neoplasm>
</details> 

* 11. What are the targets of Axitinib?
 <details>
  <summary>
**Hint ** (click here)
  </summary>
Overlay Cancer Drugs>Fetch Cancer Drugs. Maybe apply filters? 
</details> 

* 12. How many modules are statistically significant in the CoxPH analysis? 
 <details>
  <summary>
**Hint ** (click here)
  </summary>
Analyze Module Functions>Survival Analysis>Upload/Browse “OVCA_clinical.txt”. Click OK.
</details> 

* 13. What does the Kaplan-Meyer plot show for the most clinically significant module?
 <details>
  <summary>
**Hint ** (click here)
  </summary>
	
* Click the most statistically significant module link [blue line] from the CoxPH results panel. Click OK. Click #_plot.pdf to display Kaplan-Meyer plot. Repeat this for the other significant module links. KM plot: samples having genes mutated in a module (green line), and samples having no genes mutated in the module (red line).
* There is a bug with the Windows version to view the .pdf file. You may want to search for the CytoscapeConfiguration folder.
* On Mac, go to: "/Users/[name]/CytoscapeConfiguration/3/karaf_data/tmp/#########_plot.pdf"

</details> 	 
* 14. Taking into what you have learned about module 4, what is your hypothesis?

### Example 4: Illuminating the interactions and functions of an understudied protein
* Point the web browser to https://idg.reactome.org. 
* On the home page that appears, in the search box in the middle of the page, click the text box, type [PRKY or another gene symbol], then press the Enter key (or click the Search button). 
* This brings up the search results in the bottom of the page.
* Refer to the User Docs: https://idg.reactome.org/documentation/userguide to guide you through the site.



