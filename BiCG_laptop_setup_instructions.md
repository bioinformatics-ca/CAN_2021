---
layout: tutorial_page
permalink: /BiCG_laptop_setup_instructions
title: BiCG
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
home: https://bioinformaticsdotca.github.io/BiCG_2019
---

1) Install latest version of R which can be downloaded from http://probability.ca/cran/.

1b) Download and install the most recent version of [R Studio desktop](http://www.rstudio.com/).  If prompted to install git, select yes.

2) Install the BioConductor core packages. If you have installed R version 3.5.0 or higher, open R and at the '>' prompt, paste the commands:
 
```
install.packages("BiocManager");
library(BiocManager);
BiocManager::install();
```

If you already have an older version of R installed (3.4.4 or lower), open R and at the '>' prompt, paste the commands:

```
source("http://bioconductor.org/biocLite.R");
biocLite();
```

If you are unsure which version you have installed, open R and at the '>' prompt, enter the command:

```
version
```

3) A robust text editor.   

* For Windows/PC - [notepad++](http://notepad-plus-plus.org/)  
* For Linux - [gEdit](http://projects.gnome.org/gedit/)  
* For Mac – [TextWrangler](http://www.barebones.com/products/textwrangler/download.html)

4) A file decompression tool.  

* For Windows/PC – [7zip](http://www.7-zip.org/).  
* For Linux – [gzip](http://www.gzip.org).   
* For Mac – already there.

5) A robust internet browser such as Firefox or Safari (Internet Explorer and Chrome are not recommended because of Java issues).

6) Java -The visualization program that we will be using (IGV) requires Java. Check if you have Java installed: https://www.java.com/verify/ and download Java if you do not have it installed (You need Java 8. Do NOT install Java 10).

7) Integrative Genomics Viewer 2.4 (IGV) - Once java is installed, go to http://www.broadinstitute.org/igv/ and register in order to get access to the downloads page. Once you have gained access to the download page, click on the appropriate launch button that matches your computer's operating system   

8) SSH client - Mac and Linux users already have a command line ssh program that can be run from the terminal. For Windows users, please download [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html).  

9) SCP/SFTP client - We will be moving data from the servers to the student laptops for visualization. Mac and Linux users already have a command line scp and sftp program. For Windows users, please install [WinSCP](http://winscp.net/eng/download.php).

10) A PDF viewer (Adobe Acrobat or equivalent).

11) Install [Cytoscape 3.7.1](https://cytoscape.org/download-platforms.html).  

Choose the version corresponding to your operating system (OS, Windows or UNIX) 
Cytoscape requires Java8: check your version at  https://www.java.com/verify/ and download Java8 if you do not have it installed. Contact your system administrator if you have trouble with Java installation. 

12) Within the Cytoscape program, install the following Cytoscape apps.  

From the menu bar, select ‘Apps’, then ‘App Manager’.
 
Within 'all apps', search for the following and install:  

 * EnrichmentMap (3.2.0)
 * EnrichmentMap Pipeline Collection (it will install ClusterMaker2, WordCloud and AutoAnnotate) 
 * GeneMANIA (3.5.0)
 * Iregulon  
 * ReactomeFIPlugin - http://apps.cytoscape.org/apps/reactomefiplugin  
 * stringApp
 
 
13) Install the data set within GeneMANIA app.

From the menu bar, select 'Apps', hover over 'GeneMANIA', then select 'Choose Another Data Set'.  
From the list of available data sets, select the most recent (2017-07-13/13 July 2017) and under ‘Include only these networks:' select ‘all’. Click on ‘Download’.  
An ‘Install Data' window will pop-up. Select H.Sapiens Human (2589 MB). Click on ‘Install’.  
This requires time and a good network connection to download completely, so be patient (around 15mins).  

  
14) Install GSEA.  

Go to the [GSEA page](http://www.broadinstitute.org/gsea/index.jsp)    
Register  
Login  
In menu, choose Downloads  
Go to the javaGSEA Java Jar file section and download the gsea-3.0.jar file and save in your Documents folder (do not leave it in the “Downloads”folder).  
 
To run GSEA during the workshop, you must use the command line. You will need to open a terminal and execute the install commands. Since we will need to run GSEA this same way each time, it will be a good idea to save this information on how to run GSEA.
 
**MAC/Linux Computer** 

* On a MAC, the Terminal window is located in Applications/Utilities. Tip: add the terminal window to your dock so it is easy to open when needed.  
* At the prompt, type the command in your terminal window and hit enter:

```
java -Xmx4G -jar ~/Documents/gsea-3.0.jar
```

**PC/Windows Computer** 

* On Windows, go to the start icon and type cmd (for command prompt) in the search box.  
* At the prompt, type the following commands, hitting enter in between each command and waiting for the prompt before the next command:

```
cd Documents
java -Xmx4G -jar gsea-3.0.jar
```
