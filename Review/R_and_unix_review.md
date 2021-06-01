---
layout: tutorial_page
permalink: /BiCG_Unix_and_R
title: BiCG Unix and R Review
header1: Workshop Pages for Students
header2: Bioinformatics for Cancer Genomics 2019
image: /site_images/CBW_cancerDNA_icon-16.jpg
author: ?, edited by Heather Gibling
home: https://bioinformaticsdotca.github.io/BiCG_2019
---


# CBW BiCG Unix and R Review Session  

## Introduction

Welcome to the review for **Unix and R**! This lab will introduce you to the the command line and using R.

After this lab, you will be able to:

* use common commands in Unix
* read data into R
* produce basic plots in R

## Before We Begin

Read these [directions](http://bioinformaticsdotca.github.io/AWS_setup) for information on how to log in to your assigned Amazon node.

## Unix Review

The unix terminal/shell/bash console is a text interface for interacting with a computer. It lets you see the folders and files saved on a computer, read and write new files, and other activities.

A lot of biological data files are very large, and running them on your laptop is not feasible. Researchers use very large computer clusters or cloud computing in order to generate, manipulate, store, and view data files. These clusters are accessed remotely from another computer, such as your laptop. Because these clusters generally require use of the unix terminal in order to perform these actions, a lot of bioinformatics software is written specifically for use in a unix terminal. Therefore it is helpful to learn the basics of navigating and using a terminal.

You can't use a mouse and click to move the cursor--everything is done by typing. If you make a typo while typing a command, you need to use the arrow keys to navigate to the position in order to edit. A short cut to jump to the beginning of the line is `ctrl + a`, and `ctrl + e` will jump to the end of the line.

You can use the mouse to select and highlight text to copy/paste. To copy/paste within the terminal, you can do the following:
* Mac: `cmd + c` and `cmd + v`
* Windows: highlight text and click the right mouse button
* Linux: `ctrl + shift + c` and `ctrl + shift + v`

### Basic commands:

* ls: list the files in the directory:  
```
ls
```

* cd: change directories: 
```
cd workspace
```

* mkdir: make a new directory:  
```
mkdir Review_Session
```

* `cd ..` will take you up one directory. Give it a try.

* pwd: print working directory (where you are currently located): 
```
pwd
```
You should be in `/home/ubuntu`. 

* Let's move into the new directory we created:
```
cd workspace/Review_Session
```

* echo: print what was typed: 
```
echo 'Hello World!'
```

* `>` redirects the output from a command into a new file:  
```
echo 'Hello World!' > test.txt
```

* cat: print the contents of a file:  
```
cat test.txt
```

* curl: used to get contents from a URL:  
```
curl https://raw.githubusercontent.com/bioinformaticsdotca/Genomic_Med_2017/master/test.fasta > test.fasta
curl https://raw.githubusercontent.com/bioinformaticsdotca/BiCG_2017/master/Review_session/Gene_R_example.csv > Gene_R_Example.csv
```

* head and tail: get the beginning or end of a file:  
```
head test.fasta
tail test.fasta
```

* less: look at the contents of a file in a new window (handy for looking at large files when you don't want to loose sight of your previous commands in the terminal):  
```
less test.fasta
```
To exit `less`, press `q`.

* cp: copy a file
```
cp test.fasta test2.fasta
```

* mv: move a file
```
mkdir test-files
mv test2.fasta test-files
ls
ls test-files
```
You should only see test2.fa listed in the test-files directory. 

* Let's move it back:
```
mv test-files/test2.fasta .
```
The single `.` refers to the current directory.

* `mv` can also be used to rename files
```
mv test2.fasta test3.fasta
ls
```
You should see test.fasta and test3.fasta, but not test2.fasta

* rm: remove files  
```
rm test3.fasta
```
Note: Removing files with `rm` **permanently deletes** them, so be careful! You can using `rm -i` to avoid accidentally deleting a wanted file. This forces you to confirm by typing `y`.  

* grep: search for a matching pattern in a file:
```
grep TTT test.fasta
```

* You can use the pipe character `|` to string commands together:
```
head test.fasta | grep TTT
```

* Flags are additional options for commands. They are usually indicated with a dash and a single letter:
```
head -n 5 test.fasta
```
`-n` is the flag to specify the number of lines to be shown for the `head` command, and `5` is the option we chose.

* man: open a help page for a command:
```
man head
```
Scroll through the usage and options available for the `head` command. Press `q` to exit when done.

* history: prints out the \~1000 most recently commands that you entered:
```
history
```

**Useful tips:**

* You can get commands that you previously used to reappear in the terminal by pressing the up arrow on your keyboard.
* If you want to quickly delete everything that you typed in the terminal, or if you want to cancel a command that is taking a long time to run, hold `ctrl` and type `c`.


***

## R Review

### Connecting To RStudio

* Open an internet browser.  
* In the URL bar, enter *http://##.oicrcbw.ca:8080* replacing *##* with your provided student number.  
* Enter the supplied username and password.

### Navigating RStudio

RStudio is an interactive program for reading and writing R code. It is divided into four panes that can be resized and rearranged to your preference:
* Source
* Environment, History, and Connections
* Console
* Files, Plots, Help, and Viewer

The *console* pane is similar to the bash/unix/terminal console. You can enter code here and see output, but you can't save the code you've written. It's useful for checking quick things, like the size of a list of genes.

The *source* pane is where you write and edit scripts or RStudio notebooks (where you write code that you want to save). You can select code to run and the results will show up in the console (as well as the source, if you are using a notebook).

The *environment, history, and connections* pane gives an overview of what items you have entered in your R session (*environment*) and the history of commands you've typed (*history*). (We don't need to worry about the *connections* tab)

The *files, plots, packages, help, and viewer* pane lets you navigate folder and load files (*files*), view images you've created (*plots*), see what packages have been installed and loaded in your R session (*packages*), and look at help files for specific packages and functions (*help*). (We don't need to worry about the *viewer* tab)

### RStudio Notebooks

Information on RStudio Notebooks can be found [here](http://rmarkdown.rstudio.com/r_notebooks.html).  

RStudio notebooks are written in R Markdown and contain text that can be executed independently.  

To run code chunks, place your cursor within the code chunk and press *Cmd+Shift+Enter* on Mac or *Crtl+Shift+Enter* on Windows or Linux, or click the green triangle run button at the top right of the grey code chunk.  The output of the chunk appears below the code chunk.  

To start a new notebook, "File" > "New File" > "R Notebook". It will be populated with some information about notebooks--you can delete this info if you want.

Copy the code [here](https://raw.githubusercontent.com/bioinformaticsdotca/BiCG_2019/master/Review/R_review_notebook.Rmd) into the notebook you have created. 

The code below is the same as the link above. If you'd prefer to work in an R script instead of an RStudio notebook, you can copy/paste what's in the code boxes. Note that if you don't want text notes to be interpreted as code, they have to start with `#`

### Intro to R

R is used to manipulate various kinds of data. Data objects can be given names, called variables:

```r
5 + 6

number1 <- 5
number2 <- 6
number1 + number2
number.sum <- number1 + number2

# print the content of a data object or variable to screen
print(number.sum)

# typing just the name of the object also prints it to screen
number.sum
```

Data in R are manipulated primarily with functions. Functions are called by typing the name of the function followed by round brackets. Most functions require arguments, or options, to be specified in the brackets. For `print(number)` used above, `print` is the name of the function, and `number.sum` is the variable option that we want to print.

If you are not sure what a specific function does, or if you need a reminder on what the arguments for a funtion are, you can view the help page by typing `?function.name` without the round brackets:

```r
?print
```

Note: The help page will open in the "Help" tab in the Files/Plots/Packages/Help/Viewer pane.

### Getting Around

#### The Hard Way

```r
# get the current working directory:
getwd()

# set a new working directory:
setwd("~/workspace/Review_Session")

# to set a working directory on your own computer, you can do one of the following:
# (these commands are commented out so that these code chunks will run on AWS RStudio)

# setwd("C:/myPATH") # on Windows
# setwd("~/myPATH") # on Mac
# setwd("/Users/david/myPATH") # on Mac

# list the files in the current working directory:
list.files()

# list objects in the current R session:
ls()
```

#### The Easy Way

In RStudio we can use "Session" > "Set Working Directory" > "Choose Directory".  
  
### Data Types

#### Vectors

Vectors contain multiple pieces of data. The elements of a vector must all be of the same type (numeric, logical, character).
The `c()` function combines the items inside the round brackets and can be used to create a new vector:


```r
numeric.vector <- c(1,2,3,4,5,6,2,1)
numeric.vector

character.vector <- c("Fred", "Barney", "Wilma", "Betty")
character.vector

logical.vector <- c(TRUE, TRUE, FALSE, TRUE)
logical.vector
```

To refer to specific elements in the vector, use square brackets.
Because vectors are one dimensional (unlike a two-dimensional matrix), only a single number can be specified in the brackets.
If you want more than one element, you can specify a range (with `:`) or a vector (with `c()`):

```r
character.vector
character.vector[2]
character.vector[2:3]
character.vector[c(2,4)]
```

#### Matrices

You can create a 3x4 numeric matrix with:

```r
matrix.example <- matrix(1:12, nrow = 3, ncol=4, byrow = FALSE)
matrix.example

matrix.example <- matrix(1:12, nrow = 3, ncol=4, byrow = TRUE)
matrix.example
```

Alternatively, you can create a matrix by combining vectors:

```r
dataset.a <- c(1,22,3,4,5)
dataset.b <- c(10,11,13,14,15)
dataset.a
dataset.b

rbind.together <- rbind(dataset.a, dataset.b)
rbind.together

cbind.together <- cbind(dataset.a, dataset.b)
cbind.together
```

To get elements of the matrix, square brackets are used again, but two dimensions must be specified.
To do this, a comma is used, where the number before the comma is the row number and the number after the comma is the column number.
If you want to specify all of the rows, leave the space before the comma blank.
If you want to specify all of the columns, leave the space after the comma blank:

```r
matrix.example[2,4]
matrix.example[2,]
matrix.example[,4]
```

You can add column and row names to the matrix and use the new names to get the elements of the matrix:

```r
colnames(matrix.example) <- c("Sample1","Sample2","Sample3","Sample4")
rownames(matrix.example) <- paste("gene",1:3,sep="_")
matrix.example

matrix.example[,"Sample2"]
matrix.example[1,"Sample2"]
matrix.example["gene_1","Sample2"]
```

Note that all columns in a matrix must have the same type(numeric, character, etc.) and the same length.


#### Dataframes

Dataframes are similar to matrices, but different columns can have different types (numeric, character, factor, etc.): 

```r
people.summary <- data.frame(
                             age = c(30,29,25,25),
                             names = c("Fred", "Barney", "Wilma", "Betty"),
                             gender = c("m", "m", "f", "f")
                             )
people.summary
```

Getting elements of a dataframe is similar to getting elements of a matrix:

```r
people.summary[2,1]
people.summary[2,]
people.summary[,1]
```

An easier way to specify a dataframe column by name is by using `$`:

```r
people.summary$age
```

#### Lists

Lists gather together a collection of objects under one name:

```r
together.list <- list(
                      vector.example = dataset.a, 
                      matrix.example = matrix.example,
                      data.frame.example = people.summary
                      )
together.list
```

There are several ways to get elements of a list:

```r

together.list$matrix.example
together.list$matrix.example[,3]
together.list["matrix.example"]
together.list[["matrix.example"]]
together.list[["matrix.example"]][,2]
```

### Reading Data In

We use `read.data` or `read.csv` to read in data:

```r
gene_example <- read.csv("Gene_R_Example.csv")
```

In RStudio, we can use the "File" navigation window instead.  Navigate to the directory containing the Gene_R_example.csv that we downloaded previously. Click on the file name then click "Import Dataset."  A new window appears allowing you to modify attributes of your file.  Rename the file to the object "gene_example".

Commands like `head` and `tail` also work in R. `View` will open the dataframe in a new window:

```r
head(gene_example)
View(gene_example)
```

### Basic Plotting

A very basic plot:

```r
plot(x=gene_example$Control, y=gene_example$Treated)
```

You can specify how you want your plot to look, in terms of color, shape, labels, etc.
A nicer plot:

```r
plot(x=gene_example$Control, y=gene_example$Treated,
    xlab = "Control",
    ylab = "Treated",
    cex.lab = 1.5,
    main = "A nice scatter plot",
    pch = 16,
    bty = "n",
    col = "dark blue",
    las = 1
    )

# las
# To change the axes label style, use the graphics option las (label style). This changes the orientation angle of the labels:
# 0: The default, parallel to the axis
# 1: Always horizontal
# 2: Perpendicular to the axis
# 3: Always vertical

# bty
# To change the type of box round the plot area, use the option bty (box type):
# "o": The default value draws a complete rectangle around the plot.
# "n": Draws nothing around the plot.
```

Connecting the dots:

```r
plot(x=gene_example$Control, y=gene_example$Treated,
    xlab = "Control",
    ylab = "Treated",
    cex.lab = 1.5,
    main = "A nice scatter plot",
    pch = 16,
    bty = "n",
    col = "dark blue",
    type = "b",
    las = 1
    )

# type
# The 'type' option changes what the data points are plotted with.
# "p": points
# "l": lines
# "b": both (points and lines)
# for more options, see the plot() help page
```

#### Histograms

```r
hist(gene_example$Control)

hist(gene_example$Control,
    xlab = "Expression",
    ylab = "Number of Genes",
    cex.lab = 1.5,
    main = "A nice histogram",
    col = "cyan",
    breaks = 10,
    las = 1
    )
```

#### Boxplots

```r
 boxplot(gene_example[,2:3])
```

```r
boxplot(gene_example[,2:3],
	width = c(3,1),
	col = "red",
	border = "dark blue",
	names = c("Control", "Treatment"),
	main = "My boxplot",
	notch = TRUE,
	horizontal = TRUE
	)
```

### Saving your plots as PDFs

```r
pdf("myfigure.pdf", height=10, width=6)
par(mfrow=c(2,1))

plot(x=gene_example$Control, y=gene_example$Treated,
    xlab = "Control",
    ylab = "Treated",
    cex.lab = 1.5,
    main = "A nice scatter plot",
    pch = 16,
    bty = "n",
    col = "dark blue",
	type = "b",
	las = 1
	)
	
boxplot(gene_example[,2:3],
	width = c(3,1),
	col = "red",
	border = "dark blue",
	names = c("Control", "Treatment"),
	main = "My boxplot",
	notch = TRUE,
	horizontal = TRUE
	)
  
  dev.off()
```

If you don't specify a path in the pdf() function like we did above, the pdf will be saved in your working directory. In a new internet tab, go to http://##.oicrcbw.ca/Review_Session/myfigure.pdf to view the saved image.

To save your RStudio notebook, go to "File" > "Save As..." and enter a name.