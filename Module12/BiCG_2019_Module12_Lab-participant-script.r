
install.packages(pkgs = c('SNFtool','RColorBrewer','reshape2','ggplot2',"survival",'rms',"repmis"))

# Load libraries
library('SNFtool')
library('RColorBrewer')
library('reshape2')
library('ggplot2')
library("survival")
library('rms')  
library('repmis')

# My working directories (one for my laptop and one for my lab computer), enter your own here!
# setwd("C:/Users/Owner/Desktop/Goldenberg Lab/CBW/Genomic Medicine Workshop/")
setwd("<FILL IN WITH YOUR OWN WORKING DIRECTORY>")

## Load data using AWS: 
source_data('https://github.com/bioinformaticsdotca/BiCG_2017/raw/master/module9/CBW-CancerGenomics-June2017-Lab-Data-Half1.RData')
source_data('https://github.com/bioinformaticsdotca/BiCG_2017/raw/master/module9/CBW-CancerGenomics-June2017-Lab-Data-Half2.RData')

# load('CBW-CancerGenomics-June2017-Lab-Data-Half1.RData') ## This is how you would import it on on your local machine
# load('CBW-CancerGenomics-June2017-Lab-Data-Half2.RData')

## What's in your environment now? 
ls() ## You should have loaded 7 datasets: "brca.subtype.data","clinical.data","cnv","methyl","mirna","mrna","surv.data"

## Looking at dimensions of the data we're integrating: mRNA, methylation, miRNA, and cnv
dim(mrna)
dim(methyl)
dim(mirna)
dim(cnv)

## Looking at the data we're integrating (looking at first 6 rows and first 6 columns)
mrna[1:6,1:6] # mRNA
methyl[1:6,1:6] # Methylation levels averaged over genes
mirna[1:6,1:6] # miRNA 
cnv[1:6,1:6] # Copy number variation levels averaged over genes 

## Check that each dataset has the same patients
identical(names(mrna),names(methyl))
identical(names(mrna),names(mirna))
identical(names(mrna),names(cnv))

## Put datasets into a list

data_list <- list(mrna,methyl,mirna,cnv)
names(data_list) <- c("mRNA","Methyl","miRNA","CNV")

##**************************************************##
##                                                  ## 
##          Similarity Network Fusion               ##
##                                                  ##
##**************************************************##

## First, set all the parameters:
K = 20; ##number of neighbors, must be greater than 1. usually (10~30)
alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
T = 20; ###Number of Iterations, usually (10~50)

## Standard normalize each data type

data_norm_list <- lapply(<PUT THE DATA LIST HERE>,standardNormalization)

## Create distance matrices for each data type 
  # (NOTE: THIS WILL TAKE A LITTLE TIME)
data_dist_list <- lapply(X = <NORMALIZED DATA LIST HERE>,
                         function(x){dist2(t(x),t(x))})

## Create affinity matrices for each data type
data_aff_list <-  lapply(X = <DISTANCE LIST HERE>,
                         function(x){affinityMatrix(x,K,alpha)})

## Get clusters for each data type (we're doing this simply/quickly since it's just a check)
  # Estimate number of clusters
(cluster_nums <- lapply(X = <AFFINITY MATRIX LIST HERE>,
                       function(x){estimateNumberOfClustersGivenGraph(W = x,NUMC = 2:10)[[1]]}))

  # Cluster data
clusters <- sapply(1:4,function(i){spectralClustering(affinity = data_aff_list[[i]],
                                                             K = cluster_nums[[i]])})
colnames(clusters) <- c("mRNA","Methyl","miRNA","CNV")

  # Looking closer at clusters
apply(clusters,2,table)

## Create heatmaps for each data type

displayClustersWithHeatmap(W = data_aff_list[["mRNA"]], group = clusters[,"mRNA"], col = brewer.pal(name = "RdGy",n = 10))
displayClustersWithHeatmap(W = data_aff_list[["Methyl"]], group = clusters[,"Methyl"], col = brewer.pal(name = "RdGy",n = 10))
displayClustersWithHeatmap(W = data_aff_list[["miRNA"]], group = clusters[,"miRNA"], col = brewer.pal(name = "RdGy",n = 10))
displayClustersWithHeatmap(W = data_aff_list[["CNV"]], group = clusters[,"CNV"], col = brewer.pal(name = "RdGy",n = 10))

##  Running SNF algorithm
W = SNF(<AFFINITY MATRIX LIST HERE>,K,T)
  # Add column and row names to SNF affinity matrix
colnames(W) <- rownames(W) <- colnames(mirna)

## Choosing number of clusters using eigen-gaps and rotation cost algorithms 
estimateNumberOfClustersGivenGraph(W = W,2:10)
## Both eigen-gaps [[1]] and rotation cost [[2]] 
##    indicate that 2 is the best number of clusters

# Perform clustering on the fused network
clustering2 = spectralClustering(W,2)

# Look at distribution of group membership
  # With table
table(clustering2)
  # With barplot
barplot(table(clustering2),col = c('darkorchid4','dodgerblue4'))

# Create dataframe for cluster 2
cluster2.df <- data.frame(cbind(colnames(mrna),clustering2))
names(cluster2.df) <- c("id","cluster")

# Heatmap of fused matrix
displayClustersWithHeatmap(W = W,group = clustering2,
                           col = brewer.pal(name = "Spectral",n = 10))

# Can check out how clustering of 3 as well
clustering3 = spectralClustering(W,3)
displayClustersWithHeatmap(W = W,group = clustering3,col = brewer.pal(name = "RdGy",n = 10))

##******************************************************************##
##                                                                  ## 
##    CONDUCTING SURVIVAL ANALYSIS USING CLUSTERS AS PREDICTORS     ##
##                                                                  ##
##******************************************************************##

###########################
#     GENERATING SURVIVAL OUTCOME
###########################

### Now we're going to look at the survival data we have for our patients
head(surv.data) # These are the columns in our survival data

## Let's take a look at the structure of variables we're going to use
str(surv.data$patient.follow_ups.follow_up.vital_status)
str(surv.data$patient.days_to_death)
str(surv.data$patient.days_to_last_followup)
str(surv.data$patient.age_at_initial_pathologic_diagnosis)

# Make age at initial pathological diagnosis a numeric variable
surv.data$patient.age_at_initial_pathologic_diagnosis <- as.numeric(surv.data$patient.age_at_initial_pathologic_diagnosis)

surv.data$time.to.event <- surv.data$patient.days_to_last_followup
surv.data$time.to.event[is.na(surv.data$patient.days_to_death) == FALSE] <- 
  surv.data$patient.days_to_death[is.na(surv.data$patient.days_to_death) == FALSE]

surv.data$event <- NA
surv.data$event[surv.data$patient.follow_ups.follow_up.vital_status == "alive"] <- 0
surv.data$event[surv.data$patient.follow_ups.follow_up.vital_status == "dead"] <- 1

surv.data$survival.outcome <- Surv(surv.data$time.to.event,
                                   surv.data$event)

#######################
#   SUMMARIZING SURVIVAL WITHOUT COVARIATE
#######################
(survival.fit <- survfit(surv.data$survival.outcome ~ 1, 
                              conf.type = "log-log"))

#######
#     CREATING BASIC KM CURVE
#######
plot(<SURVIVAL FIT OBJECT HERE>,col="blue4")

#######
#     COMPARING SURVIVAL ACROSS GROUPS
#######

## merging SNF groups and clinic.data (recall we made 'ids.groups2' when we clustered SNF)
clinic.snf.group.df <- merge(x = surv.data, ## clinic dataframe
                               y = cluster2.df, ## SNF group dataframe
                               by="id")  ## survival and SNF group ID column

clinic.snf.group.df$cluster.fac <- factor(clinic.snf.group.df$cluster,
                                        levels = c(1,2))

###
#   TESTING THE DIFFERENCE IN SURVIVAL TIME USING PETO&PETO MODIFICATION ON THE 
#       GEHAN-WILCOXON TEST, USING THE survdiff FUNCTION
###

## Cluster Assignment
survdiff(clinic.snf.group.df$survival.outcome ~ 
           clinic.snf.group.df$cluster.fac, rho=1)

# npsurv (non-parametric survival fit) function is a work around/replacement 
#   for survfit since survfit no longer works with survplot which we want to use below

## looking at marginal survival difference by SNF generated cluster

survival.fit.by.snf.group <- npsurv(survival.outcome ~ cluster.fac,
                                         data = clinic.snf.group.df)

#### 
#   ADDING CONFIDENCE BOUNDS AND COLORS TO KM CURVES PLOTTED BY GROUP
####  

survplot(fit = survival.fit.by.snf.group,col=c('forestgreen','darkorchid4'),lwd=2.5,
         col.fill = sapply(c('forestgreen','darkorchid4'),
                           function(x){adjustcolor(x, alpha.f = 0.4)}),
         xlab="Days to Event")

###########
#   Cox proportional hazards analysis
###########

# Want to fit model with survival as an outcome

## Efron method
(coxph.fit <- coxph(survival.outcome ~ 
                      cluster.fac + patient.age_at_initial_pathologic_diagnosis, 
                    data = clinic.snf.group.df,method = "efron"))
summary(coxph.fit)

## Exact method
(coxph.fit <- coxph(survival.outcome ~ 
                      cluster + patient.age_at_initial_pathologic_diagnosis, 
                    data = clinic.snf.group.df,method = "exact"))
summary(coxph.fit)

## Breslow method

(coxph.fit <- coxph(survival.outcome ~ 
                      cluster + patient.age_at_initial_pathologic_diagnosis, 
                    data = clinic.snf.group.df,method = "breslow"))
summary(coxph.fit)

## Extracting results
str(summary(coxph.fit))

(coxph.coefs <- summary(coxph.fit)$coef)
(coxph.confint <- summary(coxph.fit)$conf.int)

(coxph.results <- cbind(coxph.coefs,coxph.confint))
colnames(coxph.results)

write.csv(coxph.results[,c("coef","exp(coef)","se(coef)","Pr(>|z|)","lower .95","upper .95" )],
          file = "Cox-PH-model-results.csv")

###
#   Using cox.zph to test for covariate-specific 
#     and global proportional hazards as well as
#     plotting scho residuals to check for 
#     non-proportional hazards -- significance implies non-proportionality
###

cox.zph(fit = coxph.fit)  
par(mfrow=c(2,1))
plot(cox.zph(fit = coxph.fit))
par(mfrow=c(1,1))

### 
#   Checking for influential observations (outliers)
###

dfbeta <- residuals(coxph.fit, type = 'dfbeta') ## Dataframe of change in coefficients as each individual removed

par(mfrow=c(2,1))
for(j in 1:2){
  plot(dfbeta[,j],ylab=names(coef(coxph.fit))[j],
       pch=19,col='blue')
  abline(h=0,lty=2)
}
par(mfrow=c(1,1))
##  No extremely influential points


###
#   Checking for linearity in the covariates using plots of 
#   martingale residuals against the individual covariates
#     NOTE: This is not necessary for binary variables
#         so we only check it in our age of initial diagnosis 
#           covariate
###
martingale.resids <- residuals(coxph.fit,type = 'martingale')

head(martingale.resids)

head(names(martingale.resids))

head(as.numeric(names(martingale.resids)))

non.missing.inds <- as.numeric(names(martingale.resids)) 
## only including individuals who have all the data and therefore a residual value

par(mfrow=c(2,1))
## Null plot for residuals: 
plot(y = martingale.resids,
     x = clinic.snf.group.df$patient.age_at_initial_pathologic_diagnosis[non.missing.inds],
     ylab = 'Residuals', xlab = 'Age of Initial Diag', pch= 19, col = 'blue')
abline(h=0,lty=2,col='red')
lines(lowess(x = clinic.snf.group.df$patient.age_at_initial_pathologic_diagnosis[non.missing.inds],
             y = martingale.resids, iter = 0))

## Component-plus-residual plot

b <- coef(coxph.fit)[2]
x <- clinic.snf.group.df$patient.age_at_initial_pathologic_diagnosis[non.missing.inds]
plot(x, b*x + martingale.resids, 
     xlab='Age of Initial Diag',
     ylab="component+residual", 
     pch = 19, col = 'blue')
abline(lm(b*x + martingale.resids ~ x), lty=2, col = 'red')
lines(lowess(x, b*x + martingale.resids, iter = 0))

## deviation of lowess line from 0-line and fit slope are 
# extremely small therefore linearity seems to hold

par(mfrow=c(1,1))


#============================================================#
#************************************************************#
#                                                            #      
#                     BONUS MATERIAL                         #
#                                                            #    
#************************************************************#
#============================================================#

##*********************************************************##
##                                                         ## 
##    BONUS: ASSESSING CLINICAL ATTRIBUTES OF CLUSTERS     ##
##                                                         ##
##*********************************************************##

## 
##    LOOKING AT DISTRIBUTION OF Basal, LumA, and LumB subtypes 
##

head(brca.subtype.data)

brca.subtype.data$basal <- "Not Basal"
brca.subtype.data$basal[substr(x = brca.subtype.data$V2,start = 1,
                               stop = 9) == "ER-/HER2-"] <- "Basal"

brca.subtype.data$luma <- "Not LumA"
brca.subtype.data$luma[substr(x = brca.subtype.data$V2,start = 11,
                              stop = nchar(brca.subtype.data$V2)) == "Low Prolif"] <- "LumA"

brca.subtype.data$lumb <- "Not LumB"
brca.subtype.data$lumb[substr(x = brca.subtype.data$V2,start = 11,
                              stop = nchar(brca.subtype.data$V2)) == "High Prolif"] <- "LumB"

brca.subtype.cluster2.merge <- merge(brca.subtype.data,
                                     cluster2.df,
                                     by.x = 'V1',
                                     by.y = 'id')

table(brca.subtype.cluster2.merge$basal,
      brca.subtype.cluster2.merge$cluster)
table(brca.subtype.cluster2.merge$luma,brca.subtype.cluster2.merge$cluster)
table(brca.subtype.cluster2.merge$lumb,brca.subtype.cluster2.merge$cluster)

chisq.test(table(brca.subtype.cluster2.merge$basal,
                 brca.subtype.cluster2.merge$cluster))
chisq.test(table(brca.subtype.cluster2.merge$luma,brca.subtype.cluster2.merge$cluster))
chisq.test(table(brca.subtype.cluster2.merge$lumb,brca.subtype.cluster2.merge$cluster))

## Simple barplot
barplot(table(brca.subtype.cluster2.merge$basal,
              brca.subtype.cluster2.merge$cluster))

## Colorful barplots side-by-side
par(mfrow = c(1,3)) # this plots 1 row and 3 columns of graphs (3 graphs side by side)
barplot(table(brca.subtype.cluster2.merge$basal,brca.subtype.cluster2.merge$cluster), main = "Basal", xlab = 'cluster',
        col = c('tomato',rgb(red=255, green=99, blue=71, alpha=90, names = NULL, maxColorValue = 255)))
barplot(table(brca.subtype.cluster2.merge$luma,brca.subtype.cluster2.merge$cluster), main = "LumA", xlab = 'cluster',
        col = c('darkturquoise',rgb(red=0, green=206, blue=209, alpha=60, names = NULL, maxColorValue = 255)))
barplot(table(brca.subtype.cluster2.merge$lumb,brca.subtype.cluster2.merge$cluster), main = "LumB", xlab = 'cluster',
        col = c('seagreen4',rgb(red=46, green=139, blue=87, alpha=90, names = NULL, maxColorValue = 255)))
## cluster 2 appears to be primarily basal with some lumB

par(mfrow = c(1,1)) ## change back to single frame graph

## 
##      BONUS: CREATE PIE CHART OF SUBTYPES IN CLUSTER 2
## 

# First, create data frame of cluster 2 only
cluster1only.df <- brca.subtype.cluster2.merge[brca.subtype.cluster2.merge$cluster == 1,] 
# Have a look at it
head(cluster1only.df)
# Create new variable indicating the subtype
cluster1only.df$Subtype <- 'No subtype'
cluster1only.df$Subtype[cluster1only.df$basal == "Basal"] <- "Basal" # Use the three subtype variables to name subtypes in a single variable
cluster1only.df$Subtype[cluster1only.df$luma == "LumA"] <- "LumA"
cluster1only.df$Subtype[cluster1only.df$lumb == "LumB"] <- "LumB"

# Frquencies of each subtype in cluster 2
table(cluster1only.df$Subtype)

# Relative frequencies
table(cluster1only.df$Subtype)/sum(table(cluster1only.df$Subtype))

# Put frequencies in a pie chart
pie(table(cluster1only.df$Subtype),labels = names(table(cluster1only.df$Subtype)),col = brewer.pal("Dark2",n = 4))

## 
##    BONUS: LOOKING AT DISTRIBUTION OF OTHER CLINICAL VARIABLES  
##

# Merge our other clinical set with the cluster data frame
clinical.clust.df <- merge(clinical.data,cluster2.df,by.x = 'sample_id',by.y = 'id')
head(clinical.clust.df)

# Look at counts for different histological types
table(clinical.clust.df$histological_type,clinical.clust.df$cluster)

# Test if age of diagnosis confers with cluster membership
t.test(age_at_initial_pathologic_diagnosis ~ cluster,clinical.clust.df)
# Plot age of onset for each cluster
ggplot(clinical.clust.df,aes(age_at_initial_pathologic_diagnosis,fill = cluster)) + xlab('Age at Initial Pathological Diagnosis') + 
  geom_density(alpha = 0.6) + scale_fill_manual( values = c('darkorchid4','dodgerblue4')) + theme_bw()

t.test(tumor_nuclei_percent ~ cluster,clinical.clust.df)
ggplot(clinical.clust.df,aes(tumor_nuclei_percent,fill = cluster)) + xlab('Tumor Nuclei Percentage') + 
  geom_density(alpha = 0.6) + scale_fill_manual( values = c('darkorchid4','dodgerblue4')) + theme_bw()


##                                                   
##    BONUS: DETERMINING FEATURES THAT DRIVE CLUSTERS (this can take a bit of time)
##                                                  

  # miRNA
kw.mirna <- data.frame(matrix(nrow = nrow(mirna),ncol = 2))
names(kw.mirna) <- c('miRNA','kw_pval')

for(k in 1:nrow(mirna)){
  kw.mirna[k,'miRNA'] <- rownames(mirna)[k]
  kw.mirna[k,'kw_pval'] <- kruskal.test(unlist(mirna[k,]) ~ clustering2)$p.value
}

  # Methylation
kw.methyl <- data.frame(matrix(nrow = nrow(methyl),ncol = 2))
names(kw.methyl) <- c('probe','kw_pval')

for(k in 1:nrow(methyl)){
  kw.methyl[k,'probe'] <- rownames(methyl)[k]
  kw.methyl[k,'kw_pval'] <- kruskal.test(unlist(methyl[k,]) ~ clustering2)$p.value
}

  # mRNA
kw.mrna <- data.frame(matrix(nrow = nrow(mrna),ncol = 2))
names(kw.mrna) <- c('gene','kw_pval')

for(k in 1:nrow(mrna)){
  kw.mrna[k,'gene'] <- rownames(mrna)[k]
  kw.mrna[k,'kw_pval'] <- kruskal.test(unlist(mrna[k,]) ~ clustering2)$p.value
}

# CNV
kw.cnv <- data.frame(matrix(nrow = nrow(cnv),ncol = 2))
names(kw.cnv) <- c('gene','kw_pval')

for(k in 1:nrow(cnv)){
  kw.cnv[k,'gene'] <- rownames(cnv)[k]
  kw.cnv[k,'kw_pval'] <- kruskal.test(unlist(cnv[k,]) ~ clustering2)$p.value
}

# Put kruskal wallis data frames into a list 
kw.list <- list(kw.mrna,kw.methyl,kw.mirna,kw.cnv)
names(kw.list) <- c('mRNA','Methyl','miRNA','CNV')

# Sort kw.list
head(kw.methyl)
kw.list.sorted <- lapply(kw.list,function(x)x[order(x = x[,'kw_pval']),])

## Look at top associated genes and miRNA for each of the data sets
lapply(kw.list.sorted,head) ## Note: numeric values for mRNA gene names -- important to label rows and columns in your data! 

## 
##  BONUS: Create cytoscape file
## 
W_cyto <- W

# Set lower triangle and diagonal to missing since matrix is symmetric
W_cyto[lower.tri(W_cyto, diag = F)] <- NA

W_cyto[1:6,1:6] ## look at first six rows and columns, lower triangle and diagonal are now 'NA'

# Make matrix into a dataframe
W_cyto_mlt <- melt(W_cyto)

# Check it out -- lots of missing values
head(W_cyto_mlt)

# Let's overwrite the object with the missing values removed
W_cyto_mlt <- W_cyto_mlt[!(is.na(W_cyto_mlt$value)),]

# That's better! 
head(W_cyto_mlt)

# Write a text file with the affinities for each of your patients relative to each other
write.table(W_cyto_mlt,file = "SNF-affinity-matrix-cytoscape-infile.txt",quote=FALSE,row.names=FALSE,col.names=TRUE)
