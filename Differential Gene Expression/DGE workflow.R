## DeSeq2
# Goal: evaluating the transcripts obtained in 4 cell lines treated and not treated with a corticoid

# Set working directory
setwd("/Users/danielaalejandragonzalez/Library/CloudStorage/OneDrive-Personal/5.2024/1.Bioinformagician/7. DeSeq2")

#Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("airway")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


install.packages("tidyverse")


# Load libraries
library(DESeq2)
library(tidyverse)
library(airway)

#Retrieving the data
# script to get data from airway package

library(airway)

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)


# STEP 1: PREPARING COUNT DATA
# read in counts data
counts_data <-read.csv('counts_data.csv')
head(counts_data)


# read in sample info to know which information gives each column because we do not known by looking at the colnames
colData<- read.csv('sample_info.csv') #reading the csv file with information

# making sure the row names in colData matches to column names in counts_data
# and in the same order
#check if the row names matchs with the columns names in sample_info
all(colnames(counts_data) %in% rownames(colData))

# are in the same order in the file data vs sample_info?
all(colnames(counts_data)==rownames(colData))

# STEP 2
# Construct a DESeqDataSet object
dds <-DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design= ~dexamethasone)

dds

# Prefiltering : removing rows with low gene counts
# Keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,] #we define the filter and then applied
dds

# Set the factor level Comparison between treated and untreated
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

#NOTE : if the assay have technical replicates, then you have to collapse all of them
# Only the technical, not the biological(never collapse the last ones)

#STEP 3 : Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# Explore results 
#log2 fold change (MLE): dexamethasone treated vs untreated : always treated
#(positive are upregulated genes in the treated condition )
# Wald test p-value: dexamethasone treated vs untreated 
# basemean:average of the normalized counts taking all samples 
# lfcSE : standard error
# Take into consider the p-adjusted value for the false positives

summary(res)
res0.01<-results(dds, alpha=0.01) #adjusting p-value
summary(res0.01)

#Contrasts
resultsNames(dds)

# e.g: treated_4hs, treated_8hr, untreated : Comparing between level
#results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated" ))

#MA plot
plotMA(res)



