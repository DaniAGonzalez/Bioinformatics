# Script to Manipulating gene expression data GSE18394

# Setting WD
# setwd("~/OneDrive/5.2024/2.Bioinformagician")

# Loading the libraries
install.packages("dplyr")
install.packages("tidyverse")
install.packages("GEOquery")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(dplyr)
library(tidyverse)
library(GEOquery)

#Reading the data
data <-read.csv(file= "~/OneDrive/5.2024/2.Bioinformagician/GSE183947_fpkm.csv")
dim(data)

#Retrieving the metadata
gse<-getGEO(GEO='GSE183947', GSEMatrix = TRUE)
gse

# pData to get the phenodata
metadata<-pData(phenoData(gse[[1]]))
head(metadata)

# Retrieving the columns from the dataframe
metadata.subset<-select(metadata,c(1, 10, 11, 17))


# Additional operations to the dataset
metadata %>%
  select(1,10,11,17) %>%
  head()

# Renaming columns
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue=characteristics_ch1) %>%
  rename(metastasis=characteristics_ch1.1)%>%
  mutate(tissue=gsub("tissue: ", " ",tissue))%>%
  mutate(metastasis=gsub('metastasis: ', '',metastasis))

# Reshaping data 
data.long<-data %>%
  rename(gene=X)%>%
  gather(key="samples", value= "FPKM", -gene)

# Join the df data.long+metadata.modified
data.long <-data.long%>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# Explore data in datalong
data.long%>%
  filter(gene=="BRCA1"| gene == "BRCA2" ) %>%
  group_by(gene, tissue)%>%
  summarize(mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)


## Script to visualize gene expression (RNA-seq)
# Load libraries
library(tidyverse)
library(ggplot2)

# Baisc format for ggplot2
# aes= axis, geo = types of plots
# ggplot(dataframe, aes(x = variable, y= variable)) + geo
data.long %>%
  filter(gene=='BRCA1')%>%
  ggplot(.,aes(x=samples, y= FPKM, fill= tissue )) + geom_col() #with the dot as paramater
#we are saying "use the previous output"

# Density plot
# alpha value to deal with opaque
data.long %>%
  filter(gene=='BRCA1')%>%
  ggplot(., aes(x=FPKM, fill = tissue))+
  geom_density(alpha = 0.5)

# Comparing the gene expression between samples w/ different metastasis, 2 groups
data.long %>%
  filter(gene=='BRCA1')%>%
  ggplot(., aes(x=metastasis, y = FPKM))+
  geom_violin()

# Scatterplot: comparing the expression of 2 genes
# brca1 brca2
data.long %>%
  filter(gene=='BRCA1'|gene=='BRCA2' )%>%
  # Changing the shape of tha data to do the scatterplot
  spread(key=gene, value=FPKM)%>% #will return a df with the gene of interest as columns
  ggplot(., aes(x=BRCA1, y = BRCA2))+
  geom_point()
# The scatter plot shows a positive trend between both gene expression
# But it must be apply a corr test to see if it is indeed and try to fit a line
data.long %>%
  filter(gene=='BRCA1'|gene=='BRCA2' )%>%
  # Changing the shape of tha data to do the scatterplot
  spread(key=gene, value=FPKM)%>% #will return a df with the gene of interest as columns
  ggplot(., aes(x=BRCA1, y = BRCA2, color = tissue))+
  geom_point()+
  geom_smooth(method='lm', se = FALSE)




# Heatmap: to compare expression of multiple genes between multiple samples
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
  p<-data.long %>%
  filter(gene%in%genes.of.interest)%>%
  ggplot(., aes(x=samples, y= gene, fill=FPKM))+
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

#Saving the plot
  ggsave(p, filename = 'heatmap_save1.png', width = 10, height = 8)









