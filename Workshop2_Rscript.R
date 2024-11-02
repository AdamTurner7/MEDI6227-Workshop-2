###MEDI6227 WORKSHOP 2
###Differential Gene Expression Analysis
###R version 4.3.2
#install and load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("recount3")
install.packages("tidyverse")
library("recount3")
library("edgeR")
library("tidyverse")

#import data set
load("Workshop1_output(2).Rdata")

#Convert chr data to cat data
dim(sample.meta.data)
view(sample.meta.data)
class(sample.meta.data$Celltype)
class(sample.meta.data$Donor)
Donor <- factor(sample.meta.data$Donor)
Celltype <- factor(sample.meta.data$Celltype)
class(Celltype)
class(Donor)

#Create model/design matrix based on the different groups e.g. M1 & M2
design<-model.matrix(~Celltype)
design

#add dispersion estimates to dgelist
dgelist.filtered.norm<-estimateDisp(dgelist.filtered.norm, design=design)

#fit linear model using the design matrix and dgelist
fit <- glmQLFit(dgelist.filtered.norm, design)
head(fit$coefficients)

#statistical testing to assess differences of expression between M1 and M2 celltypes
M2.v.M1<-glmQLFTest(fit, coef="CelltypeM2")
head(M2.v.M1$table)

#create a data frame of genes with descending P-value, join gene IDs to gene names (in gene meta data)
DE <- topTags(M2.v.M1, n = Inf)
DE <- left_join(rownames_to_column(DE$table, "gene_id"), gene.meta.data)

#produce volcano plot of Log10 and -Log10 pvalues
ggplot() + geom_point(data=DE, aes(x=logFC, y=PValue))
ggplot() + geom_point(data=DE, aes(x=logFC, y=-log10(PValue)))

#set colour of data point based on FDR value 
ggplot() + 
  geom_point(data=filter(DE, FDR<0.05), aes(x=logFC, y=-log10(PValue))) +
  geom_point(data=filter(DE, FDR>=0.05), aes(x=logFC, y=-log10(PValue)), colour="grey") +
  theme_bw()

#add labels for top20 differentially expressed genes
#install ggrepel package - this helps to separate the labels from the data points making them easier to read
install.packages("ggrepel")
library(ggrepel)
ggplot() + 
  geom_point(data=filter(DE, FDR<0.05), aes(x=logFC, y=-log10(PValue))) +
  geom_point(data=filter(DE, FDR>=0.05), aes(x=logFC, y=-log10(PValue)), colour="grey") +
  geom_text_repel(data = slice_head(DE, n=20), aes(x=logFC, y=-log10(PValue), label = gene_name))

###Clustering analysis and heatmaps