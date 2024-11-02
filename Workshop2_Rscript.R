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
BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")

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

#shorten name of file to cpm
 cpm <- cpm.filtered.norm
 head(cpm)
 
 #rename colnames of cpm to their respective celltype and donor
 sample.meta.data<-sample.meta.data |> 
  mutate("sample_name"=paste0(Celltype, "_Donor", Donor))
 colnames(cpm) <- sample.meta.data$sample_name
 
 #set rownames of cpm as gene ID
 filtered.gene.meta.data<-left_join(data.frame("gene_id"=rownames(cpm)), gene.meta.data)
 head(filtered.gene.meta.data)
 rownames(cpm)<-filtered.gene.meta.data$gene_name
 
 #z score scaling on genes
 z.scaled.genes <- t(cpm) |> 
   scale() |> 
   t()
 
 #find euclidian distance between samples
 sample.scaled_distances <- dist(t(z.scaled.genes), method ="euclidean")
 
# Cluster samples using hierarchecal clustering
 sample.scaled.hclust <- hclust(sample.scaled_distances, method = "complete")
plot(sample.scaled.hclust) 

#cluster gene-wise scaled cpm values
gene_distances<-dist(z.scaled.genes, method="euclidean")
gene_hclust<- hclust(gene_distances, method = "average")
plot(gene_hclust, labels = FALSE)

#cut tree into 8 clusters
clusters.genes.k8<-cutree(gene_hclust, k=8)
head(clusters.genes.k8)
table(clusters.genes.k8)

#save rownames of cluster 3 to CSV file
z.scaled.genes.cluster3<-z.scaled.genes[clusters.genes.k8==3, ]
dim(z.scaled.genes)
dim(z.scaled.genes.cluster3)
write.csv(rownames(z.scaled.genes.cluster3), file="cluster3_genenames.csv")

#plot heatmap
Heatmap(matrix=z.scaled.genes.cluster3, cluster_rows=FALSE, cluster_columns = FALSE, show_row_names = FALSE)

Heatmap(matrix=z.scaled.genes.cluster3, 
        cluster_rows=TRUE, 
        cluster_columns = TRUE, 
        show_row_names = FALSE)

#export to png

png(filename="Cluster3_z.score_heatmap.png", height=30, width=10, units="cm", res=200)
ht<-Heatmap(matrix=z.scaled.genes.cluster3, 
            cluster_rows=TRUE, 
            cluster_columns = TRUE, 
            show_row_names = FALSE,
            height=unit(20, "cm"))

draw(ht)
dev.off()
