library(magrittr)
library(tidyverse)
library(gplots)

#https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
#https://scienceparkstudygroup.github.io/rna-seq-lesson/08-cluster-analysis/index.html#3-visualization-of-a-of-the-de-genes

#load dataframe of abundances. rows = genes, columns = samples
hclust_data_frame <- read.csv('./mcf7_degs_with_abundances.csv', row.names=1)

#first, clustering by gene

#Scaling the abundance data by gene (transposing into columns, running scale  function
#on columns, then transposing again)
hclust_data_by_gene <- hclust_data_frame %>% t() %>% scale() %>% t()

#Generating a distance matrix to be inputted into the hierarchical clustering function using
#Euclidean method
gene_dist <- dist(hclust_data_by_gene, method='euclidean')

#Clustering using complete linkage method
gene_hclust <- hclust(gene_dist, method = "complete")

#Clustering by sample

#Transforming and then scaling by gene. however, we do-not retranspose the data
hclust_data_by_sample <- hclust_data_frame %>% t() %>% scale() 

#Generating a distance matrix by sample
sample_dist <- dist(hclust_data_by_sample, method='euclidean')

#Clustering using the complete linkage method
sample_hclust <- hclust(sample_dist, method="complete")
sampleTree = as.dendrogram(sample_hclust, method='average')

#Generate heatmap
hotmap <- heatmap.2(as.matrix(hclust_data_by_gene),scale='row',
                    col=redgreen(100),
                    trace='none', srtCol=20, adjCol=c(1, 0))
