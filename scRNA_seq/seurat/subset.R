library(Seurat)
library(SeuratObject)
library(dplyr)
library(xlsx)
library(ggplot2)

clustered_filtered_atlas<- readRDS("clustered_filtered_atlas_integrated_rpca.rds")

# Subclustering
clustered_filtered_atlas <- subset(clustered_filtered_atlas, idents = c("1","3","5","9"))


saveRDS(clustered_filtered_atlas, file="rpca_immune_subset.rds")


