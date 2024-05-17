#!/bin/env Rscript


library(Seurat)
library(dplyr)

GSE155512_clustered_filtered <- readRDS("GSE155512_clustered_filtered.rds")
GSE159677_clustered_filtered <- readRDS("GSE159677_clustered_filtered.rds")
GSE224273_clustered_filtered <- readRDS("GSE224273_clustered_filtered.rds")

GSE155512_clustered_filtered <- JoinLayers(GSE155512_clustered_filtered)
GSE159677_clustered_filtered <- JoinLayers(GSE159677_clustered_filtered)
GSE224273_clustered_filtered <- JoinLayers(GSE224273_clustered_filtered)



GSE155512_clustered_filtered$dataset <- "GSE155512"
GSE159677_clustered_filtered$dataset <- "GSE159677"
GSE224273_clustered_filtered$dataset <- "GSE224273"


all_human_data <- merge(GSE155512_clustered_filtered, y = c(GSE159677_clustered_filtered, GSE224273_clustered_filtered), add.cell.ids = c("GSE155", "GSE159", "GSE224"))


all_human_data <- NormalizeData(all_human_data, normalization.method = "LogNormalize", scale.factor = 10000)
all_human_data <- FindVariableFeatures(all_human_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all_human_data)
all_human_data <- ScaleData(all_human_data, features = all.genes)
all_human_data <- RunPCA(all_human_data, features = VariableFeatures(object = all_human_data))
all_human_data <- FindNeighbors(all_human_data, dims = 1:10)


all_human_data <- FindClusters(all_human_data, resolution = 0.3)
all_human_data <- RunUMAP(all_human_data, dims = 1:10)


DimPlot(all_human_data, reduction = "umap", pt.size = 3, raster = FALSE)

saveRDS(all_human_data, file = "all_human_data.rds")


