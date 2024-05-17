library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)

#Load data
clustered_filtered_atlas <- readRDS("clustered_filtered_atlas.rds")

#clustered_filtered_atlas

#clustered_filtered_atlas <- subset(clustered_filtered_atlas, downsample=1000)

#join layers as needed
clustered_filtered_atlas <- JoinLayers(clustered_filtered_atlas)


#split data by variable of choice [in meta.data]
clustered_filtered_atlas.list <- SplitObject(clustered_filtered_atlas, split.by = "GSE_numbers")

clustered_filtered_atlas.list <- lapply(X = clustered_filtered_atlas.list, FUN = function(x) {
     x <- NormalizeData(x)
     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
 })

#integrate data using anchors
features <- SelectIntegrationFeatures(object.list = clustered_filtered_atlas.list)
spatialandsinglecell.anchors <- FindIntegrationAnchors(object.list = clustered_filtered_atlas.list, anchor.features = features)
clustered_filtered_atlas <- IntegrateData(anchorset = spatialandsinglecell.anchors)

#perform analysis on integrated seurat object
DefaultAssay(clustered_filtered_atlas) <- "integrated"
clustered_filtered_atlas <- ScaleData(clustered_filtered_atlas, verbose = FALSE)
clustered_filtered_atlas <- RunPCA(clustered_filtered_atlas, npcs = 30, verbose = FALSE)
clustered_filtered_atlas <- RunUMAP(clustered_filtered_atlas, reduction = "pca", dims = 1:30)
clustered_filtered_atlas <- FindNeighbors(clustered_filtered_atlas, reduction = "pca", dims = 1:30)
clustered_filtered_atlas <- FindClusters(clustered_filtered_atlas, resolution = 0.1)

#plot data
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, group.by = "group", raster = FALSE)
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, raster = FALSE, label = TRUE, label.size = 10)
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, group.by = "GSE_numbers", raster = FALSE)

#Save data
saveRDS(clustered_filtered_atlas, file = "clustered_filtered_atlas_integrated.rds")








