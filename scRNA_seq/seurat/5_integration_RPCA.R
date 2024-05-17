library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)

#Load data
clustered_filtered_atlas <- readRDS("clustered_filtered_atlas.rds")


#join layers as needed
clustered_filtered_atlas <- JoinLayers(clustered_filtered_atlas)


#split data by variable of choice [in meta.data]
clustered_filtered_atlas.list <- SplitObject(clustered_filtered_atlas, split.by = "GSE_numbers")

clustered_filtered_atlas.list <- lapply(X = clustered_filtered_atlas.list, FUN = function(x) {
     x <- NormalizeData(x)
     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
 })

#select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = clustered_filtered_atlas.list)
clustered_filtered_atlas.list <-lapply(X = clustered_filtered_atlas.list, FUN = function(x) {
	x <- ScaleData(x, features = features, verbose = FALSE)
	x <- RunPCA(x, features = features, verbose = FALSE)
})


atlas.anchors <- FindIntegrationAnchors(object.list = clustered_filtered_atlas.list, anchor.features = features, reduction = "rpca")

# creating integrated data assay
atlas.combined <- IntegrateData(anchorset = atlas.anchors)

#perform analysis on integrated seurat object
DefaultAssay(atlas.combined) <- "integrated"
atlas.combined <- ScaleData(atlas.combined, verbose = FALSE)
atlas.combined <- RunPCA(atlas.combined, npcs = 30, verbose = FALSE)
atlas.combined <- RunUMAP(atlas.combined, reduction = "pca", dims = 1:30)
atlas.combined <- FindNeighbors(atlas.combined, reduction = "pca", dims = 1:30)
atlas.combined <- FindClusters(atlas.combined, resolution = 0.1)

#plot data
DimPlot(atlas.combined, reduction = "umap", pt.size = 0.01, group.by = "group", raster = FALSE)
DimPlot(atlas.combined, reduction = "umap", pt.size = 0.01, raster = FALSE, label = TRUE, label.size = 10)
DimPlot(atlas.combined, reduction = "umap", pt.size = 0.01, group.by = "GSE_numbers", raster = FALSE)

#Save data
saveRDS(atlas.combined, file = "clustered_filtered_atlas_integrated_rpca.rds")








