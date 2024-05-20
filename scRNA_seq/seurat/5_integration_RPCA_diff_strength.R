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


# k.anchor <- strength of integration. Default is 5

atlas.anchors10 <- FindIntegrationAnchors(object.list = clustered_filtered_atlas.list, anchor.features = features, reduction = "rpca",
	k.anchor = 10)

atlas.anchors20 <- FindIntegrationAnchors(object.list = clustered_filtered_atlas.list, anchor.features = features, reduction = "rpca",
	k.anchor = 20)

# creating integrated data assay
atlas.combined10 <- IntegrateData(anchorset = atlas.anchors10)
atlas.combined20 <- IntegrateData(anchorset = atlas.anchors20)

#perform analysis on integrated seurat object
DefaultAssay(atlas.combined10) <- "integrated"
atlas.combined10 <- ScaleData(atlas.combined10, verbose = FALSE)
atlas.combined10 <- RunPCA(atlas.combined10, npcs = 30, verbose = FALSE)
atlas.combined10 <- RunUMAP(atlas.combined10, reduction = "pca", dims = 1:30)
atlas.combined10 <- FindNeighbors(atlas.combined10, reduction = "pca", dims = 1:30)
atlas.combined10 <- FindClusters(atlas.combined10, resolution = 0.1)

DefaultAssay(atlas.combined20) <- "integrated"
atlas.combined20 <- ScaleData(atlas.combined20, verbose = FALSE)
atlas.combined20 <- RunPCA(atlas.combined20, npcs = 30, verbose = FALSE)
atlas.combined20 <- RunUMAP(atlas.combined20, reduction = "pca", dims = 1:30)
atlas.combined20 <- FindNeighbors(atlas.combined20, reduction = "pca", dims = 1:30)
atlas.combined20 <- FindClusters(atlas.combined20, resolution = 0.1)



#plot data
DimPlot(atlas.combined10, reduction = "umap", pt.size = 0.01, group.by = "group", raster = FALSE)
DimPlot(atlas.combined20, reduction = "umap", pt.size = 0.01, group.by = "group", raster = FALSE)
DimPlot(atlas.combined10, reduction = "umap", pt.size = 0.01, raster = FALSE, label = TRUE, label.size = 10)
DimPlot(atlas.combined20, reduction = "umap", pt.size = 0.01, raster = FALSE, label = TRUE, label.size = 10)
DimPlot(atlas.combined10, reduction = "umap", pt.size = 0.01, group.by = "GSE_numbers", raster = FALSE)
DimPlot(atlas.combined20, reduction = "umap", pt.size = 0.01, group.by = "GSE_numbers", raster = FALSE)

#Save data
saveRDS(atlas.combined10, file = "clustered_filtered_atlas_integrated_rpca_10.rds")
saveRDS(atlas.combined20, file = "clustered_filtered_atlas_integrated_rpca_20.rds")



