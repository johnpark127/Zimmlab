library(Seurat)
library(SeuratObject)
library(dplyr) # Provide pipe operator %>%. See notes at bottom of this script.
library(harmony)
library(Rcpp)

#Load data
clustered_filtered_atlas <- readRDS("clustered_filtered_atlas.rds")

#Filter out low quality cells
#clustered_filtered_atlas[["percent.mt"]] <- PercentageFeatureSet(clustered_filtered_atlas, pattern = "^MT-")

#clustered_filtered_atlas <- subset(clustered_filtered_atlas, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)


#split data by variable of choice [in meta.data]
clustered_filtered_atlas.list <- SplitObject(clustered_filtered_atlas, split.by = "GSE_numbers")

clustered_filtered_atlas.list <- lapply(X = clustered_filtered_atlas.list, FUN = function(x) {
     x <- NormalizeData(x)
     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
 })


#packageVersion("harmony")

#clustered_filtered_atlas <- IntegrateLayers(object = clustered_filtered_atlas, method = HarmonyIntegration, orig.reduction ="pca",
#				new.reduction = 'harmony', verbose = FALSE)

#Run harmony
clustered_filtered_atlas <- RunHarmony(clustered_filtered_atlas, 
				group.by.vars = c("GSE_numbers"), 
				reduction = "pca", reduction.save = "harmony")

clustered_filtered_atlas <- RunUMAP(clustered_filtered_atlas, reduction = "harmony", dims = 1:20)
clustered_filtered_atlas <- FindNeighbors(object = clustered_filtered_atlas, reduction = "harmony")
clustered_filtered_atlas <- FindClusters(clustered_filtered_atlas, resolution = 0.2)


#plot data
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, raster = FALSE, label = TRUE, label.size = 10)
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, group.by = "GSE_numbers", raster = FALSE)
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, group.by = "group", raster = FALSE)
DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 1, group.by = "type_of_injury", raster = FALSE)

#clustered_filtered_atlas.markers <- FindAllMarkers(clustered_filtered_atlas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#clustered_filtered_atlas.markers %>%
 #   group_by(cluster) %>%
  #  top_n(n = 5, wt = avg_log2FC) -> top5
#DoHeatmap(clustered_filtered_atlas, features = top5$gene) + NoLegend()



#Save harmonized data
saveRDS(clustered_filtered_atlas, file = "clustered_filtered_atlas_harmonized.rds")






