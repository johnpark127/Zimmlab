library(Seurat)
library(SeuratObject)
library(dplyr) # Provide pipe operator %>%. See notes at bottom of this script.
library(harmony)
library(Rcpp)

#Load data
clustered_filtered_atlas_harmonized_immune_harmonized_immune <- readRDS("clustered_filtered_atlas_harmonized_immune.rds")

#Filter out low quality cells
#clustered_filtered_atlas_harmonized_immune_harmonized_immune[["percent.mt"]] <- PercentageFeatureSet(clustered_filtered_atlas_harmonized_immune, pattern = "^MT-")

#clustered_filtered_atlas_harmonized_immune_harmonized_immune <- subset(clustered_filtered_atlas_harmonized_immune, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)


#split data by variable of choice [in meta.data]
clustered_filtered_atlas_harmonized_immune_harmonized_immune.list <- SplitObject(clustered_filtered_atlas_harmonized_immune, split.by = "GSE_numbers")

clustered_filtered_atlas_harmonized_immune_harmonized_immune.list <- lapply(X = clustered_filtered_atlas_harmonized_immune.list, FUN = function(x) {
     x <- NormalizeData(x)
     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
 })


#packageVersion("harmony")

#clustered_filtered_atlas_harmonized_immune_harmonized_immune <- IntegrateLayers(object = clustered_filtered_atlas_harmonized_immune, method = HarmonyIntegration, orig.reduction ="pca",
#				new.reduction = 'harmony', verbose = FALSE)

#Run harmony
clustered_filtered_atlas_harmonized_immune_harmonized_immune <- RunHarmony(clustered_filtered_atlas_harmonized_immune, 
				group.by.vars = c("GSE_numbers"), 
				reduction = "pca", reduction.save = "harmony")

clustered_filtered_atlas_harmonized_immune_harmonized_immune <- RunUMAP(clustered_filtered_atlas_harmonized_immune, reduction = "harmony", dims = 1:20)
clustered_filtered_atlas_harmonized_immune_harmonized_immune <- FindNeighbors(object = clustered_filtered_atlas_harmonized_immune, reduction = "harmony")
clustered_filtered_atlas_harmonized_immune_harmonized_immune <- FindClusters(clustered_filtered_atlas_harmonized_immune, resolution = 0.2)


#plot data
DimPlot(clustered_filtered_atlas_harmonized_immune_harmonized_immune, reduction = "umap", pt.size = 1, raster = FALSE, label = TRUE, label.size = 10)
DimPlot(clustered_filtered_atlas_harmonized_immune_harmonized_immune, reduction = "umap", pt.size = 1, group.by = "GSE_numbers", raster = FALSE)
DimPlot(clustered_filtered_atlas_harmonized_immune_harmonized_immune, reduction = "umap", pt.size = 1, group.by = "group", raster = FALSE)
DimPlot(clustered_filtered_atlas_harmonized_immune_harmonized_immune, reduction = "umap", pt.size = 1, group.by = "type_of_injury", raster = FALSE)

#clustered_filtered_atlas_harmonized_immune_harmonized_immune.markers <- FindAllMarkers(clustered_filtered_atlas_harmonized_immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#clustered_filtered_atlas_harmonized_immune_harmonized_immune.markers %>%
 #   group_by(cluster) %>%
  #  top_n(n = 5, wt = avg_log2FC) -> top5
#DoHeatmap(clustered_filtered_atlas_harmonized_immune_harmonized_immune, features = top5$gene) + NoLegend()



#Save harmonized data
saveRDS(clustered_filtered_atlas_harmonized_immune_harmonized_immune, file = "clustered_filtered_atlas_harmonized_immune_harmonized.rds")






