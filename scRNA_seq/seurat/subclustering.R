library(Seurat)
library(SeuratObject)
library(dplyr)
library(xlsx)
library(ggplot2)

clustered_filtered_atlas<- readRDS("clustered_filtered_atlas_integrated_rpca.rds")

# Subclustering
clustered_filtered_atlas <- subset(clustered_filtered_atlas, idents = c("1","3","5","9"))
DimPlot(clustered_filtered_atlas, reduction="umap", pt.size = 0.01, raster=FALSE, label = TRUE, label.size = 5)
VlnPlot(clustered_filtered_atlas, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filtering again for immune cells
clustered_filtered_atlas <- subset(clustered_filtered_atlas, subset=nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

clustered_filtered_atlas <- NormalizeData(clustered_filtered_atlas, normalization.method = "LogNormalize", scale.factor = 10000)
clustered_filtered_atlas <- FindVariableFeatures(clustered_filtered_atlas, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(clustered_filtered_atlas), 10)
all.genes <- rownames(clustered_filtered_atlas)
clustered_filtered_atlas <- ScaleData(clustered_filtered_atlas, features = all.genes)
clustered_filtered_atlas <- RunPCA(clustered_filtered_atlas, features = VariableFeatures(object = clustered_filtered_atlas))

DimPlot(clustered_filtered_atlas, reduction = "pca")

clustered_filtered_atlas <- FindNeighbors(clustered_filtered_atlas, reduction = "pca", dims = 1:10)
clustered_filtered_atlas <- FindClusters(clustered_filtered_atlas, resolution = 0.3)
clustered_filtered_atlas <- RunUMAP(clustered_filtered_atlas, dims = 1:10)


DimPlot(clustered_filtered_atlas, reduction="umap", pt.size = 1, group.by = "GSE_numbers", raster = NULL)
DimPlot(clustered_filtered_atlas, reduction="umap", pt.size = 1, split.by = "GSE_numbers", raster = NULL, combine = FALSE)
DimPlot(clustered_filtered_atlas, reduction="umap", pt.size = 1, group.by = "group", raster = NULL)
DimPlot(clustered_filtered_atlas, reduction="umap", pt.size = 1, split.by = "group", raster = NULL)

Trem2_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Trem2"), split.by = "group")
Spp1_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Spp1"), split.by = "group")
Igf1_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Igf1"), split.by = "group")
Gpnmb_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Gpnmb"), split.by = "group")

# Save immune subcluster
saveRDS(clustered_filtered_atlas, file="clustered_filtered_atlas_rpca_immune.rds")


clustered_filtered_atlas.markers <- FindAllMarkers(clustered_filtered_atlas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

clustered_filtered_atlas.markers %>% 
	group_by(cluster) %>% 
	dplyr::filter(avg_log2FC > 1) %>%
	slice_head(n = 10) %>%
	ungroup-> top10

DoHeatmap(clustered_filtered_atlas, features = top10$gene) + NoLegend()


write.xlsx(clustered_filtered_atlas.markers, file = "RPCA_immune_sub_gene_lists.xlsx", sheetName = "immune_sub_marker_genes", col.names = TRUE, row.names = TRUE, showNA = TRUE)
write.xlsx(top10, file = "immune_sub_top10_genes.xlsx", sheetName = "RPCA_immune_sub_top10_genes", col.names = TRUE, row.names = TRUE, showNA = TRUE)





