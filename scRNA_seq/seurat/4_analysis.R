#!/bin/env Rscript

# It requires ids.csv file which contains GSM numbers and GSE numbers
# From this script, you will get:
#				1) .rds files for: (i) merged atlas (ii) clustered and filtered atlas
#				2) .pdf files for: (i) plots (ii)

library(Seurat)
library(dplyr)
library(ggplot2)


merge_atlas <- readRDS("Merged_atlas.rds")
merge_atlas <- JoinLayers(merge_atlas)

# Pipe merged data through Seurat functions
clustered <-
	merge_atlas %>%
	NormalizeData() %>%
	FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
	ScaleData(features = row.names(.)) %>%
	RunPCA(features = VariableFeatures(.)) %>% # dim reduction, see "PCA NOTE"
	JackStraw(num.replicate = 100) %>%
	ScoreJackStraw(dims = 1:20) %>%
	FindNeighbors(reduction = "pca", dims = 1:10) %>%
	FindClusters(resolution = 0.5) %>%
	RunUMAP(reduction = "pca", dims = 1:10)


# Save integrated, clustered data
saveRDS(clustered, "clustered_filtered_atlas.rds")
cat("Saved clustered and filtered atlas","\n")


clustered.markers <- FindAllMarkers(clustered,
                        only.pos = TRUE,
                        min.pct = 0.25,
                        logfc.threshold = 0.25)

cat("Finding gene expression markers","\n")

top5 <- clustered.markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC)
cat("Top 5 genes found","\n")

# Open PDF device to write to
pdf("plots.pdf")

	# Jack Straw plot page
	print(JackStrawPlot(clustered, dims = 1:20))

	# UMAP of clusters
	print(DimPlot(clustered,
		reduction = "umap",
		group.by = "seurat_clusters",
 		label = TRUE,
		repel = TRUE)
	)

	# UMAP of samples
	print(DimPlot(clustered,
		reduction = "umap",
		group.by = "sample")
	)

	print(DimPlot(clustered,
		reduction = "umap",
	group.by = "group")
	)


	print(DoHeatmap(clustered, features = top5$gene) + NoLegend())

	print(table(clustered$group))
	print(table(clustered$sample))
	print(table(Idents(clustered)))
	print(table(Idents(clustered), clustered$group))
	print(table(Idents(clustered), clustered$sample))

	for (i in c(10, 20)) {
		topN <- head(VariableFeatures(clustered), i)

		print(
			DotPlot(clustered,
			features = topN) +
			RotatedAxis()
		)
	}

dev.off()

top20 <- head(VariableFeatures(clustered), 20)

pdf("DotPlot_long.pdf", width = 8, height = 16)
	print(
	DotPlot(clustered,
		features = top20,
		cols     = c("blue", "red", "purple", "orange"),
		split.by = "sample"
	) +
	RotatedAxis()
	)
dev.off()





