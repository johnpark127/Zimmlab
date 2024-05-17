library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)

clustered_filtered_atlas_integrated_rpca <- readRDS("clustered_filtered_atlas_integrated_rpca.rds")

umap <- DimPlot(clustered_filtered_atlas_integrated_rpca, reduction="umap", pt.size = 0.001, group.by = "GSE_numbers", raster = FALSE)
umap_GSE_split <- DimPlot(clustered_filtered_atlas_integrated_rpca, reduction="umap", pt.size = 0.001, split.by = "GSE_numbers", raster = FALSE, combine = FALSE, ncol = 5)
umap_group <- DimPlot(clustered_filtered_atlas_integrated_rpca, reduction="umap", pt.size = 0.001, group.by = "group", raster = NULL)
umap_group_split <- DimPlot(clustered_filtered_atlas_integrated_rpca, reduction="umap", pt.size = 1, split.by = "group", raster = NULL, ncol = 3)
umap_injury_split <- DimPlot(clustered_filtered_atlas_integrated_rpca, reduction="umap", pt.size = 1, split.by = "type_of_injury", raster = NULL, ncol = 3)

Trem2_pos <- FeaturePlot(clustered_filtered_atlas_integrated_rpca, features = c("Trem2"), split.by = "group")
Spp1_pos <- FeaturePlot(clustered_filtered_atlas_integrated_rpca, features = c("Spp1"), split.by = "group")
Igf1_pos <- FeaturePlot(clustered_filtered_atlas_integrated_rpca, features = c("Igf1"), split.by = "group")
Gpnmb_pos <- FeaturePlot(clustered_filtered_atlas_integrated_rpca, features = c("Gpnmb"), split.by = "group")






#umap <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, raster = FALSE, label = TRUE, label.size = 5)
#umap_GSE <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, group.by = "GSE_numbers", raster = NULL)
#umap_GSE_split <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, raster = FALSE, split.by = "GSE_numbers", combine = FALSE)
#umap_group <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, group.by = "group", raster = NULL)
#umap_group_split <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, raster = FALSE, split.by = "group", combine = FALSE) 
#umap_injury <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, group.by = "type_of_injury", raster = NULL)
#umap_injury_split <- 
#DimPlot(clustered_filtered_atlas, reduction = "umap", pt.size = 0.01, raster = FALSE, split.by = "type_of_injury", combine = FALSE)

#Trem2_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Trem2"), split.by = "group")
#Spp1_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Spp1"), split.by = "group")
#Igf1_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Igf1"), split.by = "group")
#Gpnmb_pos <- FeaturePlot(clustered_filtered_atlas, features = c("Gpnmb"), split.by = "group")


#widths <- c[1,2,3,4,5,6,7,8,9,10]



pdf("features_RPCA.pdf", width = 15, height = 5)
print(Trem2_pos)
print(Spp1_pos)
print(Igf1_pos)
print(Gpnmb_pos)
dev.off()

#pdf("umap_RPCA.pdf", width = 20, height = 20)
#print(umap)
#print(umap_GSE)
#print(umap_GSE_split)
#print(umap_group)
#dev.off()

pdf("umap_RPCA_GSE_split.pdf", width = 24, height = 20)
print(umap_GSE_split)
dev.off()


#pdf("umap_RPCA_split.pdf", width = 15, height = 6)
#print(umap_group_split)
#dev.off()

pdf("umap_RPCA_injury.pdf", width = 25, height = 5)
print(umap_injury_split)
dev.off()

#umap_GSE + ggtitle("GSE numbers")
#umap_GSE_split + guide_area() + plot_layout(ncol = 5, guides = "collect")  + ggtitle("GSE numbers")
#umap_group + ggtitle("Group")
#umap_group_split + guide_area() + plot_layout(ncol = 3, guides = "collect") + ggtitle("Groups")
#umap_injury + ggtitle("Injury")
#umap_injury_split + guide_area() + plot_layout(ncol = 5, guides = "collect") + ggtitle("Injury")

#widths <- c[1,2,3,4,5,6,7,8,9,10]

#pdf("")
#Trem2_pos + guide_area() + plot_layout(ncol = 3, guides = "collect") + ggtitle("Trem2 positive")
#Spp1_pos + guide_area() + plot_layout(ncol = 3, guides = "collect") + ggtitle("Spp1 positive")
#Igf1_pos + guide_area() + plot_layout(ncol = 3, guides = "collect") + ggtitle("Igf1 positive")
#Gpnmb_pos + guide_area() + plot_layout(ncol = 3, guides = "collect") + ggtitle("Gpnmb positive")
