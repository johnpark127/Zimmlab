library(Seurat)
library(SeuratObject)
library(dplyr)
library(xlsx)
library(ggplot2)
library(cowplot)

clustered_filtered_atlas <- readRDS("clustered_filtered_atlas_integrated_rpca.rds")

clustered_filtered_atlas.markers <- FindAllMarkers(clustered_filtered_atlas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

clustered_filtered_atlas.markers %>% 
	group_by(cluster) %>% 
	dplyr::filter(avg_log2FC > 1) %>%
	slice_head(n = 10) %>%
	ungroup-> top10

heatmap <- DoHeatmap(clustered_filtered_atlas, features = top10$gene) + NoLegend()

pdf("RPCA_heatmap.pdf", width = 15, height = 15)
print(heatmap)
dev.off()


write.xlsx(clustered_filtered_atlas.markers, file = "RPCA_gene_lists.xlsx", sheetName = "marker_genes", col.names = TRUE, row.names = TRUE, showNA = TRUE)
write.xlsx(top10, file = "top10_genes.xlsx", sheetName = "RPCA_top10_genes", col.names = TRUE, row.names = TRUE, showNA = TRUE)


