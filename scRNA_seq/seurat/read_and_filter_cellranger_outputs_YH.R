#!/bin/env Rscript


library(Seurat)
library(dplyr) # Provide pipe operator %>%. See notes at bottom of this script.
library(ggplot2)

# Get list of samples and their groups
sample_info <- read.csv("ids.csv")

all_samples <- list()

#{Process all of the samples individually}
for (row_name in row.names(sample_info)) {

    row <- sample_info[row_name, ]

    # Force sample name to be text, even if it looks like a number
    samp_name <- as.character(row$GSM)

    # All outs files are stored in "outs_cellranger" with GSM numbers as their name
    matrix_dir <- paste0("outs_cellranger/",
			 samp_name,
                         "/outs/filtered_feature_bc_matrix/")

    all_samples[[samp_name]] <- Read10X(matrix_dir) %>%
                                         CreateSeuratObject(
                                           counts = .,
                                           project = "human_data",
                                           min.cells = 3,
                                           min.features = 200
                                         )

    # Store all information in meta data
    all_samples[[samp_name]]$sample <- samp_name
    all_samples[[samp_name]]$group <- row$group
    all_samples[[samp_name]]$control_type <- row$control_type
    all_samples[[samp_name]]$type_of_injury <- row$type_of_injury
    all_samples[[samp_name]]$injury_dose <- as.character(row$injury_dose)
    all_samples[[samp_name]]$surgery_time <- as.character(row$surgery_time)
    all_samples[[samp_name]]$surgery_temp <- as.character(row$surgery_temp)
    all_samples[[samp_name]]$nephrectomy <- row$nephrectomy
    all_samples[[samp_name]]$treatment <- row$treatment
    all_samples[[samp_name]]$anesthetics <- row$anesthetics
    all_samples[[samp_name]]$age_at_injury <- as.character(row$age_at_injury)
    all_samples[[samp_name]]$days_post_injury <- as.character(row$days_post_injury)
    all_samples[[samp_name]]$species <- row$species
    all_samples[[samp_name]]$strain <- as.character(row$strain)
    all_samples[[samp_name]]$sex <- row$sex
    all_samples[[samp_name]]$sorting <- as.character(row$sorting)
    all_samples[[samp_name]]$unique_ann <- as.character(row$unique_ann)
    all_samples[[samp_name]]$GSE_numbers <- as.character(row$GSE_numbers)

    #{Determine mitochondrial reads per cell}
    # updated pattern to match genome annotation (e.g. Humans may be "^MT-")
    all_samples[[samp_name]][["percent.mt"]] <-
        PercentageFeatureSet(all_samples[[samp_name]], pattern = "^mt-")

    #{Save original data set, before selection}
    saveRDS(all_samples[[samp_name]], file = paste0(samp_name, "_orig_obj.rds"))

    # Adjust the following filters, if not per-sample at least per batch
    filtered_cells <- all_samples[[samp_name]] %>%
                        subset(subset = nFeature_RNA > 200 &
                          nFeature_RNA < 3000 &
                          percent.mt < 50)

     #{Save filtered cells}
     saveRDS(filtered_cells,
             file = paste0(samp_name,
                           "_mt_200_to_3K_features_filtered_obj.rds")
     )
}

# Start off with the first item as "merged". Others will be merged into it.
merged <- all_samples[[1]] # critical: use double square brackets

# Get all samples after the first
remaining_samples <- unlist(all_samples[2:length(all_samples)])

for (remaining_sample in remaining_samples) {
    merged <- merge(merged, y = remaining_sample)
}

merged <- JoinLayers(merged)

# Pipe merged data through Seurat functions
clustered <-
    merged %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(features = row.names(.)) %>%
    RunPCA(features = VariableFeatures(.)) %>% # dim reduction, see "PCA NOTE"
    JackStraw(num.replicate = 100) %>%
    ScoreJackStraw(dims = 1:20) %>%
    FindNeighbors(reduction = "pca", dims = 1:10) %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(reduction = "pca", dims = 1:10)

clustered <- JoinLayers(clustered)

# Save integrated, clustered data
saveRDS(clustered, "clustered_filtered.rds")



clustered.markers <- FindAllMarkers(clustered,
                        only.pos = TRUE,
                        min.pct = 0.25,
                        logfc.threshold = 0.25)

top5 <- clustered.markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC)

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


# PCA NOTE: Single cell gene expression profiles are so multidimensional that we
# need to reduce the dimensionality to be able to compute on them. For example,
# if you have expression data for 30,000 genes, that's 30,000 dimensions.
# Principle Component Analysis (PCA) provides a way for us to represent gene
# expression profiles in fewer dimensions. This is called dimensionality
# reduction. Once we have this lower dimensional space (PCA space), we can then
# cluster cells based on the similarity of their gene expression profiles within
# this PCA space and cluster them based on the closeness in similarity
# distances", cluster the cells. UMAP provides an imperfect two-dimensional
# representation of the position of each cell within PCA space.

# NOTE: In this sample data, there was no "integration" (a.k.a.
# batch-correction) performed because the experiment was carefully designed to
# avoid batch effects.

#
# ABOUT dplyr pipes %>%
#
# Think of the %>% operator taking what is on the left and placing it as the
# first argument on the right. For example, if we take the first line of code
# below and apply the pipes, one at a time, we get the following:
#
# 100 %>% sqrt() %>% sum(26) %>%  sqrt()
#
#      sqrt(100) %>% sum(26) %>%  sqrt()
#
#         10     %>% sum(26) %>%  sqrt()
#
#                        36  %>%  sqrt()
#
#                                sqrt(36)
#
#                                   6
#
# The following three code blocks are equivalent. Which is cleanest?
#
# # nested function calls
# a <- sqrt(sum(sqrt(100),26))
#
# # Reassign result to original variable
# a <- 100
# a <- sqrt(a)
# a <- sum(a,26)
# a <- sqrt(a)
#
# # pipes, one operation per line
# a <- 100     %>%
#      sqrt()  %>%
#      sum(26) %>%
#      sqrt()
#
# # same use of pipes, all operations in one line
# a <- 100 %>% sqrt() %>% sum(26) %>% sqrt()
