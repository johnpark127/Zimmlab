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

