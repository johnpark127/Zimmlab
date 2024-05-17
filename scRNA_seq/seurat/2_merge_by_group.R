#!/bin/env Rscript

# This script is for merging each individual seurat objects by GSE numbers
# Each seurat object was filtered by nFeature_RNA and percent.mt then saved as .rds file
# It requires ids.csv file which contains GSM numbers and GSE numbers

library(Seurat)

# Get list of samples and their groups
sample_info <- read.csv("ids.csv")

# Get unique GSE numbers from ids.csv
GSE_numbers <- unique(sample_info$GSE_numbers)

# Loop getting .rds files by GSM numbers, merge them one by one
merged_seurat <- list()

# Pull up each GSE numbers and merge them by groups 
for (GSE in GSE_numbers) {
	cat("Merging: ", GSE, "\n")
	
	# Subset sample_info for the current GSE; only samples having matched GSE numbers are kept in merged_GSE
	merged_GSE <- sample_info[sample_info$GSE_numbers == GSE, ] 
	
	
	# Initialize merged_seurat for the current GSE
	merged <- NULL
	
	# Another loop getting each row of merged_GSE subset and read corresponding filtered .rds files
	for (i in 1:nrow(merged_GSE)) {
		file_name <- paste0("filtered_RDS/", merged_GSE$GSM[i],"_mt_200_to_3K_features_filtered_obj.rds")
		filtered_rds <- readRDS(file_name)	
	
		# Merge them together
		if (is.null(merged)){
		merged <- filtered_rds
		} else {
			merged <- merge(merged, filtered_rds, add.cell.ids = TRUE) 
		}
	}

	# Store merged files in the list
	saveRDS(merged, file = paste0(GSE, "_merged.rds"))
	cat("Merged and saved: ", GSE, "\n")
}



