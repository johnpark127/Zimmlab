#!/bin/env Rscript

# This script is for merging seurat objects which were already merged by GSE numbers
# Check 'merge_analysis_by_groups.R' script for the first merge step
# It requires ids.csv file which contains GSM numbers and GSE numbers


library(Seurat)
library(dplyr)
library(ggplot2)


# Get list of samples and their groups
sample_info <- read.csv("ids.csv")

# Get unique GSE numbers from ids.csv
GSE_numbers <- unique(sample_info$GSE_numbers)


# Loop getting .rds files by GSM numbers, merge them one by one
merge_atlas <- NULL

# Pull up each GSE numbers and merge them by groups 
for (GSE in GSE_numbers) {
	file_name <- paste0("merged_RDS/", GSE, "_merged.rds")
	merged_GSE <- readRDS(file_name)	
	
	# Merge them together
	if (is.null(merge_atlas)){
	merge_atlas <- merged_GSE
	} else {merge_atlas <- merge(merge_atlas, merged_GSE, add.cell.ids = TRUE) 
	}
	
	
	cat("Merged",GSE,"\n")
	# Store merged files in the list
	saveRDS(merge_atlas, file = "Merged_atlas.rds")
	cat("Saved merged atlas","\n")
}

