
library(Seurat)
library(SeuratObject)

# Load your Seurat object
clustered_filtered_atlas <- readRDS("clustered_filtered_atlas_integrated_rpca.rds")

# Function to check and correct the features alignment
check_and_correct_features <- function(seurat_object, assay_name = "RNA") {
  assay <- seurat_object@assays[[assay_name]]
  
  # Check and correct features in data
  data_features <- rownames(assay@data)
  
  # Check and correct features in var.features
  var_features <- assay@var.features
  missing_var_features <- setdiff(var_features, data_features)
  if(length(missing_var_features) > 0) {
    cat("Missing variable features:\n")
    print(missing_var_features)
    seurat_object@assays[[assay_name]]@var.features <- intersect(var_features, data_features)
  }
  
  # Check if scale.data slot exists and correct features
  if("scale.data" %in% slotNames(assay)) {
    scale_data_features <- rownames(assay@scale.data)
    missing_scale_features <- setdiff(scale_data_features, data_features)
    if(length(missing_scale_features) > 0) {
      cat("Missing features in scale.data. Removing...\n")
      seurat_object@assays[[assay_name]]@scale.data <- assay@scale.data[data_features, , drop = FALSE]
    }
  }
}

# Run the checking and correction function
check_and_correct_features(clustered_filtered_atlas)




