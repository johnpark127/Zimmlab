#!/usr/bin/bash -l
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=50G
#SBATCH --partition large_mem 

# Make Seurat available
source /home/kzimmer1/conda/bin/activate
conda activate /home/kzimmer1/conda/envs/harmony

Rscript harmonize_a_merged_and_integrated_seurat_object.R
