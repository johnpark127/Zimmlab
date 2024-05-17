#!/usr/bin/bash -l
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=30G
#SBATCH --partition large_mem

# Make Seurat available
source /home/kzimmer1/conda/bin/activate
conda activate /home/kzimmer1/conda/envs/harmony

Rscript seurat_on_merged_object.R

