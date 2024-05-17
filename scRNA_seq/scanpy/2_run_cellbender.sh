#!/usr/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --partition gpu


source /home/ypark/conda/bin/activate
conda activate cellbender

for file in unfiltered_for_cb/*_raw.h5ad; do
    cellbender remove-background \
        --cuda \
        --input $file \
        --output cellbender/$(basename $file .h5ad)_cellbender.h5ad;

done

