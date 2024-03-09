#!/usr/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition normal # remove this line when finished debugging
partition=normal
memory=20 # 16G+ for use with GRCh38-2020-A (human reference)
threads=5 # about 40% faster than 2 cores. Cost of using 5 cores
          # included with use of 20GB of RAM (1 core per 4GB)

script_dir=$(dirname $0)

module load parallel

# Make cellranger available
MODULEPATH="/home/molecules/module:$MODULEPATH"
module load cellranger

# transcriptome reference
transcriptome=/ourdisk/hpc/zimmlab/kzimmer1/dont_archive/ref/refdata-gex-mm10-2020-A/
fastq_dir=/scratch/ypark/GSE198621_mouse

# Run cellranger for each sample. This assumes that sample IDs match FASTQ file names
# --wait prevents this script from continuing until all the samples have been processed
cat ids.txt | parallel "sbatch --wait \
			--partition $partition \
                        --cpus-per-task $threads \
                        --mem ${memory}G \
                        `# This is the basic cellranger count command.` \
                        --wait \
                        --wrap 'cellranger count --fastqs $fastq_dir \
                           --id {} \
                           --sample {} \
                           --transcriptome $transcriptome \
                           --localcores $threads \
                           --localmem $memory \
                        ' "

# Move all of the cellranger outputs under one directory
mkdir -p outs_cellranger
cat ids.txt | parallel "mv {} outs_cellranger/{}"

# If your ids.txt file had columns like:
# sample_nameA,accessionA
# sample_nameB,accessionB
# Then do this instead:
# cat ids.txt | parallel --colsep=, "sbatch cellranger_script.sh $transcriptome $fastq_dir {1} {2}"

# If your ids.txt file refers to multiple FASTQ files per sample:
# sample_nameA<TAB>ERR12,ERR13
# sample_nameB<TAB>ERR45,ERR46
# Then do this:
# cat ids.txt | parallel --colsep=$'\t' "sbatch cellranger_script.sh $transcriptome $fastq_dir {1} {2}" 
