#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=1:mem=2gb

# Need to export the feature table from QIIME 2 to a biom file. 
# Can't run this in the same script as the Picrust2 pipeline because it requires a different conda environment.
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate qiime2-amplicon-2024.2

PROJECT="$HOME/Project9-PROFOUND/URT/Merged-runs-nonv4-16S"

# Export feature table to a directory
qiime tools export \
    --input-path $PROJECT/filtered-table-gg2.qza \
    --output-path $PROJECT/exported-table

# (Optional) Copy or move the exported biom file for the next step
cp $PROJECT/exported-table/feature-table.biom $PROJECT/filtered-table-gg2.biom