#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=10:mem=20gb

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate picrust2

PROJECT="$WORK/Project9-PROFOUND/URT/Merged-runs-nonv4-16S"

picrust2_pipeline.py \
    -s $PROJECT/filtered-metadata-table-sequences/dna-sequences.fasta \
    -i $PROJECT/filtered-table-gg2.biom \
    -o $PROJECT/picrust2-output \
    -p 10 \
    --stratified