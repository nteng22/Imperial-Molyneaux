#!/bin/sh

## ask PBS for time (format hh:mm:ss)
#PBS -l walltime=04:00:00
## memory (per node)
#PBS -lselect=1:ncpus=4:mem=10gb

set -e

##load application module
module load anaconda3/personal
source activate qiime2-amplicon-2024.2
# THINGS TO CHANGE IN THE SCRIPTS #
  # List the runs you want to include in the analysis in the RUN variable.
      # The format will be "RUN1 RUN2 RUN3"
RUN="NorthUmbria_stool" # name of the folder with the artefacts
PROJECT="$WORK/Project9-PROFOUND/stool/test-SRA" # path to the project folder
MSQ="$WORK/Processed-runs" # contains the processed sequencing files

cd $WORK

for f in $RUN
do

cd $PROJECT

# retain samples from the manifest file, specific to your study
# this needs the sample-id and paths to the files. 
# Use previous manifest file to create a new one with the sample-id and paths to the files. 

qiime demux filter-samples \ 
  --i-demux $MSQ/"$f"/demux/"$f"_demux.qza \
  --m-metadata-file $PROJECT/Manifest-baseline.txt \
  --o-filtered-demux $PROJECT/"$f"_filtered-demux.qza

mkdir $PROJECT/fastq-NCBI

# Get fasta files from the filtered demux
# This is what you need for the NCBI submission
qiime tools export \
  --input-path $PROJECT/"$f"_filtered-demux.qza \
  --output-path $PROJECT/fastq-NCBI

# Create filtered table with the samples from the manifest file
qiime feature-table filter-features \
  --i-table $MSQ/"$f"/denoise/"$f"_table-dada2.qza \
  --m-metadata-file $PROJECT/Manifest-baseline.txt \
  --o-filtered-table $PROJECT/"$f"_filtered-table.qza \
  --p-exclude-ids

# Filter the sequences based on the filtered table generated above
qiime feature-table filter-seqs \
  --i-data $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --i-table $PROJECT/"$f"_filtered-table.qza \
  --o-filtered-data $PROJECT/"$f"_filtered-rep-seqs.qza

# Map the sequences to the Greengenes database
qiime greengenes2 non-v4-16s \
  --i-table $PROJECT/"$f"_filtered-table.qza \
  --i-sequences $PROJECT/"$f"_filtered-rep-seqs.qza \
  --i-backbone $WORK/database/2022.10.backbone.full-length.fna.qza \
  --o-mapped-table $PROJECT/filtered-table-gg2.qza \
  --o-representatives $PROJECT/filtered-rep-seqs-gg2.qza

# Classify the sequences using the Greengenes database
qiime greengenes2 taxonomy-from-table \
    --i-reference-taxonomy $WORK/database/2022.10.taxonomy.asv.nwk.qza \
    --i-table $PROJECT/filtered-table-gg2.qza \
    --o-classification $PROJECT/filtered-taxonomy-gg2.qza

# Create the midrooted tree for the sequences
qiime alignment mafft \
  --i-sequences $PROJECT/filtered-rep-seqs-gg2.qza \
  --o-alignment $PROJECT/final-aligned-seqs-gg2.qza

qiime alignment mask \
  --i-alignment $PROJECT/final-aligned-seqs-gg2.qza \
  --o-masked-alignment $PROJECT/aligned-masked-seqs.qza

qiime phylogeny fasttree \
  --i-alignment $PROJECT/aligned-masked-seqs.qza \
  --o-tree $PROJECT/unrooted-tree.qza

qiime phylogeny midpoint-root \
  --i-tree $PROJECT/unrooted-tree.qza \
  --o-rooted-tree $PROJECT/midrooted-tree.qza

done  
