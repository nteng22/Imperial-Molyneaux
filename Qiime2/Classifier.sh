#!/bin/sh

## ask PBS for time (format hh:mm:ss)
#PBS -l walltime=02:00:00
## memory (per node)
#PBS -lselect=1:ncpus=32:mem=32gb

##load application module
module load anaconda3/personal
source activate qiime2-amplicon-2023.9

# Need to clone the repository to get the data, uncomment line 14 if this still neds to be done.
    # https://github.com/BenKaehler/readytowear
# git clone https://github.com/BenKaehler/readytowear.git 

# # Follow instructions in the READ ME file to get the ref-seqs data
# qiime tools import \
#     --input-path $WORK/database/readytowear/data/gg_13_8/full_length/gg_13_8_otus/rep_set/99_otus.fasta \
#     --type FeatureData[Sequence] \
#     --output-path $WORK/database/readytowear/data/gg_13_8/full_length/ref-seqs.qza

# # Now train the classifier to get weighted taxonomies based on the enivronment i.e., human oral microbiota. 
# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads $WORK/database/readytowear/data/gg_13_8/515f-806r/ref-seqs.qza \
#   --i-reference-taxonomy $WORK/database/readytowear/data/gg_13_8/515f-806r/ref-tax.qza \
#   --i-class-weight $WORK/database/readytowear/data/gg_13_8/515f-806r/human-oral.qza \
#   --o-classifier $WORK/database/human-oral_classifier.qza

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $WORK/database/99_otus.fasta \
  --output-path $WORK/database/99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $WORK/database/99_otu_taxonomy.txt \
  --output-path $WORK/database/99_ref-taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences $WORK/database/99_otus.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 0 \
  --p-min-length 30 \
  --o-reads $WORK/database/99_ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $WORK/database/99_ref-seqs.qza \
  --i-reference-taxonomy $WORK/database/99_ref-taxonomy.qza \
  --o-classifier $WORK/database/99_classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier $WORK/database/99_classifier.qza \
  --i-reads $WORK/database/99_ref-seqs.qza \
  --o-classification $WORK/database/99_taxonomy.qza