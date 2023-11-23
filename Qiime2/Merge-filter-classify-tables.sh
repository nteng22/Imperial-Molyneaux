#!/bin/sh

## ask PBS for time (format hh:mm:ss)
#PBS -l walltime=02:00:00
## memory (per node)
#PBS -lselect=1:ncpus=2:mem=2gb

##load application module
module load anaconda3/personal
source activate qiime2-amplicon-2023.9

PROJECT="$WORK/Project2-IPF"
MSQ="$WORK/Processed-runs"

cd $WORK
#Merge tables
mkdir $PROJECT/Merged-runs

  # Prints out the paths to each of the files that will be used in the command line arguments for qiime feature-table.
find $MSQ/MSQ*/quality-control/ -name "*_filtered-table.qza" -type f -exec readlink -f {} \; > $PROJECT/Merged-runs/filtered-table_paths.txt
find $MSQ/MSQ*/quality-control/ -name "*_sequence-hits.qza" -type f -exec readlink -f {} \; > $PROJECT/Merged-runs/filtered-seqs_paths.txt

# Read the paths from the file and store them in arrays
mapfile -t table_paths < $PROJECT/Merged-runs/filtered-table_paths.txt
mapfile -t seq_paths < $PROJECT/Merged-runs/filtered-seqs_paths.txt

# Build the command line arguments for qiime feature-table merge
merge=""
for path in "${table_paths[@]}"; do
  merge+="--i-tables $path "
done

# Build the command line arguments for qiime feature-table merge-seqs
MergeSeq=""
for path in "${seq_paths[@]}"; do
  MergeSeq+="--i-data $path "
done

# Run the qiime feature-table merge command with the generated arguments
qiime feature-table merge $merge --o-merged-table $PROJECT/Merged-runs/merged_table.qza --p-overlap-method sum
qiime feature-table merge-seqs $MergeSeq --o-merged-data $PROJECT/Merged-runs/merged_sequences.qza

qiime feature-table summarize \
  --i-table $PROJECT/Merged-runs/merged_table.qza \
  --o-visualization $PROJECT/Merged-runs/merged_table.qzv

qiime feature-table tabulate-seqs \
  --i-data $PROJECT/Merged-runs/merged_sequences.qza \
  --o-visualization $PROJECT/Merged-runs/merged_sequences.qzv 

# FILTER SAMPLES 
# FILTER BY SAMPLES I.E., METADATA

# Samples kept
qiime feature-table filter-samples \
  --i-table $PROJECT/Merged-runs/merged_table.qza \
  --m-metadata-file $PROJECT/Mapping/Manifest-map.txt \
  --o-filtered-table $PROJECT/Merged-runs/filtered-metadata-table.qza 

# Samples dropped
qiime feature-table filter-samples \
  --i-table $PROJECT/Merged-runs/merged_table.qza \
  --m-metadata-file $PROJECT/Mapping/Manifest-map.txt \
  --p-exclude-ids \
  --o-filtered-table $PROJECT/Merged-runs/filtered-dropped-metadata-table.qza \

# Create a feature table artifact based on the filtered metadata table.
qiime tools export \
  --input-path $PROJECT/Merged-runs/merged_sequences.qzv \
  --output-path $PROJECT/Merged-runs/filtered-metadata-table-sequences

qiime tools import \
  --input-path $PROJECT/Merged-runs/filtered-metadata-table-sequences/sequences.fasta \
  --output-path $PROJECT/Merged-runs/filtered-metadata-table-sequences.qza \
  --type 'FeatureData[Sequence]'

# Tabular summary of the sequences in a feature data object
qiime feature-table tabulate-seqs \
  --i-data $PROJECT/Merged-runs/filtered-metadata-table-sequences.qza \
  --o-visualization $PROJECT/Merged-runs/filtered-metadata-table-sequences.qzv
    
# Summary of the filtered table
qiime feature-table summarize \
  --i-table $PROJECT/Merged-runs/filtered-metadata-table.qza \
  --o-visualization $PROJECT/Merged-runs/merged-filtered-metadata-table.qzv 

# Summary of the dropped table
qiime feature-table summarize \
  --i-table $PROJECT/Merged-runs/filtered-dropped-metadata-table.qza \
  --o-visualization $PROJECT/Merged-runs/filtered-dropped-metadata-table.qzv

# Weighted taxanomic classification
  # This is a machine learning algorithm that will be trained on the Greengenes database.
  # Can be downloaded from: https://data.qiime2.org/2023.9/common/gg-13-8-99-nb-weighted-classifier.qza
  # The classfier used has been 'trained' based off of 14 habitat types includin human-oral

# Can either use the human oral classifier: https://zenodo.org/records/6395539
# Or build your own by using the instructions on Classifier.sh. 

qiime feature-classifier classify-sklearn \
  --i-reads $PROJECT/Merged-runs/filtered-metadata-table-sequences.qza \
  --i-classifier $WORK/database/99_classifier.qza \
  --o-classification $PROJECT/Merged-runs/taxonomy-merged-final.qza
  
qiime metadata tabulate \
  --m-input-file $PROJECT/Merged-runs/taxonomy-merged-final.qza \
  --o-visualization $PROJECT/Merged-runs/taxonomy-merged-final.qzv

qiime taxa barplot \
  --i-table $PROJECT/Merged-runs/filtered-metadata-table.qza \
  --i-taxonomy $PROJECT/Merged-runs/taxonomy-merged-final.qza \
  --m-metadata-file $PROJECT/Mapping/Manifest-map.txt \
  --o-visualization $PROJECT/Merged-runs/taxa-bar-plots-merged-final.qzv

# Phylogenetic tree generation
qiime alignment mafft \
  --i-sequences $PROJECT/Merged-runs/merged_sequences.qza \
  --o-alignment $PROJECT/Merged-runs/merged-final-aligned-seqs.qza

qiime alignment mask \
  --i-alignment $PROJECT/Merged-runs/merged-final-aligned-seqs.qza \
  --o-masked-alignment $PROJECT/Merged-runs/merged-aligned-masked-seqs.qza

qiime phylogeny fasttree \
  --i-alignment $PROJECT/Merged-runs/merged-aligned-masked-seqs.qza \
  --o-tree $PROJECT/Merged-runs/merged-unrooted-tree.qza

qiime phylogeny midpoint-root \
  --i-tree $PROJECT/Merged-runs/merged-unrooted-tree.qza \
  --o-rooted-tree $PROJECT/Merged-runs/merged-midrooted-tree.qza

# Export files
mkdir $PROJECT/Merged-runs/exported-files

qiime tools export \
  --input-path $PROJECT/Merged-runs/merged-unrooted-tree.qza \
  --output-path $PROJECT/Merged-runs/exported-files/merged-unrooted-tree

qiime tools export \
  --input-path $PROJECT/Merged-runs/taxonomy-merged-final.qza \
  --output-path $PROJECT/Merged-runs/exported-files/taxonomy-merged-final

qiime tools export \
  --input-path $PROJECT/Merged-runs/filtered-metadata-table.qza \
  --output-path $PROJECT/Merged-runs/exported-files/merged_table

