#!/bin/sh

## ask PBS for time (format hh:mm:ss)
#PBS -l walltime=24:00:00
## memory (per node)
#PBS -lselect=1:ncpus=20:mem=20gb

##load application module
module load anaconda3/personal
source activate qiime2-amplicon-2024.2
# THINGS TO CHANGE IN THE SCRIPTS #
  # List the runs you want to include in the analysis in the RUN variable.
      # The format will be "RUN1 RUN2 RUN3"
RUN="NorthUmbria_stool-baseline"
PROJECT="$WORK/Project9-PROFOUND/stool/ms_baseline-samples"
MSQ="$WORK/Processed-runs"

cd $WORK

# Things you need to have before running this. This needs to be created in the HPC environment. 
  # Specify the path to your project folder:

# Change the directories for the reference databases if the structure is different:
  # For 99_otu.fasta, just before quality control
  # For gg-13-8-99-515-806-nb-weighted-classifier.qza, just after weighted taxonomy classification
for f in $RUN
do

mkdir $MSQ/"$f"
mkdir $MSQ/"$f"/demux

# Demultiplex samples, using demultiplexed samples follows a slightly different protocol
  # More info: https://docs.qiime2.org/2024.2/tutorials/read-joining/ 
  # More infor: https://docs.qiime2.org/2024.2/tutorials/atacama-soils/
  # Import the reads into qiime2
# qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path $PROJECT/../fastq/PFND-stool \
#   --output-path $MSQ/"$f"/demux/"$f"_demux.qza \
#   --input-format CasavaOneEightSingleLanePerSampleDirFmt

# qiime demux filter-samples \
#     --i-demux $MSQ/"$f"/demux/"$f"_demux.qza \
#     --m-metadata-file $PROJECT/Mapping/Manifest-Map_stool.txt \
#     --o-filtered-demux $MSQ/"$f"/demux/"$f"_demux-filtered.qza

# # # Obtain a summary of the demultiplexed data
# qiime demux summarize \
#    --i-data $MSQ/"$f"/demux/"$f"_demux-filtered.qza  \
#    --o-visualization $MSQ/"$f"/"$f"_demux-filtered.qzv

# DENOISE SAMPLES i.e., remove chimeras and artefacts from sequencing process. 
mkdir $MSQ/"$f"/denoise 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $MSQ/NorthUmbria_stool/demux/NorthUmbria_stool_demux-filtered.qza \
  --p-trim-left-f 15 \
  --p-trim-left-r 15 \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 230 \
  --o-representative-sequences $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --o-table $MSQ/"$f"/denoise/"$f"_table-dada2.qza \
  --o-denoising-stats $MSQ/"$f"/denoise/"$f"_dada-2_stats.qza

# Merge duplicated samples
qiime feature-table merge-seqs \
 --i-data $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
 --o-merged-data $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza

qiime metadata tabulate \
  --m-input-file $MSQ/"$f"/denoise/"$f"_dada-2_stats.qza \
  --o-visualization $MSQ/"$f"/denoise/"$f"_denoise-stats.qzv
    
  # QUALITY CONTROL
mkdir $MSQ/"$f"/quality-control
    
    # Matches the query sequences to a reference database of sequences, in this case the full-length.fna file.
    # Saves the files into two files, sequence-hits.qza and sequence-misses.qza, where sequence-misses.qza contains the sequences that did not match the reference database.
      # This can be used to remove sequences that are not 16S rRNA gene sequences, or not at a 99% match to the reference database.
      # Reference database is based off of GreenGenes, which is a curated, cleaned-up database based off of SILVA.

qiime quality-control exclude-seqs \
  --i-query-sequences $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --i-reference-sequences $WORK/database/2022.10.backbone.full-length.fna.qza \
  --p-method vsearch \
  --p-perc-identity 0.97 \
  --p-perc-query-aligned 0.97 \
  --p-threads 4 \
  --o-sequence-hits $MSQ/"$f"/quality-control/"$f"_sequence-hits.qza \
  --o-sequence-misses $MSQ/"$f"/quality-control/"$f"_sequence-misses.qza

# Filtering sequences to remove misses from the table, specified by the sequence-misses.qza file and --p-exclude-ids.
qiime feature-table filter-features \
  --i-table $MSQ/"$f"/denoise/"$f"_table-dada2.qza \
  --m-metadata-file $PROJECT/Manifest-baseline.txt \
  --o-filtered-table $MSQ/"$f"/quality-control/"$f"_filtered-table.qza \
  --p-exclude-ids

qiime feature-table summarize \
  --i-table $MSQ/"$f"/quality-control/"$f"_filtered-table.qza \
  --o-visualization $MSQ/"$f"/quality-control/"$f"_filtered-table.qzv \
  --m-sample-metadata-file $PROJECT/Manifest-baseline.txt

qiime feature-table tabulate-seqs \
  --i-data $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --o-visualization $MSQ/"$f"/quality-control/"$f"_rep-seqs.qzv

cd $WORK

mkdir $PROJECT/Merged-runs-nonv4-16S

# Create a feature table artifact based on the filtered metadata table.
qiime tools export \
  --input-path $MSQ/"$f"/quality-control/"$f"_sequence-hits.qza \
  --output-path $PROJECT/Merged-runs-nonv4-16S/filtered-metadata-table-sequences

# Using the non-V4 classifier method
qiime greengenes2 non-v4-16s \
  --i-table $MSQ/"$f"/quality-control/"$f"_filtered-table.qza \
  --i-sequences $MSQ/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --i-backbone $WORK/database/2022.10.backbone.full-length.fna.qza \
  --o-mapped-table $PROJECT/Merged-runs-nonv4-16S/filtered-table-gg2.qza \
  --o-representatives $PROJECT/Merged-runs-nonv4-16S/filtered-rep-seqs-gg2.qza

# taxonomy from table outputs  
qiime greengenes2 taxonomy-from-table \
    --i-reference-taxonomy $WORK/database/2022.10.taxonomy.asv.nwk.qza \
    --i-table $PROJECT/Merged-runs-nonv4-16S/filtered-table-gg2.qza \
    --o-classification $PROJECT/Merged-runs-nonv4-16S/filtered-taxonomy-gg2.qza

qiime metadata tabulate \
  --m-input-file $PROJECT/Merged-runs-nonv4-16S/filtered-taxonomy-gg2.qza \
  --o-visualization $PROJECT/Merged-runs-nonv4-16S/filtered-taxonomy-gg2.qzv

qiime taxa barplot \
  --i-table $PROJECT/Merged-runs-nonv4-16S/filtered-table-gg2.qza \
  --i-taxonomy $PROJECT/Merged-runs-nonv4-16S/filtered-taxonomy-gg2.qza \
  --m-metadata-file $PROJECT/Manifest-baseline.txt \
  --o-visualization $PROJECT/Merged-runs-nonv4-16S/FeatureTable_barplot-features.qzv

# Phylogenetic tree generation
qiime alignment mafft \
  --i-sequences $PROJECT/Merged-runs-nonv4-16S/filtered-rep-seqs-gg2.qza \
  --o-alignment $PROJECT/Merged-runs-nonv4-16S/final-aligned-seqs-gg2.qza

qiime alignment mask \
  --i-alignment $PROJECT/Merged-runs-nonv4-16S/final-aligned-seqs-gg2.qza \
  --o-masked-alignment $PROJECT/Merged-runs-nonv4-16S/aligned-masked-seqs.qza

qiime phylogeny fasttree \
  --i-alignment $PROJECT/Merged-runs-nonv4-16S/aligned-masked-seqs.qza \
  --o-tree $PROJECT/Merged-runs-nonv4-16S/unrooted-tree.qza

qiime phylogeny midpoint-root \
  --i-tree $PROJECT/Merged-runs-nonv4-16S/unrooted-tree.qza \
  --o-rooted-tree $PROJECT/Merged-runs-nonv4-16S/midrooted-tree.qza

done