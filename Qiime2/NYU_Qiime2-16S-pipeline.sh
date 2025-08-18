#!/bin/sh

## ask PBS for time (format hh:mm:ss)
#PBS -l walltime=47:55:00
## memory (per node)
#PBS -lselect=1:ncpus=32:mem=32gb

##load application module
module load anaconda3/personal
source activate qiime2-amplicon-2023.9
# THINGS TO CHANGE IN THE SCRIPTS #
  # List the runs you want to include in the analysis in the RUN variable.
      # The format will be "RUN1 RUN2 RUN3"
RUN="MSQ82 MSQ92 MSQ113 MSQ114"

cd $WORK

# Things you need to have before running this. This needs to be created in the HPC environment. 
  # Specify the path to your project folder:
PROJECT="$WORK/Processed-runs/with-comp-barcodes"

# Change the directories for the reference databases if the structure is different:
  # For 99_otu.fasta, just before quality control
  # For gg-13-8-99-515-806-nb-weighted-classifier.qza, just after weighted taxonomy classification
for f in $RUN
do

mkdir $PROJECT/"$f"
mkdir $PROJECT/"$f"/demux

# Demultiplex samples, i.e., assign sequences to samples based on their barcodes. 
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path $WORK/16S-fastq/"$f"/fastq  \
  --output-path $PROJECT/"$f"/"$f"_emp-paired-end-sequences.qza

qiime demux emp-paired \
    --i-seqs $PROJECT/"$f"/"$f"_emp-paired-end-sequences.qza \
    --m-barcodes-file $WORK/16S-fastq/MSQ-maps/"$f"-Map.txt \
    --m-barcodes-column BarcodeSequence \
    --o-per-sample-sequences $PROJECT/"$f"/"$f"_demux.qza \
    --p-rev-comp-mapping-barcodes \
    --p-rev-comp-barcodes \
    --o-error-correction-details $PROJECT/"$f"/"$f"_demux-details.qza
    
    
qiime demux summarize \
   --i-data $PROJECT/"$f"/"$f"_demux.qza  \
   --o-visualization $PROJECT/"$f"/demux/"$f"_demux.qzv

  # DENOISE SAMPLES i.e., remove chimeras and artefacts from sequencing process. 
mkdir $PROJECT/"$f"/denoise 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $PROJECT/"$f"/"$f"_demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 140 \
  --p-trunc-len-r 140 \
  --o-representative-sequences $PROJECT/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --o-table $PROJECT/"$f"/denoise/"$f"_table-dada2.qza \
  --o-denoising-stats $PROJECT/"$f"/denoise/"$f"_dada-2_stats.qza

qiime metadata tabulate \
    --m-input-file $PROJECT/"$f"/denoise/"$f"_dada-2_stats.qza \
    --o-visualization $PROJECT/"$f"/denoise/"$f"_denoise-stats.qzv
    
  # 99_otu.fasta is the reference database for the 16S rRNA gene sequences, http://ftp.microbio.me/greengenes_release/gg_13_8_otus/rep_set/ 
  # Download the database into a folder called database in your work directory.
    # Alternatively change the path to the database in the command below.
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $WORK/database/99_otus.fasta \
  --output-path $WORK/database/99_otus.qza
    
  # QUALITY CONTROL
mkdir $PROJECT/"$f"/quality-control
    
    # Matches the query sequences to a reference database of sequences, in this case the 99_otus.qza file.
    # Saves the files into two files, sequence-hits.qza and sequence-misses.qza, where sequence-misses.qza contains the sequences that did not match the reference database.
      # This can be used to remove sequences that are not 16S rRNA gene sequences, or not at a 99% match to the reference database.
      # Reference database is based off of GreenGenes, which is a curated, cleaned-up database based off of SILVA.
qiime quality-control exclude-seqs \
  --i-query-sequences $PROJECT/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --i-reference-sequences $WORK/database/99_otus.qza \
  --p-method vsearch \
  --p-perc-identity 0.99 \
  --p-perc-query-aligned 0.99 \
  --p-threads 4 \
  --o-sequence-hits $PROJECT/"$f"/quality-control/"$f"_sequence-hits.qza \
  --o-sequence-misses $PROJECT/"$f"/quality-control/"$f"_sequence-misses.qza

# Filtering sequences to remove misses from the table, specified by the sequence-misses.qza file and --p-exclude-ids.
qiime feature-table filter-features \
  --i-table $PROJECT/"$f"/denoise/"$f"_table-dada2.qza \
  --m-metadata-file $PROJECT/"$f"/quality-control/"$f"_sequence-misses.qza \
  --o-filtered-table $PROJECT/"$f"/quality-control/"$f"_filtered-table.qza \
  --p-exclude-ids

qiime feature-table summarize \
  --i-table $PROJECT/"$f"/quality-control/"$f"_filtered-table.qza \
  --o-visualization $PROJECT/"$f"/quality-control/"$f"_filtered-table.qzv \
  --m-sample-metadata-file $WORK/16S-fastq/MSQ-maps/"$f"-Map.txt

qiime feature-table tabulate-seqs \
  --i-data $PROJECT/"$f"/denoise/"$f"_rep-seqs-dada2.qza \
  --o-visualization $PROJECT/"$f"/quality-control/"$f"_rep-seqs.qzv
done 
