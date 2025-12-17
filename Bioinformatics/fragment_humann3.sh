#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:mem=80gb
#PBS -N humann3_batch

set -euo pipefail 
IFS=$'\n\t' 

# This script is for PacBio reads. 
# It includes a part where you're fragmenting reads into Illumina sized chunks. 

# --- Configuration ---
WORKDIR="${PBS_O_WORKDIR:-/rds/general/user/mteng/ephemeral/Project11-Nanopore-test/humann3}"
INPUT_DIR="../PacBio-amplified-reads/02_host_removed"
OUTPUT_DIR="humann3"
FASTQ_LIST="../PacBio-amplified-reads/basename.txt"
# Number of threads to give HUMAnN3 (independent of PBS ncpus)
HUMANN3_THREADS=10
# Fragment long reads into smaller pseudo-reads so Bowtie2 (used by HUMAnN3)
# can align them and produce stratified outputs. Adjust as needed.
FRAG_LEN=200      # fragment length in bp
MIN_FRAG=50       # minimum fragment length to retain

# --- Setup ---
cd "$WORKDIR" || { echo "FATAL ERROR: Failed to cd to $WORKDIR" >&2; exit 1; }

# Activate Conda environment
eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate biobakery_env

mkdir -p "$OUTPUT_DIR"

# --- PBS ARRAY JOB CHECK ---
if [ -z "${PBS_ARRAY_INDEX:-}" ]; then
    echo "ERROR: PBS_ARRAY_INDEX is not set. Run this script using qsub -t 1-N." >&2
    exit 1
fi

# Get the base name (e.g., PFND111) from the list
BASE_NAME=$(awk -v idx="$PBS_ARRAY_INDEX" 'NR==idx{print; exit}' "$FASTQ_LIST")

if [ -z "$BASE_NAME" ]; then
    echo "ERROR: No input file entry found for index $PBS_ARRAY_INDEX in $FASTQ_LIST" >&2
    exit 1
fi

# Construct the full input filename from the base name
INPUT_FILENAME="${BASE_NAME}_nonhost.fastq.gz"

# Construct the full input path
INPUT_PATH="${INPUT_DIR}/${INPUT_FILENAME}"

if [ ! -f "$INPUT_PATH" ]; then
    echo "ERROR: Input file does not exist at: $INPUT_PATH" >&2
    exit 1
fi

# Set threads variable and prepare sample output
echo "--- Starting job for: ${BASE_NAME} (Index: $PBS_ARRAY_INDEX) ---"
echo "Reading from: $INPUT_PATH"
echo "Using $HUMANN3_THREADS threads for HUMAnN3"

# Create sample-specific output directory
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${BASE_NAME}_humann3"
mkdir -p "$SAMPLE_OUTPUT_DIR" || { echo "ERROR: Failed to create output directory: $SAMPLE_OUTPUT_DIR" >&2; exit 1; }

# Decompress FASTQ.gz into sample output directory
RAW_FASTQ="${SAMPLE_OUTPUT_DIR}/${BASE_NAME}_nonhost.fastq"
echo "Decompressing $INPUT_PATH -> $RAW_FASTQ"
gzip -dc "$INPUT_PATH" > "$RAW_FASTQ" || { echo "ERROR: Failed to decompress $INPUT_PATH" >&2; rm -f "$RAW_FASTQ"; exit 1; }

# Fragment long reads into smaller pseudo-reads so Bowtie2 can align them
FRAG_FASTQ="${SAMPLE_OUTPUT_DIR}/${BASE_NAME}_frags.fastq"
echo "Fragmenting reads to ${FRAG_LEN} bp (min ${MIN_FRAG}) -> $FRAG_FASTQ"
export IN_FASTQ="$RAW_FASTQ"
export OUT_FASTQ="$FRAG_FASTQ"
export FRAG_LEN
export MIN_FRAG
python3 - <<'PY'
import os
infile = os.environ['IN_FASTQ']
outfile = os.environ['OUT_FASTQ']
frag = int(os.environ.get('FRAG_LEN', '200'))
minlen = int(os.environ.get('MIN_FRAG', '50'))
with open(infile) as fin, open(outfile, 'w') as fout:
    while True:
        header = fin.readline()
        if not header:
            break
        seq = fin.readline().rstrip('\n')
        plus = fin.readline()
        qual = fin.readline().rstrip('\n')
        if not qual:
            break
        seqlen = len(seq)
        i = 0
        part = 0
        while i < seqlen:
            sub = seq[i:i+frag]
            subq = qual[i:i+len(sub)]
            if len(sub) >= minlen:
                h = header.strip()
                if h.startswith('@'):
                    h = h[1:]
                fout.write('@' + h + '_frag' + str(part) + '\n')
                fout.write(sub + '\n')
                fout.write('+\n')
                fout.write(subq + '\n')
                part += 1
            i += frag
PY

if [ ! -s "$FRAG_FASTQ" ]; then
    echo "ERROR: Fragmented FASTQ is empty: $FRAG_FASTQ" >&2
    exit 1
fi

# --- Run Pipeline using fragmented reads ---
echo "Running HUMAnN3 on fragmented reads..."
if humann --input "$FRAG_FASTQ" --output "$SAMPLE_OUTPUT_DIR" --threads "$HUMANN3_THREADS" --resume; then
    echo "SUCCESS: HUMAnN3 analysis completed for ${BASE_NAME}."
    echo "Output saved to: $SAMPLE_OUTPUT_DIR"
    # Optionally remove intermediate decompressed file to save space
    rm -f "$RAW_FASTQ"
else
    echo "CRITICAL ERROR: Pipeline failed for ${BASE_NAME}. Check error logs in: $SAMPLE_OUTPUT_DIR" >&2
    exit 1
fi

