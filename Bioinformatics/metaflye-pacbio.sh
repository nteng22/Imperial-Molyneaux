#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=64gb

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate metaflye

cd /rds/general/user/mteng/ephemeral/Project11-Nanopore-test

flye --pacbio-hifi PFND164.fastq.gz --meta --out-dir PFND164-metaflye --resume
