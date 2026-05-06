# Imperial-Molyneaux

> Microbiome bioinformatics pipelines from the Molyneaux Imperial Lab. Scripts for amplicon sequencing (16S), metagenomics, and functional profiling.

---

## Table of Contents

- [Qiime2 Amplicon Sequencing](#qiime2-amplicon-sequencing)
  - [PE150bp Chemistry (Archived)](#pe150bp-chemistry-archived)
  - [PE250bp Chemistry](#pe250bp-chemistry)
  - [Manifest Creation](#manifest-creation)
- [Functional Profiling with PICRUSt2](#functional-profiling-with-picrust2)
- [Published Projects](#published-projects)
  - [POSTCODE](#postcode)
  - [IPF-SCFA](#ipf-scfa)
- [Bioinformatics — PacBio Metagenomics](#bioinformatics--pacbio-metagenomics)

---

## Qiime2 Amplicon Sequencing

All scripts for running [QIIME2](https://amplicon-docs.qiime2.org/en/latest/) on the HPC. Covers demultiplexing, denoising, taxonomy assignment, and filtering.

---

### PE150bp Chemistry *(Archived)*

> These scripts are archived for reproducibility. All MSQ runs to date (November 2023) have been processed with this pipeline and are stored on the Fibrosis HPC drive under `Processed-runs`.

**[Qiime2-16S-pipeline.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Qiime2-16S-pipeline.sh)**
Runs raw reads through demultiplexing, denoising, and clustering.

> ⚠️ **Barcode check:** Some barcodes are Golay-compatible — if so, you must specify `--p-rev-comp-barcodes` and `--p-rev-comp-mapping-barcodes`. In older runs, only one flag is needed.
> MSQ113 and MSQ114 are the only known runs that are **not** Golay-barcoded.

**Estimating read counts:**
```bash
wc -l <fastq_file>
# Divide result by 4 → approximate read count
# High biomass samples (saliva, stool): expect ~8,000–10,000 reads
```

**To run on HPC:**
1. Create your project folder.
2. Set the `PROJECT` variable in the script to the correct path.
3. Set the `RUN` variable to your MSQ run number.

---

**[Merge-filter-classify-tables.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Merge-filter-classify-tables.sh)**
Assigns taxonomy and filters to your samples of interest.

Steps:
1. Create a sample list from `16S-Sequencing-Summary-Runs.xlsx` (on the Fibrosis drive or in this repo).
2. Run **[Manifest-map.R](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Manifest-map.R)** to merge your sample IDs with the manifest maps submitted to the sequencing team. The script saves the output in the correct format with a generic name.
3. Copy the output to your HPC project folder.
4. Set `PROJECT` and `MSQ` variables before running.

> ⚠️ **Note:** The most critical column is `sampleID` / `#SampleID` / `index`. This script has not been fully QC'd.

---

**[Classifier.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Classifier.sh)**
Trains a classifier for weighted taxonomy assignment using a downloaded reference database and taxonomy file.

> 💡 **Recommendation:** Train your own classifier rather than using a pre-trained one — pre-trained classifiers may produce incorrect results. If training isn't feasible, use a reference database alignment instead. Adapt by downloading your preferred reference database.

---

### PE250bp Chemistry

**[Qiime2-PairedEnd_stool_nonV4.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Qiime2-PairedEnd_stool_nonV4.sh)**
Full pipeline for non-150bp reads, from raw reads to QIIME2 artefacts.

Key differences from the 150bp pipeline:
- Uses **Greengenes2** instead of Greengenes 1.
- Uses **feature-to-table** instead of Naïve Bayes classification — this preserves information from the full 250bp amplicons.
- **Does not use** `qiime feature-classifier classify-sklearn`, as that classifier is trained on 150bp V4 reads and would discard useful 250bp information.

**Required database:** [Greengenes2 (2022.10)](https://ftp.microbio.me/greengenes_release/2022.10/)

---

### Manifest Creation

**[Create-Manifest_PE.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Create-Manifest_PE.sh)**
Creates the manifest map needed to import demultiplexed sequences. Requires absolute paths to FASTQ files. Update the folder name accordingly.

**[NuOmics_Manifest-creation.R](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/NuOmics_Manifest-creation.R)**
Creates manifests primarily for the PROFOUND study (sent to NuOmics). Incorporates Picogreen values from the sequencing team.

> ⚠️ **Note:** There are known duplicate samples sent across different plates, which were not flagged by the sequencing team.

---

## Functional Profiling with PICRUSt2

The QIIME2 and PICRUSt2 environments have conflicting Python dependencies, so the pipeline is split into two scripts:

| Script | Purpose | Runtime |
|---|---|---|
| [Picrust2-prep.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Picrust2-prep.sh) | Prepares the biom feature table and ASV FASTA sequences | < 15 min |
| [Picrust2-run-pipeline.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Picrust2-run-pipeline.sh) | Runs the full PICRUSt2 pipeline | Variable |

---

## Published Projects

### POSTCODE

**Lung microbiota in post-COVID patients with long-lasting lung abnormalities.**

- **Analysis:** [POSTCODE-manuscript.R](https://github.com/molyneaux-lab/POSTCODE)
- **Repository includes:** SRA data, sample name files, metadata, and QIIME2 artefacts for full reproducibility.
- Healthy controls = non-fibrotic patients.
- CHP and IPF samples from a prior published project: [PRJNA609242](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA609242).

---

### IPF-SCFA

**Short-chain fatty acids (SCFA) in the lungs of IPF patients.**

- **Analysis:** [IPF-SCFA.R](https://github.com/molyneaux-lab/IPF-SCFA)
- **Data:** [PRJNA772278](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA772278) — SRA data and metadata sheets.

---

## Bioinformatics — PacBio Metagenomics

Scripts for long-read PacBio HiFi metagenomics on the HPC.

### Directory Structure

```text
/rds/general/user/mteng/ephemeral/pb-metagenomics-tools/
│
├── basename.txt                 # Sample IDs (e.g., PatientA) — drives all looping/naming
│
├── 01_fastq/                    # Raw PacBio HiFi sequences (input)
│   ├── sample1.fastq.gz
│   └── sample2.fastq.gz
│
├── 02_host_removed/             # Cleaned reads (human reads removed)
│   ├── sample1_clean.fastq.gz
│   └── sample2_clean.fastq.gz
│
└── 03_metaflye_assemblies/      # De novo metagenome assemblies
    ├── sample1/
    │   ├── assembly.fasta
    │   ├── assembly_info.txt
    │   └── flye.log
    └── sample2/
        ├── assembly.fasta
        └── assembly_info.txt
```

> **`basename.txt`** — contains one sample ID per line. All scripts loop over this file to name outputs consistently. Simplify FASTQ filenames in `01_fastq/` before running any pipeline step.

---

### Taxonomic Profiling Workflow

1. **Remove host reads** — [Minimap2 host removal](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Bioinformatics/Minimap2%20host%20removal)
2. **Assemble reads** — [metaflye](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Bioinformatics/metaflye-pacbio.sh)
3. **Taxonomic classification** — [PacBio Minimap2-MEGAN tutorial](https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/Tutorial-Taxonomic-Profiling-Minimap-Megan.md)

### Functional Profiling Workflow

1. **Run HUMAnN3** — [Humann3-fragment.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Bioinformatics/Humann3-fragment.sh) *(includes read fragmentation step)*
