# Imperial-Molyneaux

Projects from the Molyneaux Imperial lab

## Qiime2

This contains all the scripts needed to run qiime2 on the HPC. <br> [Qiime2-16S-pipeline.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Qiime2-16S-pipeline.sh), runs the reads from the sequencing company/machine. Demultiplexes, denoises and clustering. The tutorial with more information can be found [here](https://docs.qiime2.org/2023.9/plugins/available/diversity/alpha/).

All MSQ runs to date (November 2023) have gone through this process and is found on the Fibrosis HPC drive under "Processed-runs". For any new MSQ run, you need to run this pipeline. **Double check what type of barcodes were used, some barcodes are golay compatible. If this is the case you need to specify --p-rev-comp-barcodes and --p-rev-comp-mapping-barcodes. In the older runs only one needs to be included.** *The only runs that are not golay barcoded are MSQ113 and MSQ114 as far as I'm currently aware*

You can see how many reads are expected by going to the fastq file (forward or reverse) and then:

```         
wc -l <fastq file>
```

> Divide the number you get by 4, and you get a rough idea of how many reads are expected. For high biomass samples (saliva, stool) you'd expect 8000-10000 reads.

To run it on the HPC you need to create your project folder and change the path in the script for the variable: **PROJECT**. If you have a new MSQ run you need to specify the MSQ run/number in the variable: **RUN**. The instructions are also in the Qiime2 script.

[Merge-filter-classify-tables.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Merge-filter-classify-tables.sh) is to assign taxonomy and select for samples that you want. To do this you can create an excel spreadsheet with the samples you want taken from the "16S-Sequencing-Summary-Runs.xlsx" found on the Fibrosis drive or in this repository. **This has not been QCed properly** *e.g. if you only want Fibrosis patients you use the master sheet to copy and paste it into a new spreadsheet* **The most important thing is the sampleID/#SampleID/index column**

Then use the Manifest-map.R script to merge the sampleIDs you want with the manifest maps submitted to the sequencing team. The R script will save it with a generic name and in the correct format. Copy this to your project folder on the HPC. To run this script on the HPC, the variables **PROJECT, MSQ** need to be defined accordingly.

[Classifier.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Classifier.sh) is the script to 'train' for weighted taxonomy assignment.

[Create-Manifest_PE.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Create-Manifest_PE.sh). Will create manifest map needed for importing demultiplexed sequences. You need to provide the absolute pathways to the FASTQ files. Change the name of the folder accordingly. 

[Qiime2-PairedEnd_stool_nonV4.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Qiime2-PairedEnd_stool_nonV4.sh) is the script to run anything else than 150bp reads. The original pipeline from NYU utilises greengenes 1 and is optimised for V4 region and up to 150bp PE reads. Anything else it's recommended to use greengenes2 and use feature-to-table as opposed to NB classification. This runs it from start to finish ending with the qiime2 artefacts and a barplot. 


## POSTCODE
Lung microbiota of post-COVID patients with long lasting lung abnormalities. [POSTCODE-manuscript.R] (<https://github.com/molyneaux-lab/POSTCODE>) contains the data analysis and plots made. The repository includes: the SRA data, sample_names_files, metadata files and Qiime2 artefacts to run the analysis yourself. Healthy controls include non-fibrotic patients. CHP and IPF samples were taken from a previous published project from our lab found at [PRJNA609242](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA609242).
