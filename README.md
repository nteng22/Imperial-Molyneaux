# Imperial-Molyneaux
 Projects from the Molyneaux Imperial lab

## Qiime2
This contains all the scripts needed to run qiime2 on the HPC. <br>
[Qiime2-16S-pipeline.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Qiime2-16S-pipeline.sh), runs the reads from the sequencing company/machine. Demultiplexes, denoises and clustering. 
The tutorial with more information can be found [here](https://docs.qiime2.org/2023.9/plugins/available/diversity/alpha/).

All MSQ runs to date (November 2023) have gone through this process and is found on the Fibrosis HPC drive under "Processed-runs".
For any new MSQ run, you need to run this pipeline.
**Double check what type of barcodes were used, some barcodes are golay compatible. 
If this is the case you need to specify --p-rev-comp-barcodes and --p-rev-comp-mapping-barcodes. In the older runs only one needs to be included.**
 *The only runs that are not golay barcoded are MSQ113 and MSQ114*
 
You can see how many reads are expected by going to the fastq file (forward or reverse) and then:
```
wc -l <fastq file>
```
>Divide the number you get by 4, and you get a rough idea of how many reads are expected.

To run it on the HPC you need to create your project folder and change the path in the script, **PROJECT**. Also specify the MSQ run/number, **RUN**

[Merge-filter-classify-tables.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Merge-filter-classify-tables.sh) is to assign taxonomy and select for samples that you want. 
To do this you can create an excel spreadsheet with the samples you want taken from the "16S-Sequencing-Summary-Runs.xlsx" found on the Fibrosis drive or in this repository. **This has not been QCed properly**
*e.g. if you only want Fibrosis patients you use the master sheet to copy and paste it into a new spreadsheet* 
**The most imporant thing is the sampleID/#SampleID/index column**

Then use the Manifest-map.R script to merge the sampleIDs you want with the manifest maps submitted to the sequencing team. The R script will save it with a generic name and in the correct format. 
Copy this to your project folder on the HPC. 
To run this script on the HPC, the variables **PROJECT, MSQ** need to be defined accordingly. 

[Classifier.sh](https://github.com/nteng22/Imperial-Molyneaux/blob/main/Qiime2/Classifier.sh) is the script to 'train' for weighted taxonomy assignment. 

## COVID-BAL-OralRinse
Data from COVID patients and non-COVID patients (healthy), where healthy is defined as patients who did not have a positive COVID test at the time of sampling. <br>
<br>
All files were provided by sequencing company and can be used on the galaxy webpage of Qiime2 to visualise initially or download files e.g. <br>

## Healthy-controls
Sequencing data from historical runs defined as diagnosed with non-fibortic disease. 

## IPF-progressive
Project to look into the microbiota of stable vs. progressive IPF patients. 
