# README
### Purpose of this script: to automate the demultiplexing, denoising, removing chimeras etc. 

## Notes on using HPC
    qstat -u <username>
        ### to get a list of currently running jobs associated by your username
    qstat <job id>
        ### specific job info
    qstat <job id> -f
        ### Get full information on the job, including resources used
Always use absolute directories, just to make it fail safe. 

## Notes on fastq files
Three fastq files for each run. R1 is forward, R2 is barcode and R3 is reverse. Where there is no R3: R1 forward, I1 barcodes and R2 is reverse.
    There are maps for MSQ92 and MSQ91 but unknown with which fastq files they belong to...
Folder 16S-fastq holds all fastq files of all runs. Sorted by MSQ run.
    Copied all fastq files into one folder: fastq which is in the parent MSQ-run folder
    Structure: MSQ-run/fastq/<fastq files>
    Then compressed them
```
for f in "<ALL MSQ runs>"; do gzip -9 "$f"/fastq/*.fastq; done
```
### Renamed all files to forward.fastq.gz, reverse.fastq.gz and barcodes.fastq.gz
Qiime2 won't recognise <MSQ-run>-forward.fastq.gz 
```
for f in "<All MSQ runs>"; do cd "$f"/fastq/; mv lane1_NoIndex_L001_R1_001.fastq forward.fastq.gz; mv lane1_NoIndex_L001_R2_001.fastq barcodes.fastq.gz; mv lane1_NoIndex_L001_R3_001.fastq reverse.fastq.gz; done
```

### Script folder
Holds the scripts submitted to the cluster as for loops. 

If this error shows
    -bash: /var/spool/PBS/mom_priv/jobs/8392841.pbs.SC: /bin/sh^M: bad interpreter: No such file or directory
You need to convert it to UNIX format by doing, depends on what OS the script was written.
    dos2unix <bash script>

## Qiime2-16S-pipeline.sh
The entire pipeline from demultiplexing to denoising. Stops just before merging. 
    You need to specify the MSQ run for the RUN variable e.g.,
    ```
    RUN = "RUN1 RUN2 RUN3"
    ```
    And specify the folder of your project for the PROJECT variable e.g.,
    ```
    PROJECT="$WORK/Project1-healthy-controls"
    ```
The pipeline will take the fastq files in 16S-fastq and use this as the input file for everything. 
After this it'll work locally from your specified project folder, this is done to ensure you do not accidentally alter the original files. 

## Merge-filter-classify-tables.sh
This merges all the qiime2 artefacts, tables, feature data, sequences etc. 
You can then filter out samples based on your manifest file i.e., only healthy patients.

To create a manifest file, use Excel to just copy paste, at least, the #SampleID into an excel file. 
There's a R file which will merge the barcodes with the sample-ids you specified in the excel file. 
At the end of it the R script will create a file called: Manifest-map.txt. 

Create a folder on your HPC: 
```
mkdir $PROJECT/Mapping
```
So that the Manifest-map.txt can be found at $PROJECT/Mapping/Manifest-map.txt 

## Classifier.sh
This script creates the weighted taxonomy classifier as per Qiime2 tutorials. 
The inputs for this script can be found: https://data.qiime2.org/2023.9/common/gg-13-8-99-nb-weighted-classifier.qza
