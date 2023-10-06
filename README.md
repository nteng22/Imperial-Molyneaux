# Imperial-Molyneaux
 Projects from the Molyneaux Imperial lab

## COVID-BAL-OralRinse
Data from COVID patients and non-COVID patients (healthy), where healthy is defined as patients who did not have a positive COVID test at the time of sampling. <br>
<br>
All files were provided by sequencing company and can be used on the galaxy webpage of Qiime2 to visualise initially or download files e.g. <br>

> Used **taxa-bar-plots_norev_quality.qzv** as input file and exported as .csv file. Can use any taxa level, only needed index column as sampleID. <br>
> [Qiime2-Genus.csv](https://github.com/nteng22/Imperial-Molyneaux/blob/main/COVID-BAL-OralRinse/Quality-check/Qiime2-Genus.csv), Visualised using [galaxy](https://view.qiime2.org/). <br>
<br>

> General analysis took place, alpha and beta diversity were shannon and NMDS plots respectively.
> NMDS showed significance which was further investigated using taxa specific analysis to determine what taxa best explains the NMDS plot variance. 

### [Quality-check](https://github.com/nteng22/Imperial-Molyneaux/tree/main/COVID-BAL-OralRinse/Quality-check)
Contains all the original files to do the QC process. Most of the code is based off of Rachele's work. 
> [Metadata.csv](https://github.com/nteng22/Imperial-Molyneaux/blob/main/COVID-BAL-OralRinse/Quality-check/metadata.csv) is without Bavi's samples (pre-fixed with PH) and without mock (pre-fixed with MOC) flagged as *Sample-type == EXCLUDE* <br>
> [Metadata.R](https://github.com/nteng22/Imperial-Molyneaux/blob/main/COVID-BAL-OralRinse/Quality-check/Metadata.R), separates the sampleID submitted for sequencing to the disease types, sample types and patientID. This is done first, before normalisation and removing contaminants. <br>
> [Contamination-Normalisation.R](https://github.com/nteng22/Imperial-Molyneaux/blob/main/COVID-BAL-OralRinse/Quality-check/Contamination-Normalisation.R), script to identify and remove contaminants from samples. This was based off of 'reagent controls'. Script also includes normalising reads resulting in relative abundancse. Mapping files uses the 'index' present on the metadata file <br>
> [Metadata-to-reads.R](https://github.com/nteng22/Imperial-Molyneaux/blob/main/COVID-BAL-OralRinse/Quality-check/Metadata-to-reads.R), joins the metadata to relative abundances. Index is used as the variable to join both datasets (normalised reads with metadata).
