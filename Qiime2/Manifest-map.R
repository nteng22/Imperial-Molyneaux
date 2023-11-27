# Set working directory to where your files are found.
# You need to have made an excel sheet with the sample IDs, and other metadata you may want. 
  # Most important one is the sample ID. 

setwd("C://Documents and Settings/mteng/OneDrive - Imperial College London/2023/Experiments/16S database/IPF-stable-progressive//")

library(readxl)
library(tidyverse)
library(dplyr)

samples <- readxl::read_xlsx("HC-IPF_map.xlsx") # Read in that excel sheet with sample IDs
samples$`Sequencing run` <- as.factor(samples$`Sequencing run`)
sampleID <- as.data.frame(samples$`#SampleID`)
sampleID <- dplyr::rename("sample-id" = colnames(sampleID), sampleID)

# Set working directory to where the MSQ maps are found.
# Fibrosis drive, need to change the letter accordingly. 
setwd("Y://Phil/16S Raw Data/16S-fastq-combined/MSQ-maps/")

MSQ75 <- read_tsv("MSQ75-Map.txt")
MSQ75_barcodes <- inner_join(MSQ75, samples, by = "#SampleID")

MSQ81 <- read_tsv("MSQ81-Map.txt")
MSQ81_barcodes <- inner_join(MSQ81, samples, by = "#SampleID")

MSQ82 <- read_tsv("MSQ82-Map.txt")
MSQ82_barcodes <- inner_join(MSQ82, samples, by = "#SampleID")

MSQ91 <- read_tsv("MSQ91-Map.txt")
MSQ91_barcodes <- inner_join(MSQ91, samples, by = "#SampleID")

MSQ92 <- read_tsv("MSQ92-Map.txt")
MSQ92_barcodes <- inner_join(MSQ92, samples, by = "#SampleID")

MSQ104 <- read_tsv("MSQ104-Map.txt")
MSQ104_barcodes <- inner_join(MSQ104, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -SubjID)

MSQ113 <- read_tsv("MSQ113-Map.txt")
MSQ113_barcodes <- inner_join(MSQ113, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -SubjID)

MSQ114 <- read_tsv("MSQ114-Map.txt")
MSQ114_barcodes <- inner_join(MSQ114, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -SubjID)

MSQ139 <- read_tsv("MSQ139-Map.txt") # MSQ139 and 140 have the same mapping file.
MSQ139_barcodes <- inner_join(MSQ139, samples, by = "#SampleID")%>%
  select(-Primer_Plate, -LANE, -SubjID)

manifest <- rbind(MSQ75_barcodes,
                  MSQ81_barcodes,
                  MSQ82_barcodes,
                  MSQ91_barcodes,
                  MSQ92_barcodes,
                  MSQ104_barcodes,
                  MSQ113_barcodes,
                  MSQ114_barcodes,
                  MSQ139_barcodes)

colnames(manifest)
colnames(manifest) <- c("sample-id", #change sample-id column to avoid SQL keywords
                        "BarcodeSequence",
                        "LinkerPrimerSequence",
                        "Host",
                        "Amp_Well_Plate",
                        "MSQ",
                        "PI",
                        "Description",
                        "Level1",
                        "Sequencing-run",
                        "Diagnosis")

unique_manifest <- manifest[!duplicated(manifest$`sample-id`), ]
duplicated_manifest <- manifest[duplicated(manifest$`sample-id`),] # Check for duplicated sample-ids

manifest_ID <- manifest$`sample-id`

# Manually need to check this data frame to see if all samples add up
name_check <- full_join(sampleID, unique_manifest)

# Change working directory back to your original folder. 
setwd("C://Documents and Settings/mteng/OneDrive - Imperial College London/2023/Experiments/Bavi-Project/")

write.table(unique_manifest, file = "Manifest-map.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
