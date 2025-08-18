# Change the working directory accordingly
setwd("C://Documents and Settings/mteng/OneDrive - Imperial College London/2023/Experiments/16S database/IPF-stable-progressive/")

library(readxl)
library(tidyverse)
library(dplyr)

# Create a sheet with the sample IDs that you want to run through qiime2.
# Need to rename the column of your sheet to "#SampleID" as this is the column name that will be used to merge the barcodes later.

samples <- readxl::read_xlsx("HC-IPF_map.xlsx")
samples$`Sequencing run` <- as.factor(samples$`Sequencing run`)
sampleID <- as.data.frame(samples$`#SampleID`)
sampleID <- dplyr::rename("sample-id" = colnames(sampleID), sampleID)


MSQ75 <- read_tsv("../16S database/MSQ-maps/MSQ75-Map.txt") # May need to change the relative path to the MSQ mapping files. 
MSQ75_barcodes <- inner_join(MSQ75, samples, by = "#SampleID") %>%
  select(-Description)

MSQ81 <- read_tsv("../16S database/MSQ-maps/MSQ81-Map.txt")
MSQ81_barcodes <- inner_join(MSQ81, samples, by = "#SampleID")%>%
  select(-Description)

MSQ82 <- read_tsv("../16S database/MSQ-maps/MSQ82-Map.txt")
MSQ82_barcodes <- inner_join(MSQ82, samples, by = "#SampleID")%>%
  select(-Description)

MSQ91 <- read_tsv("../16S database/MSQ-maps/MSQ91-Map.txt")
MSQ91_barcodes <- inner_join(MSQ91, samples, by = "#SampleID")%>%
  select(-Description)

MSQ92 <- read_tsv("../16S database/MSQ-maps/MSQ92-Map.txt")
MSQ92_barcodes <- inner_join(MSQ92, samples, by = "#SampleID")%>%
  select(-Description)

MSQ104 <- read_tsv("../16S database/MSQ-maps/MSQ104-Map.txt")
MSQ104_barcodes <- inner_join(MSQ104, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -SubjID, -Description)

MSQ113 <- read_tsv("../16S database/MSQ-maps/MSQ113-Map.txt")
MSQ113_barcodes <- inner_join(MSQ113, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -SubjID, -Description)

MSQ114 <- read_tsv("../16S database/MSQ-maps/MSQ114-Map.txt")
MSQ114_barcodes <- inner_join(MSQ114, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -SubjID, -Description)

MSQ132 <- read_tsv("../16S database/MSQ-maps/MSQ132-Map.txt")
MSQ132_barcodes <- inner_join(MSQ132, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -`Sample.type`, -SubjID)

MSQ133 <- read_tsv("../16S database/MSQ-maps/MSQ133-Map.txt")
MSQ133_barcodes <- inner_join(MSQ133, samples, by = "#SampleID") %>%
  select(-Primer_Plate, -LANE, -"Sample type", -SubjID)

MSQ139 <- read_tsv("../16S database/MSQ-maps/MSQ139-Map.txt") # MSQ139 and 140 have the same mapping file.
MSQ139_barcodes <- inner_join(MSQ139, samples, by = "#SampleID")%>%
  select(-Primer_Plate, -LANE, -SubjID, -Description)

manifest <- rbind(MSQ75_barcodes,
                  MSQ81_barcodes,
                  MSQ82_barcodes,
                  MSQ91_barcodes,
                  MSQ92_barcodes,
                  MSQ104_barcodes,
                  MSQ113_barcodes,
                  MSQ132_barcodes,
                  MSQ133_barcodes,
                  MSQ114_barcodes,
                  MSQ139_barcodes)

colnames(manifest)
manifest <- dplyr::rename("sample-id" = "#SampleID", manifest)

unique_manifest <- manifest[!duplicated(manifest$`sample-id`), ] # only include unique sample IDs
duplicated_manifest <- manifest[duplicated(manifest$`sample-id`),] # Check for duplicated #SampleIDs

name_check <- full_join(samples, unique_manifest)
# Check that number of rows in samples is the same as unique manifest

write.table(unique_manifest, file = "Original-files/Manifest-map.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
