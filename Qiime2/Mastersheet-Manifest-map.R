# Change the working directory accordingly
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT001_Microbiota-analysis/Proprionate/")

library(readxl)
library(tidyverse)
library(dplyr)

# Create a sheet with the sample IDs that you want to run through qiime2.
# Need to rename the column of your sheet to "#SampleID" as this is the column name that will be used to merge the barcodes later.

samples <- readxl::read_xlsx("Original-files/Proprionate-SRA.xlsx", sheet=2)
sampleID <- as.data.frame(samples$library_ID)
sampleID <- dplyr::rename("#SampleID" = colnames(sampleID), sampleID) %>%
  drop_na()

metadata <- read_csv("../16S Master sheet.csv")

unique_manifest <- inner_join(sampleID, metadata) %>%
  select(-"#SampleID") %>%
  relocate(.after = NULL, "sample-id")

write.table(unique_manifest, file = "Original-files/Manifest-map.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
unique_manifest %>%
  group_by(`Level 1`) %>%
  count()

