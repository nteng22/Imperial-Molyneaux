setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse//")

#install.packages("phyloseq")
library("phyloseq")
#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("vegan")
library("vegan")
#install.packages("reshape2")
library("reshape2")
#remotes::install_github("jbisanz/qiime2R")
library("qiime2R")

sampleID <- read_csv("Quality-check/Qiime2-Genus.csv") %>%
  select(index)
index <- sampleID

sampleID$index <- gsub("FL.POST.", "BAL-control", sampleID$index)
sampleID$index <- gsub("NEG.OR.CON", "Oral-control", sampleID$index)
sampleID$index <- gsub("NEG.OR.", "Oral-control", sampleID$index)
sampleID$index <- gsub("OR.POST.", "Oral-disease", sampleID$index)
sampleID$index <- gsub("POST.", "BAL-disease", sampleID$index)
sampleID$index <- gsub("MOC.", "EXCLUDE-EXCLUDE", sampleID$index) # Mock group, separate analysis

metadata <- separate(sampleID, index,
                     sep = "(?<=[A-Za-z])(?=[0-9])",
                     into = c("Sample-type", "patient_ID"),
                     extra = "merge")

metadata$`Sample-type` <- gsub("FL.CON", "BAL-control", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("RC", "Reagent-control", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("ILDCON", "BAL-healthy", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("OR.CON", "Oral-disease", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("POST", "POST", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("PH", "EXCLUDE", metadata$`Sample-type`) 
# Change PH to BAL-disease if want to include
# Change PH to EXCLUDE if want to exclude

metadata <- separate(metadata, `Sample-type`,
                     sep = "-",
                     into = c("Sample-type", "Status"))
metadata$`Sample-type` <- as.factor(metadata$`Sample-type`)
metadata$Status <- as.factor(metadata$Status)

metadata$patient_ID <- gsub("[A-Za-z]", "", metadata$patient_ID)
metadata$patient_ID <- gsub("[..]", "", metadata$patient_ID)
# Note that the 02.002 nomenclature has merged
# Decided to merge this as we don't have any other information re. Site.
# Keeping nomenclature like this highlights it's from a different Study. 

# As an added safety adding the manifest sampleID to the metadata
metadata <- cbind(index,metadata)

write_csv(metadata, "Quality-check/metadata.csv")
