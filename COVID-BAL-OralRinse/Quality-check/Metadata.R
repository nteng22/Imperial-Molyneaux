setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse//")
setwd("~/Documents/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/")

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
metadata$`Sample-type` <- gsub("OR.CON", "Oral-healthy", metadata$`Sample-type`)
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

duplicates <- metadata[duplicated(metadata %>% select(!index)) &
                         duplicated(metadata %>% select(!index), fromLast = TRUE),]
unique_metadata <- metadata %>% # Manually remove the samples I do not want to keep, based of nanodrop and library size
  filter(!index == "ILDCON1064a.BAL" &
           !index == "ILDCON1065a.BAL" &
           !index == "ILDCON1066a.BAL" &
           !index == "ILDCON1069a.BAL" &
           !index == "ILDCON1071a.BAL" &
           !index == "ILDCON1076a.BAL" &
           !index == "OR.CON1064.OralRinse" &
           !index == "ILDCON1069.BAL" & 
           !index == "OR.POST.01.005a.OralRinse" &
           !index == "OR.POST.01.007a.OralRinse" &
           !index == "OR.POST.01.013a.OralRinse" &
           !index == "POST.01.003.BAL" &
           !index == "POST.01.004.BAL" &
           !index == "POST.01.007a.BAL" &
           !index == "POST.01.007.BAL" &
           !index == "POST.01.011a.BAL" &
           !index == "POST.01.013.BAL")

write_csv(unique_metadata, "Quality-check/metadata.csv")
