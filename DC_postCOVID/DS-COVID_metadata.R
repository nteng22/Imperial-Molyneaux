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

sampleID <- read_csv("Qiime2-Genus.csv") %>%
  select(index)
index <- read_csv("Qiime2-Genus.csv") %>%
  select(index)

sampleID$index <- gsub("FL.POST.", "BAL-control", sampleID$index)
sampleID$index <- gsub("NEG.OR.CON", "Oral-control", sampleID$index)
sampleID$index <- gsub("NEG.OR.", "Oral-control", sampleID$index)
sampleID$index <- gsub("OR.POST.", "Oral-disease", sampleID$index)
sampleID$index <- gsub("POST.", "BAL-disease", sampleID$index)
sampleID$index <- gsub("MOC.", "EXCLUDE-EXCLUDE", sampleID$index)

metadata <- separate(sampleID, index,
                    sep = "(?<=[A-Za-z])(?=[0-9])",
                    into = c("Sample-type", "patient_ID"),
                    extra = "merge")

metadata$`Sample-type` <- gsub("FL.CON", "BAL-control", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("RC", "Reagent-control", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("ILDCON", "BAL-disease", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("OR.CON", "Oral-disease", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("POST", "POST", metadata$`Sample-type`)
metadata$`Sample-type` <- gsub("PH", "EXCLUDE", metadata$`Sample-type`)

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

reads <- read_csv("Qiime2-Genus.csv") %>%
  select(!index)
genus_reads <- colnames(reads)
colnames_genus <- str_extract(genus_reads, "g__[A-Za-z]+")
colnames_genus <- gsub("g__", "", colnames_genus)
colnames_genus[is.na(colnames_genus)] <- "Unclassified"

colnames(reads) <- colnames_genus

Genus_metadata <- cbind(metadata, reads)
Genus_metadata$`Sample-type` <- gsub("EXCLUDE", NA, Genus_metadata$`Sample-type`)
Genus_metadata <- drop_na(Genus_metadata)

write_csv(Genus_metadata, "Genus-OTU.csv")

reads <- read_csv("Qiime2-Family.csv") %>%
  select(!index)
family_reads <- colnames(reads)
colnames_family <- str_extract(family_reads, "f__[A-Za-z]+")
colnames_family <- gsub("f__", "", colnames_family)
colnames_family[is.na(colnames_family)] = "Unclassified"
colnames(reads) <- colnames_family
Family_metadata <- cbind(metadata,reads)
Family_metadata$`Sample-type` <- gsub("EXCLUDE", NA, Family_metadata$`Sample-type`)
Family_metadata <- drop_na(Family_metadata)

write_csv(Family_metadata, "Family-OTU.csv")

