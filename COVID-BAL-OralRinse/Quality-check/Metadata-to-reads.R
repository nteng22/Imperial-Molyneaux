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

metadata <-read_csv("metadata.csv")

### Genus ###
reads_genus <- read_tsv("Genus-relative-abundance.txt")
reads_genus <- reads_genus %>%
  select(rowname:ncol(reads_genus))
reads_genus <- column_to_rownames(reads_genus, var = "rowname")
reads_genus <- as.data.frame(t(reads_genus))
reads_genus <- rownames_to_column(reads_genus, var = "index")

reads_genus_metadata <- left_join(metadata, reads_genus)
reads_genus_metadata$`Sample-type` <- gsub("EXCLUDE", NA, reads_genus_metadata$`Sample-type`)
reads_genus_metadata <- drop_na(reads_genus_metadata) # lose quite a number due to low reads

write_csv(reads_genus_metadata, "Genus-normalised-metadata.csv")

### Class ###
reads_class <- read_tsv("Class-relative-abundance.txt")
reads_class <- reads_class %>%
  select(rowname:ncol(reads_class))
reads_class <- column_to_rownames(reads_class, var = "rowname")
reads_class <- as.data.frame(t(reads_class))
reads_class <- rownames_to_column(reads_class, var = "index")

reads_class_metadata <- left_join(metadata, reads_class)
reads_class_metadata$`Sample-type` <- gsub("EXCLUDE", NA, reads_class_metadata$`Sample-type`)
reads_class_metadata <- drop_na(reads_class_metadata) # lose quite a number due to low reads

write_csv(reads_class_metadata, "Class-normalised-metadata.csv")

### Phyla ###
reads_phylum <- read_tsv("Phylum-relative-abundance.txt")
reads_phylum <- reads_phylum %>%
  select(rowname:ncol(reads_phylum))
reads_phylum <- column_to_rownames(reads_phylum, var = "rowname")
reads_phylum <- as.data.frame(t(reads_phylum))
reads_phylum <- rownames_to_column(reads_phylum, var = "index")

reads_phylum_metadata <- left_join(metadata, reads_phylum)
reads_phylum_metadata$`Sample-type` <- gsub("EXCLUDE", NA, reads_phylum_metadata$`Sample-type`)
reads_phylum_metadata <- drop_na(reads_phylum_metadata) # lose quite a number due to low reads

write_csv(reads_phylum_metadata, "Phylum-normalised-metadata.csv")
