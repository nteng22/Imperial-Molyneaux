setwd("/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/DeSeq2/")

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

metadata <-read_csv("../Quality-check/metadata.csv")
metadata <- metadata %>%
  filter(!Status == "control") %>%
  filter(!Status == "EXCLUDE")

### Genus ###
reads_genus <- read_csv("Genus-normalisedDeseq2.csv")

reads_genus <- column_to_rownames(reads_genus, var = "...1")
reads_genus <- as.data.frame(t(reads_genus))
reads_genus <- rownames_to_column(reads_genus, var = "index")

reads_genus_metadata <- inner_join(metadata, reads_genus, by = "index")

write_csv(reads_genus_metadata, "Genus-DeSeq2normalised-metadata.csv")

relabundance <- read_csv("Genus-normalisedDeseq2.csv")
relabundance <- column_to_rownames(relabundance, var = "...1")
relabundance <- as.data.frame(t(relabundance))
relabundance <- relabundance/rowSums(relabundance) * 100
rowSums(relabundance)
relabundance <- rownames_to_column(relabundance, var = "index")
relabundance_metadata <- inner_join(metadata, relabundance)

write_csv(relabundance_metadata, "Genus-DeSeq2normalised-relabundance-metadata.csv")

### Class ###
reads_Class <- read_csv("Class-normalisedDeseq2.csv")

reads_Class <- column_to_rownames(reads_Class, var = "...1")
reads_Class <- as.data.frame(t(reads_Class))
reads_Class <- rownames_to_column(reads_Class, var = "index")

reads_Class_metadata <- inner_join(metadata, reads_Class, by = "index")
reads_Class_metadata$`Sample-type` <- gsub("EXCLUDE", NA, reads_Class_metadata$`Sample-type`)
reads_Class_metadata <- drop_na(reads_Class_metadata) # lose quite a number due to low reads

write_csv(reads_Class_metadata, "Class-DeSeq2normalised-metadata.csv")

relabundance <- read_csv("Class-normalisedDeseq2.csv")
relabundance <- column_to_rownames(relabundance, var = "...1")
relabundance <- as.data.frame(t(relabundance))
relabundance <- relabundance/rowSums(relabundance) * 100
rowSums(relabundance)
relabundance <- rownames_to_column(relabundance, var = "index")
relabundance_metadata <- inner_join(metadata, relabundance)

write_csv(relabundance_metadata, "Class-DeSeq2normalised-relabundance-metadata.csv")
### Phyla ###
reads_Phylum <- read_csv("Phylum-normalisedDeseq2.csv")

reads_Phylum <- column_to_rownames(reads_Phylum, var = "...1")
reads_Phylum <- as.data.frame(t(reads_Phylum))
reads_Phylum <- rownames_to_column(reads_Phylum, var = "index")

reads_Phylum_metadata <- inner_join(metadata, reads_Phylum, by = "index")
reads_Phylum_metadata$`Sample-type` <- gsub("EXCLUDE", NA, reads_Phylum_metadata$`Sample-type`)
reads_Phylum_metadata <- drop_na(reads_Phylum_metadata) # lose quite a number due to low reads

write_csv(reads_Phylum_metadata, "Phylum-DeSeq2normalised-metadata.csv")

relabundance <- read_csv("Phylum-normalisedDeseq2.csv")
relabundance <- column_to_rownames(relabundance, var = "...1")
relabundance <- as.data.frame(t(relabundance))
relabundance <- relabundance/rowSums(relabundance) * 100
rowSums(relabundance)
relabundance <- rownames_to_column(relabundance, var = "index")
relabundance_metadata <- inner_join(metadata, relabundance)

write_csv(relabundance_metadata, "Phylum-DeSeq2normalised-relabundance-metadata.csv")
