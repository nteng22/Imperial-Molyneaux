setwd("~/GitHub/Imperial-Molyneaux/POSTCODE/Unfiltered/")

#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("vegan")
library("vegan")
#install.packages("reshape2")
library("reshape2")

metadata <- read_tsv("../Original-files/Manifest-map.txt")
metadata <- metadata %>%
  select(`sample-id`,Description,Diagnosis)
metadata <- metadata[!duplicated(metadata$`sample-id`), ] #remove duplicated sample id

normalised_data <- read_tsv("Phylum/Phylum-relative-abundance.txt")
normalised_data <- column_to_rownames(normalised_data, "...1")
normalised_data <- as.data.frame(t(normalised_data))
normalised_data_filtered <- normalised_data[rowSums(normalised_data) != 0,]
no_hits <- normalised_data[rowSums(normalised_data) == 0,] # Remove samples with a relative abundance of 0. 
normalised_data_filtered <- rownames_to_column(normalised_data_filtered, var = "sample-id")
normalised_data_filtered <- inner_join(metadata, normalised_data_filtered, by = "sample-id", relationship = "one-to-one")
write_csv(normalised_data_filtered, "Phylum/Phylum-normalised-unfiltered-metadata.csv")

normalised_data <- read_tsv("Class/Class-relative-abundance.txt")
normalised_data <- column_to_rownames(normalised_data, "...1")
normalised_data <- as.data.frame(t(normalised_data))
normalised_data_filtered <- normalised_data[rowSums(normalised_data) != 0,]
no_hits <- normalised_data[rowSums(normalised_data) == 0,]
normalised_data_filtered <- rownames_to_column(normalised_data_filtered, var = "sample-id")
normalised_data_filtered <- inner_join(metadata, normalised_data_filtered, by = "sample-id", relationship = "one-to-one")
write_csv(normalised_data_filtered, "Class/Class-normalised-unfiltered-metadata.csv")

normalised_data <- read_tsv("Family/Family-relative-abundance.txt")
normalised_data <- column_to_rownames(normalised_data, "...1")
normalised_data <- as.data.frame(t(normalised_data))
normalised_data_filtered <- normalised_data[rowSums(normalised_data) != 0,]
no_hits <- normalised_data[rowSums(normalised_data) == 0,]
normalised_data_filtered <- rownames_to_column(normalised_data_filtered, var = "sample-id")
normalised_data_filtered <- inner_join(metadata, normalised_data_filtered, by = "sample-id", relationship = "one-to-one")
write_csv(normalised_data_filtered, "Family/Family-normalised-unfiltered-metadata.csv")

normalised_data <- read_tsv("Genus/Genus-relative-abundance.txt")
normalised_data <- column_to_rownames(normalised_data, "...1")
normalised_data <- as.data.frame(t(normalised_data))
normalised_data_filtered <- normalised_data[rowSums(normalised_data) != 0,]
no_hits <- normalised_data[rowSums(normalised_data) == 0,]
normalised_data_filtered <- rownames_to_column(normalised_data_filtered, var = "sample-id")
normalised_data_filtered <- inner_join(metadata, normalised_data_filtered, by = "sample-id", relationship = "one-to-one")
write_csv(normalised_data_filtered, "Genus/Genus-normalised-unfiltered-metadata.csv")
