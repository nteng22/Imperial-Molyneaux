setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/")
setwd("~/Documents/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/") # Personal device

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
#BiocManager::install("microbiome")
library("microbiome")

metadata <- read_csv("Quality-check/metadata.csv") %>%
  select(index:patient_ID) %>%
  filter(!`Sample-type`=="EXCLUDE")

metadata <- column_to_rownames(metadata, var = "index") 
# Need to make the IDs into the rownames so physeq can 'read' it and match it.

physeq<-qza_to_phyloseq(
  features="Quality-check/table.qza",
  tree="Quality-check/rooted-tree_quality.qza",
  taxonomy="Quality-check/taxonomy.qza")
physeq 

mapping=metadata
sample_data(physeq)<-mapping

gplots::venn(list(mapping=rownames(mapping), sample_names(physeq))) # All samples match

# Contamination check
ps<-physeq
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=`Status`)) + geom_point() #somewhat rarecurve?

# Filter out samples about 1000 reads
filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x >= 1000))
otu_table<-as.data.frame(ps@otu_table) # 178 are less than 1000
library(tibble)
library(dplyr)

# Removing low reads from data
filtered_reads <- prune_taxa(!filter, ps)
filtered_reads
summarize_phyloseq(filtered_reads)

sum(taxa_sums(filtered_reads) == 0) # 366 samples equal to 0

summarize_phyloseq(filtered_reads)

# Filter out unclassified phyla, its useless
tax <- as(tax_table(ps), "matrix")
tax_df <- as.data.frame(tax)

filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla) # 26 identified Phyla

ps1 = subset_taxa(filtered_reads, !(!Phylum %in% filterPhyla))
ps1 # Only keep the Phyla in filterPhyla in the filtered reads dataset

# Check if no reads less than 1000
sum(phyloseq::genefilter_sample(ps1, filterfun_sample(function(x) x >= 1000)))

# Check if there are unique Phyla names
unique(as.data.frame(as(tax_table(ps1), "matrix"))$Phylum)

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(ps1)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), ps1)
prunedSet
summarize_phyloseq(prunedSet)

# Set function to normalise samples
normalizeSample = function(x) {
  x/sum(x)
}

Controls_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(Controls_relative)
OTU1 = as(otu_table(Controls_relative), "matrix")
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf, "Unfiltered/OTUdf.csv", )

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)
write.csv(TAXdf, "Unfiltered/tax_table.csv")

Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum') #7 phyla most likely a contaminant
otu_table<-as.data.frame(Controls_Phylum@otu_table)
tax<-as.data.frame(tax)
Controls_Phylum_final <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(Controls_Phylum_final,file="Unfiltered/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Class <- aggregate_taxa(Controls_relative, 'Class')
otu_table<-as.data.frame(Controls_Class@otu_table)
tax<-as.data.frame(tax)
Controls_Class_final <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(Controls_Class_final,file="Unfiltered/Class-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
otu_table<-as.data.frame(Controls_Genus@otu_table)
tax<-as.data.frame(tax)
Controls_Genus_final <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(Controls_Genus_final,file="Unfiltered//Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")
