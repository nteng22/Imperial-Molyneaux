setwd("~/GitHub/Imperial-Molyneaux/POSTCODE/Unfiltered/")
path <- "~/GitHub/Imperial-Molyneaux/POSTCODE/Original-files/"

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
#BiocManager::install("microbiome")
library("microbiome")

metadata <- read_tsv("../Original-files/Manifest-map.txt")

metadata <- column_to_rownames(metadata, var = "sample-id") 
# Need to make the IDs into the rownames so physeq can 'read' it and match it.

list.files(path)
physeq<-qza_to_phyloseq(
  features="../Original-files/filtered-metadata-table.qza",
  tree="../Original-files/merged-midrooted-tree.qza",
  taxonomy= "../Original-files/taxonomy-merged-final.qza")
physeq 

mapping=metadata
sample_data(physeq)<-mapping

gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq))) # All samples match
setdiff(rownames(mapping), sample_names(physeq)) #mock samples

# Contamination check
ps<-physeq
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=`Diagnosis`)) + geom_point() +
  ylim(lower = 0, upper = 80000)#somewhat rarecurve?

# Filter out taxa thta have less than 1000 reads
filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x >= 1000))
otu_table<-as.data.frame(ps@otu_table) 
library(tibble)
library(dplyr)

# Removing taxa with reads less than 1000 reads
filtered_reads <- prune_taxa(!filter, ps)
filtered_reads #96 taxa less than 1000
summarize_phyloseq(filtered_reads)

sum(taxa_sums(filtered_reads) == 0)

summarize_phyloseq(filtered_reads)

# Extract taxa information
tax <- as(tax_table(ps), "matrix")
tax_df <- as.data.frame(tax)

filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla)

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

sum(taxa_sums(prunedSet)==0)

summarize_phyloseq(prunedSet)

# Set function to normalise samples
normalizeSample = function(x) {
  x/sum(x)
}

Controls_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(Controls_relative)
OTU1 = as(otu_table(Controls_relative), "matrix")
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf, "OTUdf.csv")

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)
write.csv(TAXdf, "tax_table.csv")

Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum') #7 phyla most likely a contaminant
otu_table<-as.data.frame(Controls_Phylum@otu_table)
write.table(otu_table,file="Phylum/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Class <- aggregate_taxa(Controls_relative, 'Class')
otu_table<-as.data.frame(Controls_Class@otu_table)
write.table(otu_table,file="Class/Class-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Family <- aggregate_taxa(Controls_relative, 'Family')
otu_table<-as.data.frame(Controls_Family@otu_table)
write.table(otu_table,file="Family/Family-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
otu_table<-as.data.frame(Controls_Genus@otu_table)
write.table(otu_table,file="Genus/Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")