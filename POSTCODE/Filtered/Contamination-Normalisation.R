setwd("~/GitHub/Imperial-Molyneaux/POSTCODE/Filtered/")
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

gplots::venn(list(mapping=rownames(mapping), sample_names(physeq))) # All samples match
setdiff(rownames(mapping), sample_names(physeq)) #mock samples

# Contamination check
ps<-physeq
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=`Diagnosis`)) + geom_point() +
  ylim(lower = 0, upper = 80000)#somewhat rarecurve?

sample_data(ps)$is.neg <- sample_data(ps)$Diagnosis=="Negative Control" #Only considered reagent control as negative
contamdf.prev <- isContaminant(ps, method=c("prevalence"), neg="is.neg", threshold=0.1) 
  # Don't have DNA quantity so therefore can't do frequency
  # Prevalence "Contaminants are identified by increased prevalence in negative controls"
  # Gets data frame with a list of potential contaminants

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Diagnosis =="Negative Control", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$Diagnosis == "Negative Control", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") # Need more info

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.prev$contaminant),]

write.table(contaminants,file="../Original-files/Contaminants-list.txt", col.names=NA, row.names=T,sep="\t")

# This leaves only the contaminants
ps.contam <- prune_taxa(contamdf.prev$contaminant, ps)
summarize_phyloseq(ps.contam)

#Look at distribution
plot_bar(ps.contam)
ggsave("../Filtered/Contaminated-samples.pdf", height=150, width=400, units=c("mm"), device="pdf")

# See what the reads higher than 1000 are in the contaminants
filter <- phyloseq::genefilter_sample(ps.contam, filterfun_sample(function(x) x >= 1000))
ps.contam.1k <- prune_taxa(filter, ps.contam)
otu_table<-as.data.frame(ps.contam@otu_table)
library(tibble)
library(dplyr)

#Check if tax is now a dataframe
tax_df <- as.data.frame(tax)
is.data.frame(tax_df)

contam<-left_join(rownames_to_column(otu_table), (rownames_to_column(tax_df)))
write.table(contam,file="contamination.txt", col.names=NA, row.names=T,sep="\t")

# Removing potential contaminants from reads
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam
summarize_phyloseq(ps.noncontam)

POSTCODE <-subset_samples(ps.noncontam, !Diagnosis == "Negative Control") 
sum(taxa_sums(POSTCODE) == 0) # 111 samples have no microbiota taxa? 

summarize_phyloseq(POSTCODE)

POSTCODE

# Filter out unclassified phyla, its useless
filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla)

ps1 = subset_taxa(POSTCODE, !(!Phylum %in% filterPhyla))
ps1
summarize_phyloseq(ps1)

POSTCODE_decontaminated <-ps1 #2746 taxa

badTaxa <- contam$rowname #fits with identified number of contaminants, 66
goodTaxa <- setdiff(taxa_names(POSTCODE_decontaminated), badTaxa) # 2746
POSTCODE_decontaminated_final <- prune_taxa(goodTaxa, POSTCODE_decontaminated)
POSTCODE_decontaminated_final
summarize_phyloseq(POSTCODE_decontaminated_final)

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(POSTCODE_decontaminated_final)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), POSTCODE_decontaminated_final)
prunedSet
summarize_phyloseq(prunedSet)

# Set function to normalise samples
normalizeSample = function(x) {
  x/sum(x)
}

POSTCODE_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(POSTCODE_relative)
OTU1 = as(otu_table(POSTCODE_relative), "matrix")
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf, "OTUdf.csv")

TAXdf = as(tax_table(POSTCODE_relative), "matrix")
TAXdf = as.data.frame(TAXdf)
write.csv(TAXdf, "tax_table.csv")

POSTCODE_Phylum <- aggregate_taxa(POSTCODE_relative, 'Phylum')
otu_table<-as.data.frame(POSTCODE_Phylum@otu_table)
write.table(otu_table,file="Phylum/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

POSTCODE_Class <- aggregate_taxa(POSTCODE_relative, 'Class')
otu_table<-as.data.frame(POSTCODE_Class@otu_table)
write.table(otu_table, file="Class/Class-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

POSTCODE_Family <- aggregate_taxa(POSTCODE_relative, 'Family')
otu_table<-as.data.frame(POSTCODE_Family@otu_table)
write.table(otu_table,file="Family/Family-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

POSTCODE_Genus <- aggregate_taxa(POSTCODE_relative, 'Genus')
otu_table<-as.data.frame(POSTCODE_Genus@otu_table)
write.table(otu_table,file="Genus/Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")
