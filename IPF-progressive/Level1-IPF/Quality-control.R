setwd("~/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF/")
setwd("~/Documents/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF") # Personal device

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
#BiocManager::install("microbiome")
library("microbiome")

metadata <- read_tsv("Manifest-map.txt")

metadata <- column_to_rownames(metadata, var = "sample-id") 
# Need to make the IDs into the rownames so physeq can 'read' it and match it.

physeq<-qza_to_phyloseq(
  features="filtered-metadata-table.qza",
  tree="merged-midrooted-tree.qza",
  taxonomy="taxonomy-merged-final.qza")
physeq 

mapping=metadata
sample_data(physeq)<-mapping

gplots::venn(list(mapping=rownames(mapping), sample_names(physeq))) # All samples match
setdiff(rownames(mapping), sample_names(physeq))

# Contamination check
ps<-physeq
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, colour = Level1)) + geom_point() +  ylim(lower = 0, upper = 200000) #somewhat rarecurve?

sample_data(ps)$is.neg <- sample_data(ps)$`Level1`=="Negative" #Only considered reagent control as negative
contamdf.prev <- isContaminant(ps, method=c("prevalence"), neg="is.neg", threshold=0.1) 
# Don't have DNA quantity so therefore can't do frequency
# Prevalence "Contaminants are identified by increased prevalence in negative controls"
# Gets data frame with a list of potential contaminants

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$`Level1`=="Negative", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$`Level1`=="Negative", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") # Need more info

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.prev$contaminant),] # 24 contaminants

write.table(contaminants,file="Contaminants-list.txt", col.names=NA, row.names=T,sep="\t")

# This leaves only the contaminants
ps.contam <- prune_taxa(contamdf.prev$contaminant, ps)
summarize_phyloseq(ps.contam)

#Look at distribution
plot_bar(ps.contam)
ggsave("Contamination-in-samples.pdf", height=150, width=400, units=c("mm"), device="pdf")

# Filter out samples about 1000 reads
filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x >= 1000))
otu_table<-as.data.frame(ps@otu_table) # 71 are less than 1000
library(tibble)
library(dplyr)

# Removing low reads from data
filtered_reads <- prune_taxa(!filter, ps)
filtered_reads
summarize_phyloseq(filtered_reads)

sum(taxa_sums(filtered_reads) == 0) # Pruned correctly

summarize_phyloseq(filtered_reads)

# Filter out unclassified phyla, its useless
tax <- as(tax_table(ps), "matrix")
tax_df <- as.data.frame(tax)

filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla) # 10 identified Phyla

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
write.csv(OTUdf, "OTUdf.csv", )

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
