setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/")

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
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

sample_data(ps)$is.neg <- sample_data(ps)$`Sample-type`=="Reagent" #Only considered reagent control as negative
contamdf.prev <- isContaminant(ps, method=c("prevalence"), neg="is.neg", threshold=0.1) 
  # Don't have DNA quantity so therefore can't do frequency
  # Prevalence "Contaminants are identified by increased prevalence in negative controls"
  # Gets data frame with a list of potential contaminants

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$`Sample-type` == "Reagent", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$`Sample-type` == "Reagent", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") # Need more info

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.prev$contaminant),]

write.table(contaminants,file="Quality-check/Contaminants-list.txt", col.names=NA, row.names=T,sep="\t")

# This leaves only the contaminants
ps.contam <- prune_taxa(contamdf.prev$contaminant, ps)
summarize_phyloseq(ps.contam)

#Look at distribution
plot_bar(ps.contam)
ggsave("Quality-check/Contaminated-samples.pdf", height=150, width=400, units=c("mm"), device="pdf")

# Filter out samples about 1000 reads
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

Controls <-subset_samples(ps.noncontam) #, `Status`=="control") 
  # Have changed this to all data  as I want to remove the 'reagent' contaminant
sum(taxa_sums(Controls) == 0) # 366 samples equal to 0

summarize_phyloseq(Controls)

Controls

# Filter out unclassified phyla, its useless
filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla) # 27 identified Phyla

ps1 = subset_taxa(Controls, !(!Phylum %in% filterPhyla))
ps1
summarize_phyloseq(ps1)

Controls_decontaminated <-ps1 #2746 taxa

badTaxa <- contam$rowname #fits with identified number of contaminants, 66
goodTaxa <- setdiff(taxa_names(Controls_decontaminated), badTaxa) # 2746
Controls_decontaminated_final <- prune_taxa(goodTaxa, Controls_decontaminated)
Controls_decontaminated_final
summarize_phyloseq(Controls_decontaminated_final)

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(Controls_decontaminated_final)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), Controls_decontaminated_final)
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
write.csv(OTUdf, "Quality-check/OTUdf.csv")

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)
write.csv(TAXdf, "Quality-check/tax_table.csv")

Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum')
otu_table<-as.data.frame(Controls_Phylum@otu_table)
tax<-as.data.frame(tax)
Controls_Phylum_final <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(Controls_Phylum_final,file="Quality-check/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Class <- aggregate_taxa(Controls_relative, 'Class')
otu_table<-as.data.frame(Controls_Class@otu_table)
tax<-as.data.frame(tax)
Controls_Class_final <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(Controls_Class_final,file="Quality-check/Class-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
otu_table<-as.data.frame(Controls_Genus@otu_table)
tax<-as.data.frame(tax)
Controls_Genus_final <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(Controls_Genus_final,file="Quality-check/Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")
