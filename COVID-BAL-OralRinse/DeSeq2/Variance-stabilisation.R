library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
library("microbiome")
setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/")

metadata <- read_csv("Quality-check/metadata.csv") %>%
  select(index:patient_ID) %>%
  filter(!`Sample-type`=="EXCLUDE")
metadata <- column_to_rownames(metadata, var = "index")
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
sum(taxa_sums(Controls) == 0) # 498 samples equal to 0
summarize_phyloseq(Controls)
Controls

# Filter out unclassified phyla, its useless
filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla) # 27 identified Phyla
ps1 = subset_taxa(Controls, !(!Phylum %in% filterPhyla))
ps1
summarize_phyloseq(ps1)
Controls_decontaminated <-ps1 #2758 taxa
badTaxa <- contam$rowname #fits with identified number of contaminants, 54
goodTaxa <- setdiff(taxa_names(Controls_decontaminated), badTaxa) # 2758
Controls_decontaminated_final <- prune_taxa(goodTaxa, Controls_decontaminated)
Controls_decontaminated_final
summarize_phyloseq(Controls_decontaminated_final)

Controls_decontaminated_phylum <- aggregate_taxa(Controls_decontaminated_final, 'Phylum')
df <- as.data.frame(Controls_decontaminated_phylum@otu_table)
# Minimum sequencing depth
min(sample_sums(Controls_decontaminated_final))

# Rarecurve
rarecurve(t(as.data.frame(Controls_decontaminated_phylum@otu_table)),
          step = 200,
          sample = 20000,
          col = "blue",
          cex = 0.4,
          label = FALSE,
          lwd = 0.001,
          ylab = "No. of phyla",
          xlab = "No. of sequences"
          #xlim = c(0, 250000)
)


  
### Phyla level ###
library("DESeq2")
count_data <- as.data.frame(Controls_decontaminated_phylum@otu_table)
count_data <- count_data + 1
coldata <- metadata

nrow(coldata) == ncol(count_data) # Condition 1 for DESeq
all(colnames(count_data) %in% rownames(coldata)) # Condition 2 for DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ 1)
# Variance Stabilizing Transformation
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "local")
varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
getVarianceStabilizedData(dds)
# Turn it back into counts.
prunedSet <- counts(dds, normalized = TRUE)
# Turn it back into a dataframe.
prunedSet <- data.frame(prunedSet) -1

write.csv(prunedSet, "DeSeq2/Phylum-normalisedDeseq2.csv", row.names = TRUE)

### Class level ###
Controls_decontaminated_Class <- aggregate_taxa(Controls_decontaminated_final, 'Class')

count_data <- as.data.frame(Controls_decontaminated_Class@otu_table)
count_data <- count_data + 1
coldata <- metadata

nrow(coldata) == ncol(count_data) # Condition 1 for DESeq
all(colnames(count_data) %in% rownames(coldata)) # Condition 2 for DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ 1)
# Variance Stabilizing Transformation
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "local")
varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
getVarianceStabilizedData(dds)
# Turn it back into counts.
prunedSet <- counts(dds, normalized = TRUE)
# Turn it back into a dataframe.
prunedSet <- data.frame(prunedSet) -1

write.csv(prunedSet, "DeSeq2/Class-normalisedDeseq2.csv", row.names = TRUE)

### Genus level ###
Controls_decontaminated_Genus <- aggregate_taxa(Controls_decontaminated_final, 'Genus')

count_data <- as.data.frame(Controls_decontaminated_Genus@otu_table)
count_data <- count_data + 1
coldata <- metadata

nrow(coldata) == ncol(count_data) # Condition 1 for DESeq
all(colnames(count_data) %in% rownames(coldata)) # Condition 2 for DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ 1)
# Variance Stabilizing Transformation
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "local")
varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
getVarianceStabilizedData(dds)
# Turn it back into counts.
prunedSet <- counts(dds, normalized = TRUE)
# Turn it back into a dataframe.
prunedSet <- data.frame(prunedSet) -1

write.csv(prunedSet, "DeSeq2/Genus-normalisedDeseq2.csv", row.names = TRUE)
