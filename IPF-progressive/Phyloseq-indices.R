###### BRU vs PROFILE ########

#Load required packages 
library(devtools)
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(biomformat)
library(ape)
library(scales)
library(grid)
library(plyr)
library(vegan)
library(lfda)
library(qiime2R)
library(yaml)
library(tidyverse)
library(microbiome)
library(knitr)


################################################################### IMPORT DATA ##########################################################################

setwd(dir = "~/GitHub/Imperial-Molyneaux/IPF-progressive/")
path <- "~/GitHub/Imperial-Molyneaux/IPF-progressive/"
list.files(path)

physeq<-qza_to_phyloseq("filtered-metadata-table.qza",
                        "merged-midrooted-tree.qza",
                        "taxonomy-merged-final.qza",
                        "Manifest-map.txt")

mapping = readxl::read_xlsx("HC-IPF_map.xlsx")
sample_data(physeq) <- mapping$`#SampleID`
physeq
keep = c("PROFILE-IPF", "BRU-IPF")
physeq2<- subset_samples(physeq, Diagnosis %in% keep)

#Summarise contents of phyloseq object
summarize_phyloseq(physeq2)

################################################################### ALPHA-DIVERSITY ######################################################################
#To calculate alpha-diversity use unfiltered (1,844 ASVs) feature table

##Alpha-Plot

#BRU-IPF vs PROFILE-IPF 
p = plot_richness(physeq, x="Level1", color="Level1", measures=c("Observed","Shannon"))
newSTorder = c("BRU-IPF", "PROFILE-IPF")
p$data$Diagnosis <- as.character(p$data$Diagnosis)
p$data$Diagnosis <- factor(p$data$Diagnosis)
print(p)
p + geom_boxplot(data = p$data, aes(x = Diagnosis, y = value, color = NULL), alpha = 0.1) + scale_colour_manual(values=c("darkorange", "darkblue",
                                                                                                                         "steelblue", "orchid")) + theme_classic()+ theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank())

################################################################### IMPORT DATA ##########################################################################

setwd(dir = "/Users/rachele_invernizzi/Desktop/IPF")
path <- "/Users/rachele_invernizzi/Desktop/IPF"
list.files(path)

physeq<-qza_to_phyloseq("0.005%-filtered-IPF-table.qza","0.005%-rooted-IPF-tree.qza","taxonomy-0.005%-filtered-IPF-final.qza", "mapping.tsv")

#Read in mapping file and filter to only CHP paper samples (IPF_mapping.csv) or lymphocyte subgroups (CHP_mapping_final_lymphocytes.csv)

mapping = read.csv("IPF_mapping.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))
sample_data(physeq) <- mapping
physeq
keep = c("PROFILE-IPF", "BRU-IPF")
physeq2<- subset_samples(physeq, Diagnosis %in% keep)

#Summarise contents of phyloseq object
summarize_phyloseq(physeq2)

#Weighted BRU-IPF vs PROFILE-IPF 
ordu = ordinate(physeq2, "PCoA", "unifrac", weighted=TRUE)
#Temporary fix to colour bug is to add a "dummy variable" in sample_data:
sample_data(physeq2)[ , 2] <- sample_data(physeq2)[ ,1]
p = plot_ordination(physeq2, ordu, color="Diagnosis") + geom_point(size = 5) + ggtitle("PCoA on weighted-UniFrac distance") + (scale_colour_brewer(type="qual", palette="Set1"))
newSTorder = c("BRU-IPF", "PROFILE-IPF")
p$data$Diagnosis <- as.character(p$data$Diagnosis)
p$data$Diagnosis <- factor(p$data$Diagnosis, levels=newSTorder)
print(p)
myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + theme(panel.grid.major = element_blank(),
                                                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_colour_manual(values=c("darkorange", "darkblue"))



