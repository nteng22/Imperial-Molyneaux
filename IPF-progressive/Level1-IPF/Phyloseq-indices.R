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

setwd(dir = "~/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF/")
path <- "~/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF/"
list.files(path)

physeq<-qza_to_phyloseq("filtered-metadata-table.qza",
                        "merged-midrooted-tree.qza",
                        "taxonomy-merged-final.qza",
                        "Manifest-map.txt")

mapping = readxl::read_xlsx("HC-IPF_map.xlsx")
mapping <- dplyr::rename("index" = "#SampleID", mapping)
mapping$Level1 <- gsub("Controls", "Control", mapping$Level1)
mapping <- mapping[!duplicated(mapping$index), ]
mapping <- column_to_rownames(mapping, var = "index")

gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq))) # All samples match
setdiff(rownames(mapping), sample_names(physeq))

sample_data(physeq) <- mapping
physeq
summarize_phyloseq(physeq)

################################################################### ALPHA-DIVERSITY ######################################################################
#To calculate alpha-diversity use unfiltered (1,844 ASVs) feature table

##Alpha-Plot
p = plot_richness(physeq, x="Level1", color="Level1", measures=c("Observed","Shannon", "Chao1"))

print(p)
p + geom_boxplot(data = p$data, aes(x = Level1, y = value, color = NULL), alpha = 0.1) +
  scale_colour_manual(values=c("darkorange", "darkblue","steelblue", "orchid")) + 
  theme_classic()+ theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank())

##Beta-plot
ordu = ordinate(physeq, "PCoA", "unifrac", weighted=TRUE)
#Temporary fix to colour bug is to add a "dummy variable" in sample_data:
sample_data(physeq)[ , 2] <- sample_data(physeq)[ ,1]
p = plot_ordination(physeq, ordu, color="Level1") + geom_point(size = 3) +
  ggtitle("PCoA on weighted-UniFrac distance") + (scale_colour_brewer(type="qual", palette="Set1"))
print(p)
myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_colour_manual(values=c("darkorange", "darkblue", "steelblue", "orchid"))
