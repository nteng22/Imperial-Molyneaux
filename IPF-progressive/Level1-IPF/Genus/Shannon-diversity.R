#### Housekeeping ####
setwd("~/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF/Genus")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggplot2")
library("ggplot2")
#install.packages("reshape2")
library("reshape2")
#install.paclages("vegan")
library("vegan")
#install.packages("rstatix")
library("rstatix")

# Read in the saved csv file from previously
normalised_data <- read_csv("Genus-normalised.csv") %>%
  filter(!Level1=="Negative"&!Level1=="Mock")

abundance_table <- normalised_data %>%
  select(5:ncol(normalised_data)) # Only select the reads
rowSums(abundance_table) #check if in percentages
abundance_table <- abundance_table/rowSums(abundance_table)*100

metadata <- normalised_data %>%
  select("sample-id":"Diagnosis")
metadata$Level1 <- gsub("Controls", "Control", metadata$Level1)

#### Species Richness ####
# Calculate species number 
# Margin = 1, by row
species_number <- vegan::specnumber(abundance_table, MARGIN = 1)
species_number <- as.data.frame(species_number)
species_metadata <- cbind(metadata, species_number)
species_metadata$Level1 <- as.factor(species_metadata$Level1)

wilcox.test(species_number ~ Level1, data = species_metadata) #p=0.4

ggplot(species_metadata, aes(x=Level1, y=species_number, color=Level1)) +
  geom_boxplot(aes(fill = Level1, shape = Level1, colour = Level1), outlier.shape = NA) +
  geom_point(position = "jitter", aes(fill = Level1,
                                      shape = Level1,
                                      colour = Level1), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("deepskyblue2",
                              "hotpink3",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #facet_grid(. ~ `Status`,
  #           scales = "free", space = "free") +
  labs(title="Species richness by Family",
       y = "Species richness",
       x = "Level1 Diagnosis") +
  theme(#Set the title font size
    plot.title = element_text(size=10,
                              hjust = 0),
    #Set the legend title position
    legend.position = "bottom",
    #Set the legend title font size
    legend.title = element_text(size=8),
    #Define the size of the legend text and make it italic
    legend.text = element_text(size=8,
                               face = "italic"),
    #Remove the grey background from the legend
    legend.background = element_blank(),
    #Remove the box from the legend
    legend.key = element_blank(),
    #Remove the grey background
    panel.background = element_blank(),
    #Remove the plot border
    panel.border = element_blank(),
    #Remove the major plot grid lines
    panel.grid.major = element_blank(),
    #Remove the minor plot grid lines
    panel.grid.minor = element_blank(),
    #Change orientation of x axis labels
    axis.text.x = element_text(angle=0,
                               hjust=0.5,
                               vjust=0, # "justifying" text
                               size=9),
    axis.title.x = element_text(size=9),
    #Define the axis title text size
    axis.title = element_text(size=9),
    #Define the axis label text size
    axis.text.y = element_text(size=9),
    #Add back the x and y axis lines and define thickness (size), line type, and color
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"))

ggsave("SpeciesRichness-Genus-unfiltered.pdf", width = 150, height = 150,
       units = c("mm"), device = "pdf")

### SHANNON DIVERSITY ###
# Shannon Diversity stats
shannon_diversity <- vegan::diversity(abundance_table, index = "shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(metadata, shannon_diversity)
shannon_metadata$Level1 <- as.factor(shannon_metadata$Level1)

wilcox.test(shannon_diversity ~ Level1, data = shannon_metadata) # p=0.9

shannon_metadata <- cbind(metadata, shannon_diversity)

#Check the taxa level used
ggplot(shannon_metadata, aes(x=Level1, y=shannon_diversity, color=Level1)) +
  geom_boxplot(aes(fill = Level1, shape = Level1, colour = Level1), outlier.shape = NA) +
  geom_point(position = "jitter", aes(fill = Level1,
                                        shape = Level1,
                                        colour = Level1), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("deepskyblue2",
                              "hotpink3",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  labs(title="Shannon Diversity by Family with unfiltered 16S data",
       y = "Shannon Diversity Index",
       x = "Level1 Diagnosis") +
  theme(#Set the title font size
    plot.title = element_text(size=10,
                              hjust = 0),
    plot.subtitle = element_text(size = 9,
                                 hjust = 0),
    #Set the legend title position
    legend.position = "bottom",
    #Set the legend title font size
    legend.title = element_text(size=8),
    #Define the size of the legend text and make it italic
    legend.text = element_text(size=8,
                               face = "italic"),
    #Remove the grey background from the legend
    legend.background = element_blank(),
    #Remove the box from the legend
    legend.key = element_blank(),
    #Remove the grey background
    panel.background = element_blank(),
    #Remove the plot border
    panel.border = element_blank(),
    #Remove the major plot grid lines
    panel.grid.major = element_blank(),
    #Remove the minor plot grid lines
    panel.grid.minor = element_blank(),
    #Change orientation of x axis labels
    axis.text.x = element_text(angle=0,
                               hjust=0.5,
                               vjust=0, # "justifying" text
                               size=9),
    axis.title.x = element_text(size=9),
    #Define the axis title text size
    axis.title = element_text(size=9),
    #Define the axis label text size
    axis.text.y = element_text(size=9),
    #Add back the x and y axis lines and define thickness (size), line type, and color
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"))

ggsave("ShannonDiversity-Genus-unfiltered.pdf", width = 150, height = 150,
       units = c("mm"), device = "pdf")
