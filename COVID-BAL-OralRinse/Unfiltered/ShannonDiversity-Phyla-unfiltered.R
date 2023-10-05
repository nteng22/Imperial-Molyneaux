#### Housekeeping ####
setwd("/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/")
setwd("~/Documents/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/Unfiltered/")

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
normalised_data <- read.csv("Phylum-normalised-unfiltered-metadata.csv")

abundance_table <- normalised_data %>%
  select(5:ncol(normalised_data)) # Only select the reads
rowSums(abundance_table) #check if in percentages
abundance_table <- abundance_table/rowSums(abundance_table)*100

metadata <- normalised_data %>%
  select(Sample.type:patient_ID)

#### Species Richness ####
# Calculate species number 
# Margin = 1, by row
species_number <- specnumber(abundance_table, MARGIN = 1)
species_number <- as.data.frame(species_number)
species_metadata <- cbind(metadata, species_number)

# Shannon Diversity stats
shannon_diversity <- diversity(abundance_table, index = "shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(metadata, shannon_diversity)
shannon_metadata$Sample.type <- as.factor(shannon_metadata$Sample.type)

# Healthy vs. Disease BAL only (no HC for oral rinse)
shannon_metadata_filter <- shannon_metadata %>%
  filter(Sample.type == "BAL") %>%
  filter(Status == "disease" |
           Status == "healthy")

wilcox.test(shannon_diversity ~ Status, data = shannon_metadata_filter)
# p < 0.819 phyla


# By type of sample (oral vs. BAL), can only do disease (paired samples)
shannon_metadata_filter <- shannon_metadata %>%
  filter(shannon_metadata$Status == "disease")

wilcox.test(shannon_diversity ~ Sample.type, data = shannon_metadata_filter)
# p = 0.10, phylum

shannon_metadata <- cbind(metadata, shannon_diversity)
shannon_metadata <- shannon_metadata %>%
  filter(!Status == "control") # remove control samples


#Check the taxa level used
ggplot(shannon_metadata, aes(x=Sample.type, y=shannon_diversity, color=Sample.type)) +
  geom_boxplot(aes(fill = Sample.type, shape = Sample.type, colour = Sample.type), outlier.shape = NA) +
  geom_line(aes(group = patient_ID),
            linetype = "solid",
            size = 0.05, alpha = 0.5,
            colour = "black") +
  geom_point(position = "identity", aes(fill = Sample.type,
                                        shape = Sample.type,
                                        colour = Sample.type,
                                        group = patient_ID), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  facet_grid(. ~ `Status`,
             scales = "free", space = "free") +
  labs(title="Shannon Diversity by Phyla",
       y = "Shannon Diversity Index",
       x = "Sample type") +
  theme(#Set the title font size
    plot.title = element_text(size=8,
                              hjust = 0.5),
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
                               size=8),
    axis.title.x = element_text(size=8),
    #Define the axis title text size
    axis.title = element_text(size=8),
    #Define the axis label text size
    axis.text.y = element_text(size=8),
    #Add back the x and y axis lines and define thickness (size), line type, and color
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"))

ggsave("ShannonDiversity-Phyla-unfiltered.pdf", width = 150, height = 150,
       units = c("mm"), device = "pdf")
