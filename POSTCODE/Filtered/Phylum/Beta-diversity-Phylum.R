#### Housekeeping ####
setwd("~/Documents/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/Unfiltered/Phylum/")

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
normalised_data <- read.csv("Phylum-normalised-unfiltered-metadata.csv") %>%
  filter(!Status == "control") #remove controls

abundance_table <- normalised_data %>%
  select(5:ncol(normalised_data)) # Only select the reads
rowSums(abundance_table) #check if in percentages
abundance_table <- abundance_table/rowSums(abundance_table)*100

metadata <- normalised_data %>%
  select(Sample.type:patient_ID) 
metadata$Sample.type <- as.factor(metadata$Sample.type)
metadata$Status <- as.factor(metadata$Status)
metadata$patient_ID <- as.factor(metadata$patient_ID)

adonis2(abundance_table ~ Status,
        data = metadata, permutations = 999, 
        method = "bray")

adonis2(abundance_table ~ Sample.type, 
        data = metadata, permutations = 999,
        method = "bray")

# Does status influence sample type
adonis2(abundance_table ~ Status + Sample.type + Status*Sample.type,
        data = metadata, permutations = 999,
        method = "bray")

nmds_calc = metaMDS(abundance_table, distance =  "bray", k = 2, 
                    plot = TRUE,
                    try = 100,
                    engine = "monoMDS")

# adding columns from data set to put it into context
data.scores = as.data.frame(scores(nmds_calc)$sites)

plot.data <- cbind(metadata, data.scores)

# Plot using ggplot2
ggplot(data = plot.data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Sample.type, shape = Status, fill = Sample.type),
             size = 3) + 
  stat_ellipse(aes(colour = Sample.type)) +
  # scale_"" is used to design the plot
  scale_fill_manual(values = c("white", "white", "white")) + 
  scale_colour_manual(values = c("deepskyblue3", "hotpink3", "darkgreen")) + 
  scale_shape_manual(values = c(1, 16)) +
  labs(title = "NMDS by Phylum of unfiltered taxa",
       subtitle = "PERMANOVA . ~ Sample-type p=0.02") +
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size=8),
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
    axis.text.x = element_text(angle=90,
                               hjust=0,
                               vjust=0.5, # "justifying" text
                               size=8),
    axis.title.x = element_text(size=8),
    #Define the axis title text size
    axis.title = element_text(size=8),
    #Define the axis label text size
    axis.text.y = element_text(size=8),
    #Add back the x and y axis lines and define thickness (size), line type, and color
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = 
                               "black"),
    #Set the aspect ratio of the plot
    aspect.ratio = 1)

ggsave("NMDS by Phylum- unfiltered.pdf", width = 200, height = 200, units = c("mm"), dpi = 400)
