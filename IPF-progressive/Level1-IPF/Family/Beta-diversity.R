#### Housekeeping ####
setwd("~/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF/Family//")


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
normalised_data <- as.data.frame(read_csv("Family-normalised.csv")) %>%
  filter(!Level1=="Negative"&!Level1=="Mock")
normalised_data <- subset(normalised_data, normalised_data$`sample-id` !="BRU.04039.BAL") #purely one ASV so this be GONE

# Analysis by only IPF patients
normalised_data <- normalised_data %>%
  filter(Diagnosis=="IPF" | Level1=="Controls") 

abundance_table <- normalised_data %>%
  select(5:ncol(normalised_data)) # Only select the reads
# 77 families
rowSums(abundance_table) #check if in percentages
abundance_table <- abundance_table/rowSums(abundance_table)*100

metadata <- normalised_data %>%
  select("sample-id":Diagnosis)
metadata$Diagnosis <- as.factor(metadata$Diagnosis)

adonis2(abundance_table ~ Level1,
        data = metadata, permutations = 999, 
        method = "bray")
# p = 0.001, it's good for unequal sample sizes.
  # But as there are more taxa, more tests are done, so prone to false positive...

nmds_calc = metaMDS(abundance_table, distance =  "bray", k = 2, 
                    plot = TRUE,
                    try = 100,
                    engine = "monoMDS")

# adding columns from data set to put it into context
data.scores = as.data.frame(scores(nmds_calc)$sites)

plot.data <- cbind(metadata, data.scores)

# Plot using ggplot2
ggplot(data = plot.data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Level1, shape = Level1, fill = Level1),
             size = 3) + 
  #geom_text(aes(x = NMDS1, y = NMDS2, label = `sample-id`)) +
  stat_ellipse(aes(colour = Level1)) +
  #scale_fill_manual(values = c("white", "white", "white")) + 
  #scale_colour_manual(values = c("deepskyblue3", "hotpink3", "darkgreen")) + 
  scale_shape_manual(values = c(1:16)) +
  labs(title = "NMDS by Family of unfiltered taxa") +
  theme(#Set the title font size
    plot.title = element_text(size=10),
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

ggsave("NMDS by Family- unfiltered.pdf", width = 200, height = 200, units = c("mm"), dpi = 400)
