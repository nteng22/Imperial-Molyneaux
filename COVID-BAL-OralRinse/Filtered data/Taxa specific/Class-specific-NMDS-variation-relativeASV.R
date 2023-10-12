setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/Filtered data/Taxa specific/")

library("ggplot2")
library("tidyverse")
library("dplyr")
library("reshape2")

taxa = "Class"
data <- read.csv("../../Class-normalised-metadata.csv") %>%
  filter(!Status == "control") # remove controls


meta_table <- data %>%
  select(index:Status)
# Based off of taxa analysis (bionenv): Coriobacteriia, clostridia and Bacilli
  # were the drivers of variance, so extract this from the dataset
abund_table <- data %>%
  select("Coriobacteriia", "Clostridia", "Bacilli")

data_subset <- cbind(meta_table, abund_table)

data_subset_long <- melt(data_subset, id.vars = c("index",
                                                  "Sample.type",
                                                  "Status"), variable.name = "Class")
data_subset_long

plot <- ggplot(data_subset_long, 
              # %>% filter(Class == "Coriobacteriia" ),
               aes(x = `Sample.type`, y = value)) +
  geom_boxplot(aes(fill = Status, shape = Status, colour = Status)) + #, outlier.shape = NA) +
  geom_point(position=position_dodge(width = 0.75), aes(x = `Sample.type`, y = value,
                                             fill = Status, shape = Status, colour = Status)) +
  facet_wrap(. ~ Class) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("red",
                              "blue",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="Relative abundances of potential Phyla driving NMDS variation",
       subtitle = "Proportional ASVs",
       y = "Relative abundance (%)",
       x = "Sample Type") +
  theme(#Set the aspect ratio of the plot
    aspect.ratio = 1.5)

# Save by pdf
plot

pdf(file = "Relative ASVs Coriobacteriia driving NMDS variation.pdf", width = 11, height = 8)

print(plot)
dev.off()
