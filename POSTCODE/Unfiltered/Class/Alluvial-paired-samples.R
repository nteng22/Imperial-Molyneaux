setwd("~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/Unfiltered/Class/")

library(dplyr)
library(tidyverse)
library(reshape2)
library(vegan)
library(ggalluvial)

path <- "~/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/Unfiltered/Class/"
list.files(path)  

data <- read_csv("Class-normalised-unfiltered-metadata.csv") %>%
  filter(!Status == "control") %>% #Only select for samples that are paired
  group_by(patient_ID) %>%
  filter(n()==2) %>%
  filter(!patient_ID == "1069") # Manually remove, has two copies of OR? 

abund_table <- data %>%
  select(index,5:ncol(data))
abund_table <- column_to_rownames(abund_table, var = "index")

abund_table <- abund_table %>%
  select(where(is.numeric)) %>%
  select(where(~sum(.) > 0))

rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # Relative abundances %

meta_table <- data %>%
  select(index:patient_ID)
meta_table$`Sample-type` <- as.factor(meta_table$`Sample-type`)
meta_table$Status <- as.factor(meta_table$Status)
meta_table$patient_ID <- as.factor(meta_table$patient_ID)

## Class analysis ##
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- dim(abund_table)[2] # This needs to be the total number of colums in the dataframe.
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others

Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

sample_data <- cbind(meta_table, top_other)
sample_data_long <- melt(sample_data, id.vars = c("index",
                                                  "Sample-type",
                                                  "Status",
                                                  "patient_ID"), variable.name = "Class")
sample_data_long

sample_data_long_summarised <- sample_data_long %>%
  group_by(`Sample-type`, Status, patient_ID, Class) %>%
  summarize(value=mean(value))

taxa_list
Palette <- c(Clostridia = "darkseagreen1",
             Bacteroidia = "darkslategray2",
             Bacilli = "hotpink3",
             Fusobacteriia = "goldenrod1",
             Actinobacteria = "palegreen3",
             Gammaproteobacteria = "plum",
             Flavobacteriia = "tan2",
             Betaproteobacteria = "darkseagreen4",
             Spirochaetes = "darkslategray4",
             Erysipelotrichi = "deeppink3",
             Others = "grey")

plot <- ggplot(data = sample_data_long_summarised,
               aes(x = `Sample-type`,
                   stratum = value,
                   alluvium = Class,
                   y = value,
                   fill = Class, 
                   label = Class)) +
  geom_alluvium(aes(fill = Class), decreasing = FALSE) +
  geom_flow(decreasing = FALSE) +
  geom_stratum(alpha = 0.5, decreasing = FALSE) +
  geom_text(stat = "alluvium",
            aes(label = Class), size= 1.75,
            decreasing = FALSE) +
  facet_wrap(patient_ID ~ . ) +
  scale_fill_manual(values = Palette) +
  theme_classic()  

plot <- plot + labs(x = "Sample type", y = "Relative abundance (%)",
                    title = "Relative abundance of Classes") + 
  theme(legend.text = element_text(size=8,
                                   face = "italic"))
plot

ggsave("Paired-samples alluvial Class.pdf", width = 800, height=500, units = c("mm"), dpi=400)
