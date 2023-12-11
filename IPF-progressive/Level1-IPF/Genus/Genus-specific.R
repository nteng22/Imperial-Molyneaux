setwd("~/GitHub/Imperial-Molyneaux/IPF-progressive/Level1-IPF/Genus/")

library(tidyverse)
library(ggplot2)
library(dplyr)

data <- read_csv("Genus-normalised.csv") %>%
  filter(!Level1 == "Mock")

abund_table <- data %>%
  select(5:ncol(data))
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100

meta_table <- data %>%
  select(`sample-id`:Diagnosis)
meta_table$Level1 <- gsub("Controls", "Control", meta_table$Level1)

##### STREPTOCOCCUS ######
Streptococcus <- abund_table %>% select(Streptococcus)

df <- cbind(meta_table, Streptococcus)
df$Level1 <- as.factor(df$Level1)

ggplot(df, aes(x = Level1, y = Streptococcus)) + 
  #facet_grid(.~`Sequencing run`)+
  geom_boxplot(aes(fill = Level1, colour = Level1)) + 
  geom_point(position = "jitter", aes(fill = Level1, colour = Level1)) +
  scale_colour_manual(values=c("darkorange","darkblue", "steelblue", "orchid")) + 
  scale_fill_manual(values=c("white", "white","white", "white")) +
  theme_classic() + 
  labs(title="Streptococcus relative abundance by Level1-Fibrosis",
       y = "Streptococcus relative abundance (%)",
       x = "Level 1 Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

ggsave("Streptococcus- IPF vs HC.pdf",
       width=300,
       height=200,
       units = c("mm"), dpi=400)

IPF <- df %>% filter(Level1 == "Fibrosis")
IPF_mean <- mean(IPF$Streptococcus)

HC <- df %>% filter(Level1 == "Control")
HC_mean <- mean(HC$Streptococcus)

wilcox.test(IPF$Streptococcus, HC$Streptococcus) #0.03

##### VEILONELLA ######
Veillonella <- abund_table %>% select(Veillonella)

df <- cbind(meta_table, Veillonella)
df$Level1 <- as.factor(df$Level1)

ggplot(df, aes(x = Level1, y = Veillonella)) + 
 facet_grid(.~`Sequencing run`) +
  geom_boxplot(aes(fill = Level1, colour = Level1)) + 
  geom_point(position = "jitter", aes(fill = Level1, colour = Level1)) +
  scale_colour_manual(values=c("darkorange","darkblue", "orchid")) + 
  scale_fill_manual(values=c("white", "white","white", "white")) +
  theme_classic() + 
  labs(title="Veillonella relative abundance by Level1-Fbrosis and Sequencing run",
       y = "Veillonella relative abundance (%)",
       x = "Level 1 Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

ggsave("Veillonella- IPF vs HC by MSQ run.pdf",
       width=300,
       height=200,
       units = c("mm"), dpi=400)

IPF <- df %>% filter(Level1 == "Fibrosis")
IPF_mean <- mean(IPF$Veillonella)

HC <- df %>% filter(Level1 == "Control")
HC_mean <- mean(HC$Veillonella)

wilcox.test(IPF$Veillonella, HC$Veillonella) #0.4

#### STAPHYLOCOCCUS ####
Staphylococcus <- abund_table %>% select(Staphylococcus)

df <- cbind(meta_table, Staphylococcus)
df$Level1 <- as.factor(df$Level1)

ggplot(df, aes(x = Level1, y = Staphylococcus)) + 
  #facet_grid(.~`Sequencing run`)+
  geom_boxplot(aes(fill = Level1, colour = Level1)) + 
  geom_point(position = "jitter", aes(fill = Level1, colour = Level1)) +
  scale_colour_manual(values=c("darkorange","darkblue", "orchid")) + 
  scale_fill_manual(values=c("white", "white","white", "white")) +
  theme_classic() + 
  labs(title="Staphylococcus relative abundance by Level1-Fibrosis",
       y = "Staphylococcus relative abundance (%)",
       x = "Level 1 Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

ggsave("Staphylococcus- IPF vs HC.pdf",
       width=300,
       height=200,
       units = c("mm"), dpi=400)

IPF <- df %>% filter(Level1 == "Fibrosis")
IPF_mean <- mean(IPF$Staphylococcus)

HC <- df %>% filter(Level1 == "Control")
HC_mean <- mean(HC$Staphylococcus)

wilcox.test(IPF$Staphylococcus, HC$Staphylococcus) #0.4
