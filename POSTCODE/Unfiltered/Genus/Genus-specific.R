setwd("~/GitHub/Imperial-Molyneaux/POSTCODE/Unfiltered/Genus")

library(tidyverse)
library(ggplot2)
library(dplyr)

data <- read_csv("Genus-normalised-unfiltered-metadata.csv") %>%
  filter(!Diagnosis == "Mock" & !Diagnosis == "Negative Control")
data$Prevotella <- data$Prevotella + data$`[Prevotella]` 
data <- data %>%
  select(-"[Prevotella]")

abund_table <- data %>%
  select(4:ncol(data))
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100

meta_table <- data %>%
  select(`sample-id`:Diagnosis)

##### STREPTOCOCCUS ######
Streptococcus <- abund_table %>% select(Streptococcus)

df <- cbind(meta_table, Streptococcus)
dfDiagnosis <- as.factor(df$Diagnosis)

ggplot(df, aes(x = Diagnosis, y = Streptococcus)) + 
  #facet_grid(.~`Sequencing run`)+
  geom_boxplot(aes(fill = Diagnosis, colour = Diagnosis)) + 
  geom_point(position = "jitter", aes(fill = Diagnosis, colour = Diagnosis)) +
  scale_colour_manual(values=c("darkorange","darkblue", "steelblue", "orchid", "darkgreen")) + 
  scale_fill_manual(values=c("white", "white","white", "white", "white")) +
  theme_classic() + 
  labs(title="Streptococcus relative abundance",
       subtitle = "Kruskal-test p < 0.0001, after correction only sig to COVID or IPF.",
       y = "Streptococcus relative abundance (%)",
       x = "Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

ggsave("Streptococcus.pdf",
       width=300,
       height=200,
       units = c("mm"), dpi=400)

kruskal.test(Streptococcus ~ Diagnosis, data = df) #p<0.001
dunnTest(Streptococcus ~ Diagnosis, data = df, method = "bonferroni" )

##### VEILONELLA ######
Veillonella <- abund_table %>% select(Veillonella)

df <- cbind(meta_table, Veillonella)
df$Diagnosis <- as.factor(df$Diagnosis)

ggplot(df, aes(x = Diagnosis, y = Veillonella)) + 
  geom_boxplot(aes(fill = Diagnosis, colour = Diagnosis)) + 
  geom_point(position = "jitter", aes(fill = Diagnosis, colour = Diagnosis)) +
  scale_colour_manual(values=c("darkorange","darkblue","steelblue", "orchid", "darkgreen")) + 
  scale_fill_manual(values=c("white", "white","white", "white", "white")) +
  theme_classic() + 
  labs(title="Veillonella relative abundance",
       subtitle = "Kruskal-test p < 0.1, after correction p > 0.1",
       y = "Veillonella relative abundance (%)",
       x = "Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

ggsave("Veillonella.pdf",
       width=300,
       height=200,
       units = c("mm"), dpi=400)

kruskal.test(Veillonella ~ Diagnosis, data = df) #p=0.07
dunnTest(Veillonella~Diagnosis, data=df, method = "bonferroni")

#### STAPHYLOCOCCUS ####
Staphylococcus <- abund_table %>% select(Staphylococcus)

df <- cbind(meta_table, Staphylococcus)
df$Diagnosis <- as.factor(df$Diagnosis)

ggplot(df, aes(x = Diagnosis, y = Staphylococcus)) + 
  #facet_grid(.~`Sequencing run`)+
  geom_boxplot(aes(fill = Diagnosis, colour = Diagnosis)) + 
  geom_point(position = "jitter", aes(fill = Diagnosis, colour = Diagnosis)) +
  scale_colour_manual(values=c("darkorange","darkblue", "steelblue", "orchid")) + 
  scale_fill_manual(values=c("white", "white","white", "white")) +
  theme_classic() + 
  labs(title="Staphylococcus relative abundance",
       y = "Staphylococcus relative abundance (%)",
       x = "Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

ggsave("Staphylococcus.pdf",
       width=300,
       height=200,
       units = c("mm"), dpi=400)

kruskal.test(Staphylococcus ~ Diagnosis, data = df) #p=0.182
