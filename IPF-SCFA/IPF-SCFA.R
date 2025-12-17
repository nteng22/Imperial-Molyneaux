setwd("~/Documents/GitHub/Imperial-Molyneaux/IPF-SCFA/")
#renv::activate()
renv::restore()

setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT001_Microbiota-analysis/Propionate/")
path <- getwd()

library(ggpicrust2)
library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
#BiocManager::install("microbiome")
library("microbiome")
library(pacman)
#if(!require(pacman))install.packages("pacman")
pacman::p_load('dplyr', 'tidyr', 'gapminder',
               'ggplot2',  'ggalt',
               'forcats', 'R.utils', 'png', 
               'grid', 'ggpubr', 'scales',
               'bbplot')
library(Gmisc)
library("FSA")
library(ggbiplot)
library(readxl)
library(circlize)
library(svglite)

set.seed(2024)

Palette_color = c("#BC3C29FF", "#0072B5FF")
Palette_fill = c("#BC3C29FF", "#0072B5FF")

Palette_10 = c("#003f5a", # deep navy (blue-teal)
               "#0072B5", # blue
               "#6F99AD", # slate blue
               "#6ea6a4", # muted turquoise
               "#20854E", # green-teal
               "#8ab184", # sage
               "#C4644A", # muted pink
               "#BC3C29", # red-brown
               "#de6600", # orange
               "#E18727", # gold
               "#FFDC91", # light gold
               "#fec682" # peach
)

Palette_6 = c("#003f5a", # deep navy (blue-teal)
              "#6F99AD", # slate blue
              "#20854E", # green-teal
              "#BC3C29", #red-brown
              "#de6600", # orange
              "#FFDC91"  # light gold
)

Palette_diverging1 <- c("#003f5a", "white", "#BC3C29")
Palette_diverging1 <- colorRampPalette(Palette_diverging1)
Palette_diverging1 <- Palette_diverging1(50)

# Load the qiime2 artefacts into the environment
physeq<-qza_to_phyloseq(
  features="Original-files/filtered-metadata-table.qza",
  tree="Original-files/merged-midrooted-tree.qza",
  taxonomy= "Original-files/taxonomy-merged-final.qza")
physeq 
metadata = read_csv("Original-files/Propionate_metadata.csv")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 762 taxa and 255 samples ]
#tax_table()   Taxonomy Table:    [ 762 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 762 tips and 757 internal nodes ]

metadata = metadata %>%
  mutate(Sample = case_when(Sample == "ILDCON.1022" ~ "BRU.1022",
                            Sample == "ILDCON.1066" ~ "BRU.1066",
                            Sample == "ILDCON.1017" ~ "BRU.1017",
                            Sample == "ILDCON.1032" ~ "BRU.1032",
                            Sample == "ILDCON.1031" ~ "BRU.1031",
                            Sample == "ILDCON.1066" ~ "BRU.1066",
                            Sample == "ILDCON.1040" ~ "BRU.1040",
                            Sample == "ILDCON.1037" ~ "BRU.1037",
                            Sample == "ILDCON.1039" ~ "BRU.1039",
                            Sample == "ILDCON.1036" ~ "BRU.1036",
                            Sample == "ILDCON.1042" ~ "BRU.1042",
                            Sample == "ILDCON.1000" ~ "BRU.1000",
                            Sample == "ILDCON.1028" ~ "BRU.1028",
                            Sample == "ILDCON.1018" ~ "BRU.1018",
                            Sample == "ILDCON.1006" ~ "BRU.1006",
                            .default = Sample))
metadata <- column_to_rownames(metadata, var = "Sample")
mapping=metadata 

# Demographics
demo_table = metadata %>%
  filter(!is.na(Age)) %>%
  group_by(Diagnosis) %>%
  summarise(mean_age = round(mean(Age)),
            sd_age = sd(Age),
            sex_M = sum(Sex == "M"),
            sex_M_prop = round(sex_M/n()*100),
            FVC_pp = round(mean(FVCpp, na.rm = T)),
            FVC_pp_sd = round(sd(FVCpp, na.rm = T)),
            DLCO_pp = round(mean(`DLCO % pred`, na.rm = T)),
            DLCO_pp_sd = round(sd(`DLCO % pred`, na.rm = T)),
            N_Never_Smoked = sum(Smoking_history == "Never-smoker", na.rm = TRUE),
            N_never_percent =  round(N_Never_Smoked/n()*100),
            N_Ex_Smoker = sum(Smoking_history == "Ex-smoker", na.rm = TRUE),
            N_ex_percent = round(N_Ex_Smoker/n()*100),
            N_Current_Smoker = sum(Smoking_history == "Current-smoker", na.rm = TRUE),
            N_current_percent = round(N_Current_Smoker/n()*100),
            N_Nodata_smoker = sum(is.na(Smoking_history)),
            N_nodata_percent = round(N_Nodata_smoker/n()*100))

metadata_table1 = metadata %>%
  filter(!Diagnosis == "Negative Controls") %>%
  mutate(Sex = as.factor(Sex),
         Smoking_history = as.factor(Smoking_history))
wilcox.test(Age ~ Diagnosis, metadata_table1)
wilcox.test(FVCpp ~ Diagnosis, metadata_table1)
wilcox.test(`DLCO % pred` ~ Diagnosis, metadata_table1)
fisher.test(metadata$Smoking_history, metadata$Diagnosis)
chisq.test(metadata$Sex, metadata$Diagnosis)

sample_data(physeq)<-mapping # metadata dataframe becomes the sample_data slot of physeq object
#view(sample_data(physeq))  
gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq)))
setdiff(rownames(metadata), sample_names(physeq)) 
  # Corresponds to the number of samples not sequenced. 

# Contamination check
ps<-physeq
summarize_phyloseq(ps)
# "1] Min. number of reads = 93"
# "2] Max. number of reads = 1339062"
# "3] Total number of reads = 9414683"
# "4] Average number of reads = 36920.3254901961""
# "5] Median number of reads = 30595"
# "6] Any OTU sum to 1 or less? YES"
# "7] Sparsity = 0.95119654160877"
# "8] Number of singletons = 1"
# "9] Percent of OTUs that are singletons \n 
# "10] Number of sample variables are: 20"

## Plot to look at library sizes per sample 
# Separated by Diagnosis and sample type i.e., negatives vs. true samples. 
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Diagnosis)) + geom_point() +
  ylim(lower = 0, upper = 80000) # one sample is removed, BRU.1000. Is clearly an outlier due to the ridiculous librarysize (1339062).

#### Removal of contaminants qPCR ####
ps_qPCR <- subset_samples(ps, !is.na(sample_data(ps)$Burden))
sample_data(ps_qPCR)$is.neg <- sample_data(ps_qPCR)$Diagnosis=="Negative Control" #Only considered reagent control as negative
view(sample_data(ps_qPCR))
contamdf.freq.qPCR<- isContaminant(ps_qPCR, method="combined", conc="Burden", neg="is.neg", threshold=0.1) 
# Gets data frame with a list of potential contaminants

table(contamdf.freq.qPCR$contaminant)#FALSE=757, TRUE=5
head(which(contamdf.freq.qPCR$contaminant)) # not the highest abundant taxa
# [1] 253 391 476 654 659

plot_frequency(ps_qPCR, taxa_names(ps)[sample(which(contamdf.freq.qPCR$contaminant),5)], conc="Burden") +
  xlab("DNA Concentration (gene/mL)")

# Presence/Absence plot- qPCR
ps.pa <- transform_sample_counts(ps_qPCR, function(abund) 1*(abund>0))
# Transforms data to absence or presence 
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Diagnosis =="Negative Control", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$Diagnosis == "Negative Control", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa.qPCR <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                          contaminant=contamdf.freq.qPCR$contaminant)
ggplot(data=df.pa.qPCR, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") 

tax <- as(tax_table(ps), "matrix")
contaminants<- tax[which(contamdf.freq.qPCR$contaminant),]
write.table(contaminants,file="Original-files/Contaminants-list-qPCR.txt", col.names=NA, row.names=T,sep="\t")

#### qPCR contaminants ####
ps.contam <- prune_taxa(contamdf.freq.qPCR$contaminant, ps) # only subset contaminants from the dataset
plot_bar(ps.contam)

# Identify ASVs above 1000 reads in the list of contaminants. These are big influencers. 
filter <- phyloseq::genefilter_sample(ps.contam, filterfun_sample(function(x) x >= 1000))
ps.contam.1k <- prune_taxa(filter, ps.contam)
otu_table<-as.data.frame(ps.contam.1k@otu_table)
tax_table<-as.data.frame(ps.contam.1k@tax_table) # only lactobacillus 

library(tibble)
library(dplyr)

sum(taxa_sums(ps) == 0) # how many taxa aren't present in ANY samples, 
summarize_phyloseq(ps)
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 762 taxa and 255 samples ]
#sample_data() Sample Data:       [ 255 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 762 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 762 tips and 757 internal nodes ]

#If happy to remove all these then use following command using qPCR data
ps.noncontam <- prune_taxa(!contamdf.freq.qPCR$contaminant, ps)
ps.noncontam
summarize_phyloseq(ps.noncontam)

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 757 taxa and 255 samples ]
#sample_data() Sample Data:       [ 255 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 757 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 757 tips and 752 internal nodes ]

Controls_IPF <-subset_samples(ps.noncontam, !Diagnosis=="Negative Control")
sum(taxa_sums((Controls_IPF)) == 0) #34

summarize_phyloseq(Controls_IPF)

Controls_IPF

filterPhyla = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria") 
ps1 = subset_taxa(Controls_IPF, !(!Phylum %in% filterPhyla))
ps1

#Create final filtered object
Controls_IPF_decontaminated <-ps1

#Remove contaminants manually (i.e. Sphingomonas) and after you have removed them create a final phyloseq object (i.e. CHP_decontaminated_final)

badTaxa = c("bdb9e5f6c41c701ec4b3bc754f1fde53", "6046915cd4ef8792461b808c4ae8867c", "e84f6f8c8e575aa1495ef93bb3efbca5", "24c51df570756f1f5da7d9980208ee7c", "5165f65d3a69208c18841bda88097e51", "248470c3b996979c8d6409548617f915", "bc2e91c177dcf042f2ecde858c2e64ed", "e3d7402dfa7a45343b415c679366e347", "e3b0ad72978c9880998000ef253cb87c", "bba6da0e81f92e0b2af4e957cd69f42e", "fda3f8be3cbb25ffb216c5cfb50c61ea", "e85e34c0ba855c2e668db26708ee0d4f", "8bd8a20601c0703a841ecd04fc98a423", "b4192b3f531d68bf43083142eb473463", "8f9ca42c10f62423906999a969413bde", "12758c31900e156a57c226c81e4365cd", "67eeddd27a560751dff4bc2bd12e4acc", "894be7ef3d094e9930c61db2c4760ab2", "4f16e832bf84501943b0ba315d5b2acc", "9c1454b23d7b9e86ad391e0b238ee420", "4ef3d05269c762db631f033432d48b23", "7984f1cd22e30fdf65fd3e4c5038e853", "2b2ef93955a0e7fc9940b16759cf6c2a", "80a5cb9001fa1356ec097549337a0a64", "2974d9feee2aa60a8680500782e838c0","9f8b16f780589d38ccbb4a92b6cc258b", "ec62f5ab04d9c4ebc9aa3850f7e776c5", "066a26869b150e119af1ba42c9560827", "e78d9715d827eba81db3c87e840b3a16", "873532098886a3068e333839cff270e3", "6f05bbb7b495e6ffda5b2dacaf47c75c", "f0d3ad42c3754a9af2c18ab6f102bb28", "34b32933f7ccf504f3d591952229610c", "cc6669d31bc57687068868bc38e8e168", "654878c72a7532eac049ef9a09758f62", "74c3b8fe5db4ddf4eca09d816aedd8e2", "b1ebdc7dc82cf510768b062efab89fc9", "2e9f64ce9ec727ee1fde65c75aa71db7", "f65857c80ee4b6c9d502aaaed0ae786f", "78005887af64a050cbdfbfe834f60e1b", "14eeed0c81e495bdd16555cb4bc5d4c7", "871fcc3c84b0af8daad8418bebe9534c", "dd4e2672ac27cfe7a3d6b4fbf50c7ca5", "0b12835fa498ceaf800913bca8764435", "b9faf4c2c3f8a00a325acccd1a46ab93", "d524e847be1a6f4918baa05f5214f0a7", "09b3ffde2779e3c85da11372ef30cd1d", "98cd57089b10be0f057e407e11f06f71", "1006d42c9984f242af419f88f8f30b61", "147d6ecf30b0296147481cf4e4b58c28", "27c8a86b35998a0bc6cc4fc407a715df", "f3389b554fa53e0c16a261d364470eae", "2897f7c0c302e175164c6412f262d5b4")
goodTaxa <- setdiff(taxa_names(Controls_IPF_decontaminated), badTaxa)
Controls_IPF_decontaminated_final <- prune_taxa(goodTaxa, Controls_IPF_decontaminated)

#View feature table
Controls_IPF_decontaminated_final
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 738 taxa and 239 samples ]
#sample_data() Sample Data:       [ 239 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 738 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 738 tips and 733 internal nodes ]

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(Controls_IPF_decontaminated_final)
y = taxa_sums(ps.noncontam)

keepTaxa = which((x / sum(x)) > minTotRelAbun)
keepTaxa.y = which((y / sum(y)) > minTotRelAbun)

prunedSet = prune_taxa(names(keepTaxa), Controls_IPF_decontaminated_final)
prunedSet.all = prune_taxa(names(keepTaxa.y), ps.noncontam)

prunedSet_OTU <- as(otu_table(Controls_IPF_decontaminated_final), "matrix")
prunedSet_OTU <- as.data.frame(prunedSet_OTU)
write.csv(prunedSet_OTU, "Original-files/OTUdf_WGCNA.csv") # need count data for WGCNA

prunedSet_taxa <- as(tax_table(Controls_IPF_decontaminated_final), "matrix")
prunedSet_taxa <- as.data.frame(prunedSet_taxa)
write.csv(prunedSet_taxa, "Original-files/Taxa-df_WGCNA.csv")

#Final object with no contaminants (no negative controls) filtered at 0.005% 
prunedSet
prunedSet.all # includes negatives.
summarize_phyloseq(prunedSet)

# Normalize reads
normalizeSample = function(x) {
  x/sum(x)*100
}

Controls_IPF_relative = transformSampleCounts(prunedSet, normalizeSample)
prunedSet.all_relative = transformSampleCounts(prunedSet.all, normalizeSample)

#Extract OTU table and coerce to data.frame
OTU1 = as(otu_table(Controls_IPF_relative), "matrix")
OTU2 = as(otu_table(prunedSet.all_relative), "matrix")

#Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf.all = as.data.frame(OTU2)

#Export as csv file 
write.csv(OTUdf, "Original-files/OTUdf_qPCR.csv")

#Export tree 
tree = phy_tree(Controls_IPF_relative)
ape::write.tree(tree, "Original-files/tree_qPCR.nwk")

#Extract TAX table and coerce to data.frame
TAXdf = as(tax_table(Controls_IPF_relative), "matrix")

#Coerce to data.frame
TAXdf = as.data.frame(TAXdf)

#Export as csv file 
write.csv(TAXdf, "Original-files/tax_table_qPCR.csv")

Controls_Phylum <- aggregate_taxa(Controls_IPF_relative, 'Phylum') #7 phyla most likely a contaminant
Controls_Phylum.all <- aggregate_taxa(prunedSet.all_relative, 'Phylum') #7 phyla most likely a contaminant

otu_table<-as.data.frame(Controls_Phylum@otu_table)
otu_table.all <- as.data.frame(Controls_Phylum.all@otu_table)

write.table(otu_table,file="Phylum/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")
write.table(otu_table.all,file="Phylum/Phylum-all-samples-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Class <- aggregate_taxa(Controls_IPF_relative, 'Class')
otu_table<-as.data.frame(Controls_Class@otu_table)
write.table(otu_table,file="Class/Class-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Family <- aggregate_taxa(Controls_IPF_relative, 'Family')
otu_table<-as.data.frame(Controls_Family@otu_table)
write.table(otu_table,file="Family/Family-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Genus <- aggregate_taxa(Controls_IPF_relative, 'Genus')
Controls_Genus.all <- aggregate_taxa(prunedSet.all_relative, 'Genus')

otu_table<-as.data.frame(Controls_Genus@otu_table)
otu_table.all <- as.data.frame(Controls_Genus.all@otu_table)

write.table(otu_table,file="Genus/Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")
write.table(otu_table,file="Genus/Genus-all-samples-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

#### Joining reads with clinical metadata ####
## Have joined clinical metadata to the taxonomic reads. 
# every row is a sample and every column is a unique taxa member e.g., genus
metadata_clean <- rownames_to_column(metadata, "sample-id")

# Negatives are removed
reads_genus <- read_tsv("Genus/Genus-relative-abundance.txt")
reads_genus <- column_to_rownames(reads_genus, var = "...1")
reads_genus <- as.data.frame(t(reads_genus))
reads_genus <- rownames_to_column(reads_genus, var = "sample-id")
reads_genus_metadata <- right_join(metadata_clean, reads_genus)
write_csv(reads_genus_metadata, "Genus/Genus-normalised-metadata.csv")

reads_family <- read_tsv("family/family-relative-abundance.txt")
reads_family <- column_to_rownames(reads_family, var = "...1")
reads_family <- as.data.frame(t(reads_family))
reads_family <- rownames_to_column(reads_family, var = "sample-id")
reads_family_metadata <- right_join(metadata_clean, reads_family)
write_csv(reads_family_metadata, "family/family-normalised-metadata.csv")

reads_class <- read_tsv("class/class-relative-abundance.txt")
reads_class <- column_to_rownames(reads_class, var = "...1")
reads_class <- as.data.frame(t(reads_class))
reads_class <- rownames_to_column(reads_class, var = "sample-id")
reads_class_metadata <- right_join(metadata_clean, reads_class)
write_csv(reads_class_metadata, "class/class-normalised-metadata.csv")

reads_phylum <- read_tsv("phylum/phylum-relative-abundance.txt")
reads_phylum <- column_to_rownames(reads_phylum, var = "...1")
reads_phylum <- as.data.frame(t(reads_phylum))
reads_phylum <- rownames_to_column(reads_phylum, var = "sample-id")
reads_phylum_metadata <- right_join(metadata_clean, reads_phylum)
write_csv(reads_phylum_metadata, "phylum/phylum-normalised-metadata.csv")

#### Diversity indices ####
Controls_IPF_decontaminated_final # the decontaminated ps object but NOT filtered out low taxa
prunedSet # decontaminated AND filtered low taxa

p = plot_richness(prunedSet, x="Diagnosis", color="Diagnosis", measures=c("Observed","Shannon", "Chao1"))

alpha_diversity <- p$data
alpha_diversity_shannon <- alpha_diversity %>% filter(variable=="Shannon")
wilcox.test(value ~ Diagnosis, alpha_diversity_shannon) #p=0.13

p + geom_boxplot(data = p$data, aes(x = Diagnosis, y = value, color = NULL), alpha = 0.1) +
  scale_colour_manual(values=Palette_color) + 
  theme_classic()+ theme(axis.text.x=element_blank())

#### Beta-plot ASV ####
x = taxa_sums(Controls_IPF_decontaminated_final)
# Keep taxa seen at least twice in more than 1% of samples.
filteredset = filter_taxa(Controls_IPF_decontaminated_final, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
filteredset

minTotRelAbun = 0.00005
x = taxa_sums(filteredset)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet_beta = prune_taxa(names(keepTaxa), filteredset)

normalizeSample = function(x) {
  x/sum(x)*100}
tax <- as(tax_table(Controls_IPF_decontaminated_final), "matrix")

Controls_relative_beta = transformSampleCounts(prunedSet_beta, normalizeSample)
relative_beta_genus <- aggregate_taxa(Controls_relative_beta, "Genus")
otu_relative_beta_genus <- as.data.frame(relative_beta_genus@otu_table)

relative_beta_genus <- as.data.frame(t(otu_relative_beta_genus))
relative_beta_genus <- rownames_to_column(relative_beta_genus, var = "sample-id")
relative_beta_genus <- right_join(metadata_clean, relative_beta_genus)
write_csv(relative_beta_genus, "Genus/Beta_Genus-normalised-metadata.csv")

relative_beta_phylum <- aggregate_taxa(Controls_relative_beta, "Phylum")
otu_relative_beta_phylum <- as.data.frame(relative_beta_phylum@otu_table)

relative_beta_phylum <- as.data.frame(t(otu_relative_beta_phylum))
relative_beta_phylum <- rownames_to_column(relative_beta_phylum, var = "sample-id")
relative_beta_phylum <- right_join(metadata_clean, relative_beta_phylum)
write_csv(relative_beta_phylum, "Phylum/Beta_phylum-normalised-metadata.csv")

ordu = ordinate(filteredset, "PCoA", "unifrac", weighted=TRUE)
#Temporary fix to colour bug is to add a "dummy variable" in sample_data:
#sample_data(physeq)[ , 2] <- sample_data(physeq)[ ,1]
p = plot_ordination(filteredset, ordu, color="Diagnosis") + geom_point(size = 3) +
  ggtitle("PCoA on weighted-UniFrac distance") + (scale_colour_brewer(type="qual", 
                                                                      palette="Set1"))
print(p)
myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_colour_manual(values=Palette_color)

#### Top ten stacked bar plot ####
genus = read_csv("Genus/Genus-normalised-metadata.csv")
abund_genus <- genus %>%
  dplyr::select(1,Actinomyces:ncol(genus))
abund_genus <- column_to_rownames(abund_genus, var="sample-id")
rowSums(abund_genus)

meta_table <- genus %>%
  dplyr::select(`sample-id`, Diagnosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_genus[,order(colSums(abund_genus),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_genus[,order(colSums(abund_genus),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(genus)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

ten_genus <- cbind(meta_table, top_other)

# OPTIONAL: Fix the order of sample IDs in the order of genus proportion.
ten_genus$`sample-id` <- factor(ten_genus$`sample-id`,
                                levels=unique(ten_genus$`sample-id`))

ten_genus_long <- melt(ten_genus, id.vars = c("sample-id",
                                                  "Diagnosis"), variable.name = "Genus")
ten_genus_long
taxa_list

ggplot(ten_genus_long[ten_genus_long$Diagnosis=="IPF",],
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette_10) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  facet_grid(.~Diagnosis,
             scales = "free_x",
             drop=TRUE)+
  labs(title="Ten most abundant genera in IPF and Control subjects",
       x = "Sample ID",
       y = "Relative abundance (%)")+ theme_classic()+
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text = element_text(size=8,
                               face = "italic"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.x = element_text(size=8),
    axis.title = element_text(size=8),
    axis.text.y = element_text(size=8),
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    aspect.ratio = 0.5
  )

ggsave("Genus/Stackedbarplot-outliers.tiff", width=20, height=18)

# BRU.03754, BRU.03715 is mainly "Others"
# BRU.05380, BRU.03720, BRU.04110 "Others" and "Unknown"
# BRU.04338 Streptococcus
# BRU.03782 Staphylococcus
# BRU.05244, BRU.04052 Staphylooccus and "Others"
# BRU.04340 Streptococcus and "Others"

#### Heatmap to see if "outliers" cluster together ####
set.seed(123)
top <- abund_genus[,order(colSums(abund_genus),decreasing=TRUE)]
N <- 30
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_genus[,order(colSums(abund_genus),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(genus)[2]
taxa_list2 <- colnames(other)[31:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
rowSums(top_other)

sample_data <- cbind(meta_table, top_other)

# Heatmap will be made using a distance matrix
data_dist_rows <- vegdist(top_other, method = "bray")
row_clustering <- hclust(data_dist_rows, "average")
data_dist_columns <- vegdist(t(top_other), method = "bray")
col_clustering <- hclust(data_dist_columns, "average")

colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)

# Add annotation to heatmap AFTER distance matrix has been calculated
# the metadata does not influence the distance matrix. 
annot_df1 <- data.frame(Diagnosis = meta_table$Diagnosis)

col1 = list(Diagnosis = c("Controls" = "#BC3C29FF",
                          "IPF" = "#0072B5FF"))
library(ComplexHeatmap)
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Diagnosis", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

heatmap <- Heatmap(as.matrix(top_other), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   #cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 8), # Row name font size
                   column_names_gp = gpar(fontsize = 10, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 10), # Column title font size
                   row_title = "Samples", # Set row title
                   row_labels = rownames(top_other),
                   row_names_side = c("left"),
                   row_title_gp = gpar(fontsize = 10), # Set row title font size
                   column_title = "", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 10), # Set legend title font size
                                               labels_gp = gpar(fontsize = 10))) # Set legend label font size

p <- heatmap + sidebar_annotation1
p

svglite("Genus/Heatmap-outliers.svg", width = 8, height = 12)
draw(p)
dev.off()

sample_outliers = c("BRU.03782", "BRU.05380", "BRU.04110", "BRU.03720", "BRU.03754",
                    "BRU.04052")

##### Permanova #####
phylum <- read_csv("Phylum/phylum-normalised-metadata.csv")
abund_phylum <- phylum %>% select(1,Actinobacteria, Proteobacteria, Bacteroidetes,
                                  Firmicutes, Fusobacteria) %>%
  column_to_rownames("sample-id")
# need to drop missing clinical data for permanova
meta_table_restricted <- phylum %>%
  dplyr::select(`sample-id`, Diagnosis, Age, Sex, Smoking_history, Burden,
         `FVCpp`, `DLCO % pred`) %>%
  drop_na() # lose 30 samples
abund_phylum_restricted <- abund_phylum %>%
  filter(rownames(abund_phylum) %in% meta_table_restricted$`sample-id`)
phylum_restricted <- cbind(meta_table_restricted, abund_phylum_restricted)
colnames(meta_table_restricted) = c("sample-id", "Diagnosis", "Age", 
                         "Sex", "Smoking_history", "Burden", "FVCpp", "DLCOpp")
adonis2(abund_phylum_restricted ~ Diagnosis+Sex+Smoking_history+DLCOpp+FVCpp, permutations = 999,
        data = meta_table_restricted, method="euclidean", by = "margin")
 
#### GCMS ####
metadata_GCMS <- read_csv("Original-files/Propionate_metadata.csv")
metadata_GCMS <- metadata_GCMS %>%
  filter(!is.na(`Acetic acid`)) %>%
  mutate(Diagnosis = as.factor(Diagnosis))

metadata_GCMS %>%
  group_by(Diagnosis) %>%
  summarise(median_propionate = median(Propionate))

##### Boxplots #####
# I don't know what the detectable range is for these.
wilcox.test(`Acetic acid` ~ Diagnosis, metadata_GCMS)

ggplot(metadata_GCMS, aes(x = Diagnosis, y = `Acetic acid`, colour = Diagnosis)) + 
  geom_boxplot(outliers = F) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 2.5, position="dodge", aes(fill=Diagnosis)) +
  scale_fill_manual(values=Palette_fill) + 
  theme_classic() + 
  scale_color_manual(values=Palette_color) +
  annotate("text", x = 1.5, y=168, label= "***", size = 8) +
  annotate("segment", x=1, xend=2, y=165, yend=165)+
  labs(y = "[Metabolite in BAL] (µM)",
       title = "Acetate") + ylim(c(0,200)) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5), family = "bold")

ggsave("Genus/Acetate-boxplot.svg", width = 4, height = 6)

wilcox.test(`Butyric acid` ~ Diagnosis, metadata_GCMS)
ggplot(metadata_GCMS, aes(x = Diagnosis, y = `Butyric acid`, colour = Diagnosis)) + 
  geom_boxplot(outliers = F) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 0.8, position="dodge", aes(fill=Diagnosis)) +
  scale_fill_manual(values=Palette_fill) + 
  theme_classic() + 
  scale_color_manual(values=Palette_color) +
  annotate("text", x = 1.5, y=51, label= "**", size = 8) +
  annotate("segment", x=1, xend=2, y=50, yend=50)+
  labs(y = "[Metabolite in BAL] (µM)",
       title = "Butyrate") + ylim(c(0,60)) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5), family = "bold")

ggsave("Genus/Butyrate-boxplot.svg", width = 4, height = 6)

wilcox.test(`Lactic acid` ~ Diagnosis, metadata_GCMS)
ggplot(metadata_GCMS, aes(x = Diagnosis, y = `Lactic acid`, colour = Diagnosis)) + 
  geom_boxplot(outliers = F, aes(colour = Diagnosis)) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 0.3, position="dodge", aes(fill=Diagnosis)) +
  scale_fill_manual(values=Palette_fill) + 
  theme_classic() + 
  scale_color_manual(values=Palette_color) + 
  annotate("text", x = 1.5, y=20.5, label= "**", size = 8) +
  annotate("segment", x=1, xend=2, y=20, yend=20)+
  labs(y = "[Metabolite in BAL] (µM)",
       title = "Lactate") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5), family = "bold")
ggsave("Genus/Lactic-boxplot.svg", width = 4, height = 6)

wilcox.test(`Propionic acid` ~ Diagnosis, metadata_GCMS)
ggplot(metadata_GCMS, aes(x = Diagnosis, y = `Propionic acid`, colour = Diagnosis)) + 
  geom_boxplot(outliers = F) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 2.5, position="dodge", aes(fill=Diagnosis)) +
  scale_fill_manual(values=Palette_fill) + 
  theme_classic() + 
  scale_color_manual(values=Palette_color) + ylim(c(0,200))+
  annotate("text", x = 1.5, y=150, label= "**", size = 8) +
  annotate("segment", x=1, xend=2, y=148, yend=148)+
  labs(y = "[Metabolite in BAL] (µM)",
       title="Propionate") +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5), family = "bold")

ggsave("Genus/Propionate-boxplot.svg", width = 4, height = 6)

metadata_GCMS = metadata_GCMS %>%
  rename("Acetate" = `Acetic acid`,
         "Lactate" = `Lactic acid`,
         "Butyrate" = `Butyric acid`,
         "Propionate" = `Propionic acid`)

set.seed(2024)

##### PCA ####
PCA <- prcomp(metadata_GCMS[c("Acetate", "Propionate", "Butyrate", "Lactate")],
  scale. = F)

ggbiplot(PCA, obs.scale = 0, var.scale = 0.5,
         groups=metadata_GCMS$Diagnosis, ellipse = T, circle=F,
         varname.size = 5) +
  scale_color_discrete(name="")+
  theme(legend.direction = 'horizontal', legend.position = 'top',
        legend.text = element_text(size=16),
        axis.text = element_text(size=16),
        axis.title = element_text(size=16)) + theme_classic(base_size = 10) + 
  scale_color_manual(values=Palette_color)

ggsave("Genus/GCMS-PCA.svg", width = 8, height=8)

#### Correlation heatmap ####
abund_genus = genus %>%
  select(1,Actinomyces:ncol(genus)) %>%
  column_to_rownames("sample-id")

meta_GCMS_PFT_cor_heatmap = metadata_GCMS %>%
  filter(Sample %in% rownames(abund_genus)) %>%
  mutate(Diagnosis = as.character(Diagnosis))

meta_GCMS_PFT_cor_heatmap = metadata_GCMS %>%
  filter(Sample %in% rownames(abund_genus))

abund_genus_restricted <- abund_genus %>%
  filter(rownames(abund_genus) %in% meta_GCMS_PFT_cor_heatmap$Sample)

labels = colnames(meta_GCMS_PFT_cor_heatmap)

x<-abund_genus_restricted[,order(colSums(abund_genus_restricted),decreasing=TRUE)]
#Extract list of top N Taxa
N<-11
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)
x<-data.frame(x[,colnames(x) %in% taxa_list])
y<-meta_GCMS_PFT_cor_heatmap %>%
  column_to_rownames("Sample") %>%
  dplyr::select(Acetate, Butyrate, Propionate, Lactate)

#Let us group on diagnosis
groups<-meta_GCMS_PFT_cor_heatmap$Diagnosis

#You can use kendall, spearman, or pearson below:
method<-"kendall"
df<-NULL

for(i in colnames(x)){
  for(j in colnames(y)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
      if(is.null(df)){
        df<-tmp  
      }
      else{
        df<-rbind(df,tmp)
      }    
    }
  }
}

df<-data.frame(row.names=NULL,df)
colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
df$Pvalue<-as.numeric(as.character(df$Pvalue))
df$AdjPvalue<-rep(0,dim(df)[1])
df$Correlation<-as.numeric(as.character(df$Correlation))

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-5

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Env)){
    for(j in unique(df$Type)){
      sel<-df$Env==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Env)){
    sel<-df$Env==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for signifant values
df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
p <- p + geom_tile() + scale_fill_gradient2(low="#0072B5FF", mid="white", high="#BC3C29") 
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x") + theme(axis.text.y = element_text(face="italic"))
p

#### WGCNA ####
library(vegan)
library(WGCNA)
set.seed(2024)
metadata_GCMS = rownames_to_column(metadata, "Sample")
metadata_GCMS = metadata_GCMS %>%
  select(Sample, Diagnosis, `Acetic acid`,
         `Lactic acid`, `Propionic acid`, `Butyric acid`) %>%
  filter(!is.na(`Acetic acid`))
colnames(metadata_GCMS) = c("Sample", "Diagnosis", "Acetate", 
                        "Lactate", "Propionate", "Butyrate")

tax_table <- read_csv("Original-files/Taxa-df_WGCNA.csv")
ASV_genus <- tax_table %>%
  filter(!is.na(Genus))
ASV_list <- ASV_genus$...1

df <- read_csv("Original-files/OTUdf_WGCNA.csv")
df <- dplyr::rename("ASV" = "...1", df)
df <- df %>%
  filter(ASV %in% ASV_list)
df <- column_to_rownames(df, "ASV")
df <- as.data.frame(t(df))
df <- rownames_to_column(df, "Sample")

df <- df %>%
  filter(Sample %in% metadata_GCMS$Sample)
dropped_samples = setdiff(metadata_GCMS$Sample, df$Sample)
metadata_GCMS = data.frame(metadata_GCMS) %>%
  filter(!Sample %in% dropped_samples)

df = column_to_rownames(df, "Sample")

HellingerData<-decostand(df,method = "hellinger")
HD <- colSums(HellingerData)
HD <- data.frame(HD)
HD <- data.frame(HD[order(HD$HD),])
HD <- dplyr::rename("colSums" = "HD.order.HD.HD....", HD)
HD$Index <- seq(nrow(HD))
ggplot(data=HD, aes(x=Index, y=colSums)) + geom_point() 

HellingerData<-decostand(df,method = "hellinger")
HellingerData <- rownames_to_column(HellingerData, var="sample-id")

df_transformed <- left_join(metadata_GCMS[,c("Sample", "Acetate", "Lactate", "Propionate", "Butyrate")], 
                            HellingerData, by=c("Sample"="sample-id"))
df_transformed <- column_to_rownames(df_transformed, var="Sample")
df_transformed_abund <- df_transformed %>%
  dplyr::select(5:ncol(df_transformed))

gsg=goodSamplesGenes(df_transformed_abund, verbose = 3);

# Check for genes and samples with too many missing values:
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", 
                     paste(names(df_transformed_abund)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", 
                     paste(rownames(df_transformed_abund)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  df_transformed_abund = df_transformed_abund[gsg$goodSamples, gsg$goodGenes]
}

gsg=goodSamplesGenes(df_transformed_abund, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(df_transformed_abund), method = "average");
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

OTUSamples = rownames(df_transformed_abund)
traitRows = match(OTUSamples, metadata_GCMS$Sample);
datTraits = metadata_GCMS[traitRows,c("Acetate", "Lactate", "Propionate", "Butyrate")]
#metadata_GCMS_IPF = metadata_GCMS_IPF %>% filter(!Sample == 'BRU.05330')
rownames(datTraits) = metadata_GCMS$Sample
collectGarbage()

sampleTree2 = hclust(dist(df_transformed), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Code-Network Analysis
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()

powers = c(c(1:10), seq(from = 10, to=30, by=1))

sft = pickSoftThreshold(df_transformed_abund, 
                        powerVector = powers, verbose = 5, networkType = "signed")

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.35,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 9;
adjacency = WGCNA::adjacency(df_transformed_abund, power = softPower, type = "signed");
taxa_names = colnames(df_transformed_abund)
colnames(adjacency) = taxa_names
rownames(adjacency) = taxa_names

# Calculate degree of overlap between connections of two modules within a network
TOM_WGCNA = WGCNA::TOMsimilarity(adjacency, TOMType = "signed");
dimnames(TOM_WGCNA) = list(taxa_names, taxa_names)

dissTOM = 1-TOM_WGCNA

TaxaTree = hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 18;

dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM, method="hybrid",
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")

# Eigengenes merging if below threshold
MEList = moduleEigengenes(df_transformed_abund, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs); #dissimilarity of module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25 # from the tutorial

abline(h=MEDissThres, col = "red") # all are to be merged unsurprising
merge = mergeCloseModules(df_transformed_abund, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors;
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)

plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

nTaxa = ncol(df_transformed_abund);
nSamples = nrow(df_transformed_abund);

MEs0 = moduleEigengenes(df_transformed_abund, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(80,15)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(7, 15, 3, 3));

data_dist_rows <- vegdist(top_other, method = "bray")
row_clustering <- hclust(data_dist_rows, "average")

# save heatmap
library(svglite)
library(pheatmap)
p <- pheatmap(moduleTraitCor,
              cluster_rows = METree, # Use the ME dendrogram for row clustering
              cluster_cols = T,  # You can cluster traits if you want
              labels_row = c("Module1", "Module2", "Module3", "Module4", 
                             "Module5", "Module0"),
              labels_col = names(datTraits[1:4]),
              display_numbers = textMatrix,
              color = Palette_diverging1,
              main = "Module-trait relationships",
              legend = T,
              fontsize = 14, number_color = "black",
              treeheight_row = 100)
p

svglite("WGCNA.svg", width=12, height=8)
print(p)
dev.off()

## Detailed
Propionate = as.data.frame(datTraits$Propionate);
names(Propionate) = "Propionate"

modNames = substring(names(MEs), 3)
TaxaModuleMembership = as.data.frame(cor(df_transformed_abund, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
names(TaxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
TaxaTraitSignificance = as.data.frame(cor(df_transformed_abund, Propionate, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
names(TaxaTraitSignificance) = paste("GS.", names(Propionate), sep="");
names(GSPvalue) = paste("p.GS.", names(Propionate), sep="");

module = "green"
column = match(module, modNames);
moduleTaxa = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(TaxaModuleMembership[moduleTaxa, column]), 
                   abs(TaxaTraitSignificance[moduleTaxa, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Taxa significance for Propionate",
                   main = paste("Module membership vs. Taxa significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(df_transformed_abund)
names(df_transformed_abund)[moduleColors=="green"]

probes = names(df_transformed_abund)
probes2annot = match(probes, tax_table$...1)

sum(is.na(probes2annot))

TaxaInfo0 = data.frame(Taxon = probes,
                       TaxaSymbol = tax_table$...1[probes2annot],
                       LinkID = paste(tax_table$Genus[probes2annot],tax_table$Species[probes2annot]),
                       moduleColor = moduleColors,
                       TaxaTraitSignificance,
                       GSPvalue)
modOrder = order(-abs(cor(MEs, Propionate, use = "p")));

for (mod in 1:ncol(TaxaModuleMembership))
{
  oldNames = names(TaxaInfo0)
  TaxaInfo0 = data.frame(TaxaInfo0, TaxaModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(TaxaInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

TaxaOrder = order(TaxaInfo0$moduleColor, -abs(TaxaInfo0$GS.Propionate));
TaxaInfo = TaxaInfo0[TaxaOrder, ]

write.csv(TaxaInfo, file = "TaxaInfo_Propionate.csv")

##### WGCNA Genera #####
WGCNA_plot = TaxaInfo %>%
  #filter(p.GS.Propionate < 0.05) %>%
  filter(moduleColor == "green") %>%
  mutate(LinkID = gsub("NA", "", LinkID),
         moduleColor = case_when(moduleColor == "yellow" ~ 1,
                                 moduleColor == "green" ~ 2,
                                 moduleColor == "brown" ~ 3,
                                 moduleColor == "blue" ~ 4,
                                 moduleColor == "turquoise" ~ 5),
         moduleColor = as.factor(moduleColor)) %>%
  arrange(-log(p.GS.Propionate)) %>%
  mutate(Taxon = factor(Taxon, levels = unique(Taxon))) 

# Assuming WGCNA_plot has been correctly created, sorted, and Taxon is a factored column.

library(ggplot2)
library(dplyr)

ggplot(WGCNA_plot, aes(y = Taxon, x = -log(p.GS.Propionate))) +
  geom_segment(aes(y = Taxon, yend = Taxon, x = 0, xend = -log(p.GS.Propionate)), 
               color = "gray40", 
               linewidth = 1.5) + 
  geom_vline(xintercept = -log(0.05), linetype = "dashed", color = "gray40") +
  #geom_text(label = c("p=0.05"), x=c(3.1), y=c(3), size = 5, color = "gray40", angle=270)+
  geom_point(aes(color = GS.Propionate), 
             size = 14) + 
  scale_color_gradient2(low = "#FFDC91", mid = "#E18727", high = "#BC3C29", 
                        midpoint = median(WGCNA_plot$GS.Propionate),
                        name = "Pearson's\nCorrelation") +
  scale_y_discrete(labels = WGCNA_plot$LinkID) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 20, hjust = 1, face = "italic", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    y = "Genus", 
    x = "-log(p-value) (Significance Strength)"
  )

ggsave("Genus/WGCNA-sig-plot.svg", width = 10, height=12)

##### WGCNA Network analysis #####
# Find the names of all taxa in the green module
green_taxa = colnames(df_transformed_abund)[moduleColors == "green"]

# Subset the TOM to only include the taxa in the green module
green_TOM = TOM_WGCNA[green_taxa, green_taxa]
green_genus_names = tax_table$Genus[match(green_taxa, tax_table$...1)]

# Create an edge list
green_edges = exportNetworkToCytoscape(green_TOM,
                                       edgeFile = "green_module_edges.txt",
                                       nodeFile = "green_module_nodes.txt",
                                       threshold = 0.2, # You may need to adjust this threshold to get a clear plot
                                       nodeNames = green_taxa,
                                       altNodeNames = green_genus_names)

## SCFA concentration ##
taxaGS_Propionate = as.numeric(cor(df_transformed_abund, datTraits$Propionate, use = "p"))
taxaGS_Propionate_pvalue = corPvalueStudent(taxaGS_Propionate, nSamples)
# Create a data frame for your node attributes
node_attributes = data.frame(
  ASV = colnames(df_transformed_abund),
  ModuleColor = moduleColors,
  Propionate_GS = taxaGS_Propionate,
  Propionate_GS_pvalue = taxaGS_Propionate_pvalue
)

# Export the table as a tab-separated file
write.table(node_attributes, file = "all_taxa_attributes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### Top ten stacked plot ####
genus = read_csv("Genus/Genus-normalised-metadata.csv")
genus = genus %>%
  mutate(Prevotella = Prevotella + `[Prevotella]`) %>%
  dplyr::select(-`[Prevotella]`)

abund_table <- genus %>%
  dplyr::select(1,19:ncol(genus), -Unknown) #remove unknown 
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)

meta_table <- genus %>%
  dplyr::select(`sample-id`, Diagnosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

Others = 100 - rowSums(top) # remaining taxa grouped into Others. 

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(genus)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others_list <- data.frame(other[,colnames(other) %in% taxa_list2])

top_others = cbind(top, Others)
df_ten = cbind(meta_table, top_others)
wilcox.test(Streptococcus ~ Diagnosis, df_ten) 
wilcox.test(Prevotella ~ Diagnosis, df_ten) # sig
wilcox.test(Veillonella ~ Diagnosis, df_ten) 
wilcox.test(Haemophilus ~ Diagnosis, df_ten) 
wilcox.test(Gemella ~ Diagnosis, df_ten) 
wilcox.test(Actinomyces ~ Diagnosis, df_ten) 
wilcox.test(Actinobacillus ~ Diagnosis, df_ten) # sig
wilcox.test(Rothia ~ Diagnosis, df_ten) 
wilcox.test(Staphylococcus ~ Diagnosis, df_ten) 

df_long <- reshape2::melt(df_ten, id.vars = c("Diagnosis"), 
                          measure.vars = c(colnames(top_others)), 
                          variable.name = "Genus",
                          factorsAsStrings = TRUE, na.rm = TRUE)

df_long_summarised <- df_long %>%
  group_by(Diagnosis, Genus) %>%
  dplyr::summarise(n_samples =n(),
                   mean_value=mean(value),
                   sd=sd(value),
                   median=median(value),
                   q1 = quantile(value, 0.25),  # 1st quartile
                   q3 = quantile(value, 0.75)) # 3rd quartile

df_long_summarised$Diagnosis <- factor(df_long_summarised$Diagnosis, levels = c("Controls", "IPF"))

dat_text <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Controls", "IPF"),
  x = c(2,2), 
  y = c(40,38)
)

dat_text1 <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Controls", "IPF"),
  x = c(7,7), 
  y = c(8,0)
)

ggplot(df_long_summarised, aes(x=Genus, y=median)) + 
  geom_bar(aes(y = median, x = Genus, fill = Genus),
           stat="identity") + facet_wrap(.~Diagnosis, scales="fixed") + 
  geom_errorbar(aes(x=Genus, ymin=(q1), ymax=(q3)), width=0.3, color='black', linewidth=0.5) + 
  scale_fill_manual(values = Palette_10) +
  theme_classic() + 
  labs(x = "", y = "Median relative abundance (%)",
       title = "Top ten most abundant genera: Controls vs. IPF",
       subtitle = "") + 
  geom_text(data=dat_text, mapping= aes(x=x, y=y, label=label), size = 10) +
  geom_text(data=dat_text1, mapping= aes(x=x, y=y, label=label), size = 10) +
  theme(axis.text.x = element_text(face="italic", size=12, hjust=1, angle=45),
        legend.text = element_text(face = "italic"),
        legend.position = "none",
        strip.text = element_text(size=16),
        plot.margin = unit(c(10,0,10,0), "pt"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
        axis.title = element_text(size=18), axis.text = element_text(size=16),
        panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
        axis.line = element_line(color="grey"),    
  guides(fill = guide_legend(title = "Genus")))

ggsave("Genus/Top-ten.svg", width=10, height=7)

phylum = read_csv("Phylum/Phylum-normalised-metadata.csv")

abund_table <- phylum %>%
  dplyr::select(1,19:ncol(phylum)) #remove unknown 
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)

meta_table <- phylum %>%
  dplyr::select(`sample-id`, Diagnosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

df_phylum = cbind(meta_table, top)
wilcox.test(Firmicutes ~ Diagnosis, df_phylum) 
wilcox.test(Bacteroidetes ~ Diagnosis, df_phylum) # sig
wilcox.test(Proteobacteria ~ Diagnosis, df_phylum) 
wilcox.test(Actinobacteria ~ Diagnosis, df_phylum) 
wilcox.test(Fusobacteria ~ Diagnosis, df_phylum) 

df_long <- reshape2::melt(df_phylum, id.vars = c("Diagnosis"), 
                          measure.vars = c(colnames(top)), 
                          variable.name = "Phylum",
                          factorsAsStrings = TRUE, na.rm = TRUE)

df_long_summarised <- df_long %>%
  group_by(Diagnosis, Phylum) %>%
  dplyr::summarise(n_samples =n(),
                   mean_value=mean(value),
                   sd=sd(value),
                   median=median(value),
                   q1 = quantile(value, 0.25),  # 1st quartile
                   q3 = quantile(value, 0.75)) # 3rd quartile

df_long_summarised$Diagnosis <- factor(df_long_summarised$Diagnosis, levels = c("Controls", "IPF"))

ggplot(df_long_summarised, aes(x=Diagnosis, y=mean_value)) + 
  geom_bar(aes(y = mean_value, x = Diagnosis, fill = Phylum),
           stat="identity") +
  #geom_errorbar(aes(x=Genus, ymin=(q1), ymax=(q3)), width=0.3, color='black', linewidth=0.5) + 
  scale_fill_manual(values = c("#6F99AD","#8ab184", "#de6600", "firebrick", "goldenrod")) +
  theme_classic() + 
  labs(x = "", y = "Mean relative abundance (%)",
       title = "",
       subtitle = "") + 
  theme(axis.text.x = element_text(face="bold", size=12, hjust=1),
        legend.text = element_text(),
        strip.text = element_text(size=16),
        plot.margin = unit(c(10,0,10,0), "pt"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
        axis.title = element_text(size=18), axis.text = element_text(size=16),
        panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
        axis.line = element_line(color="grey"),    
        guides(fill = guide_legend(title = "Genus")))

ggsave("Phylum/Top-stacked.svg", width=10, height=8)

#### Burden and SCFA tertiles ####
metadata_GCMS = read_csv("Original-files/Propionate_metadata.csv")
metadata_GCMS = metadata_GCMS %>%
  filter(!is.na(`Acetic acid`)) %>% 
  filter(!is.na(Burden)) %>%
  select(Sample, Diagnosis, Burden, `Acetic acid`, `Lactic acid`, `Butyric acid`, `Propionic acid`) %>%
  rename("Acetate" = `Acetic acid`,
         "Lactate" = `Lactic acid`,
         "Butyrate" = `Butyric acid`,
         "Propionate" = `Propionic acid`) %>%
  column_to_rownames("Sample")
  

GCMS_tertiles = metadata_GCMS %>%
  mutate(Tertiles = ntile(Burden, 3),
         Tertiles = case_when(Tertiles == "1" ~ "Low",
                              Tertiles == "2" ~ "Medium",
                              Tertiles == "3" ~ "High"),
         Tertiles = factor(Tertiles, levels = c("Low", "Medium", "High"))) %>%
  dplyr::select(Diagnosis, Acetate, Propionate, Butyrate, Lactate, Tertiles) %>%
  rownames_to_column(., "Sample")

GCMS_tertiles = reshape2::melt(GCMS_tertiles, 
                     id.vars = c("Sample", "Diagnosis", "Tertiles"),
                     variable.name = "SCFA")

# Acetate
kruskal.test(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Acetate",])
dunnTest(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Acetate",])
p <- ggplot(GCMS_tertiles[GCMS_tertiles$SCFA=="Acetate",], aes(x=Tertiles, y=value, color=Tertiles)) + 
  geom_boxplot(outliers = F) + geom_point(size=3, aes(shape = Diagnosis)) + 
  theme_classic() + scale_color_manual(values = c("#BC3C29FF", "#6ea6a4", "#0072B5FF")) + 
  scale_shape_manual(values = c(17, 16))+
  labs(title = "Acetate",
       y = "Metabolite in BAL (µM)") + 
  annotate("text", x = 1.5, y=168, label= "***", size = 8) +
  annotate("segment", x=1, xend=2, y=165, yend=165) +
  annotate("text", x = 2, y=178, label = "***", size=8) +
  annotate("segment", x=1, xend=3, y=175, yend=175) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        strip.text.y = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(size = 20,
                             face = "bold")) + 
  guides(color = "none")

p

ggsave("Genus/Acetate-Tertiles.svg", width=8, height =6)

# Lactate
kruskal.test(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Lactate",])
dunnTest(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Lactate",])

p <- ggplot(GCMS_tertiles[GCMS_tertiles$SCFA=="Lactate",], aes(x=Tertiles, y=value, color=Tertiles)) + 
  geom_boxplot(outliers = F) + geom_point(size=3, aes(shape = Diagnosis)) + 
  theme_classic() + scale_color_manual(values = c("#BC3C29FF", "#6ea6a4", "#0072B5FF")) + 
  scale_shape_manual(values = c(17, 16))+
  labs(title = "Lactate",
    y = "Metabolite in BAL (µM)") + 
  annotate("text", x = 1.5, y=23, label= "***", size = 8) +
  annotate("segment", x=1, xend=2, y=20, yend=20)+
  annotate("text", x = 2, y=28, label = "***", size=8) +
  annotate("segment", x=1, xend=3, y=25, yend=25)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        strip.text.y = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 20)) + 
  guides(color = "none")
  
p

ggsave("Genus/Lactate-Tertiles.svg", width=8, height =6)

# Propionate
kruskal.test(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Propionate",])
dunnTest(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Propionate",])

p <- ggplot(GCMS_tertiles[GCMS_tertiles$SCFA=="Propionate",], aes(x=Tertiles, y=value, color=Tertiles)) + 
  geom_boxplot(outliers = F) + geom_point(size=3, aes(shape = Diagnosis)) + 
  theme_classic() + scale_color_manual(values = c("#BC3C29FF", "#6ea6a4", "#0072B5FF")) + 
  scale_shape_manual(values = c(17, 16))+
  labs(title = "Propionate",
       y = "Metabolite in BAL (µM)") + 
  annotate("text", x = 1.5, y=145, label= "***", size = 8) +
  annotate("segment", x=1, xend=2, y=142, yend=142)+
  annotate("text", x = 2, y=155, label = "***", size=8) +
  annotate("segment", x=1, xend=3, y=152, yend=152) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        strip.text.y = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size=20)) + 
  guides(color = "none")

p

ggsave("Genus/Propionate-Tertiles.svg", width=8, height =6)

# Butyrate
kruskal.test(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Butyrate",])
dunnTest(value ~ Tertiles, GCMS_tertiles[GCMS_tertiles$SCFA=="Butyrate",])

p <- ggplot(GCMS_tertiles[GCMS_tertiles$SCFA=="Butyrate",], aes(x=Tertiles, y=value, color=Tertiles)) + 
  geom_boxplot(outliers = F) + geom_point(size=3, aes(shape = Diagnosis)) + 
  theme_classic() + scale_color_manual(values = c("#BC3C29FF", "#6ea6a4", "#0072B5FF")) + 
  scale_shape_manual(values = c(17, 16))+
  labs(y = "Metabolite in BAL (µM)",
       title = "Butyrate") + 
  annotate("text", x = 1.5, y=50, label= "***", size = 8) +
  annotate("segment", x=1, xend=2, y=47, yend=47)+
  annotate("text", x = 2, y=58, label = "***", size=8) +
  annotate("segment", x=1, xend=3, y=55, yend=55)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        strip.text.y = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 20)) + 
  guides(color = "none")

p

ggsave("Genus/Butyrate-Tertiles.svg", width=8, height =6)

#### Heatmap with SCFA metabolites ####
set.seed(123)
genus = read_csv("Genus/Genus-normalised-metadata.csv")
meta_heatmap = genus %>%
  filter(Diagnosis=="IPF") %>%
  select(`sample-id`, Diagnosis, `Acetic acid`, `Butyric acid`,
         `Lactic acid`, `Propionic acid`) %>%
  filter(!is.na(`Acetic acid`)) %>%
  arrange(Diagnosis)
colnames(meta_heatmap) = c("Sample", "Diagnosis", "Acetate",
                           "Butyrate", "Lactate", "Propionate")
abund_genus <- genus %>%
  filter(Diagnosis == "IPF") %>%
  filter(`sample-id` %in% meta_heatmap$Sample) %>%
  dplyr::select(1,Actinomyces:ncol(genus)) %>%
  mutate(Prevotella = Prevotella + `[Prevotella]`) %>%
  select(-`[Prevotella]`, -Unknown)
abund_genus <- column_to_rownames(abund_genus, var="sample-id")
rowSums(abund_genus)

top <- abund_genus[,order(colSums(abund_genus),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
dim(top)

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_genus[,order(colSums(abund_genus),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(genus)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)

Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
rowSums(top_other)

colour_palette <- colorRampPalette(colors=c("#FFFFFF",
                                            "#406C7F",
                                            "#003f5a"))(100)

# Add annotation to heatmap AFTER distance matrix has been calculated
data_dist_rows <- vegdist(top, method = "bray")
row_clustering <- hclust(data_dist_rows, "average")

# the metadata does not influence the distance matrix. 
# only interested in: FVC, DLCO, CT, Fibrosis
annot_df1 <- data.frame(Diagnosis = meta_heatmap$Diagnosis)

col1 = list(Diagnosis = c("Controls" = "#BC3C29FF",
                          "IPF" = "#0072B5FF"))
library(ComplexHeatmap)
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Diagnosis", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 16), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 14))) # Sidebar legend label font size

annot_df2 <- data.frame(Propionate = meta_heatmap$Propionate)
col2 = list(Propionate= colorRamp2(c(0, 50, 100, 150), 
                                          c("#fec682","#E18727","#BC3C29", "#B22222")))
sidebar_annotation2 <- rowAnnotation(df = annot_df2,
                                     col=col2,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Propionate", 
                                                                    title_gp = gpar(fontsize = 16),
                                                                    labels_gp = gpar(fontsize = 14)))

annot_df3 <- data.frame(Lactate = meta_heatmap$Lactate)
col3 = list(Lactate= colorRamp2(c(0, 10, 15, 20), 
                                      c("#fec682","#E18727","#BC3C29", "#B22222")))
sidebar_annotation3 <- rowAnnotation(df = annot_df3,
                                     col=col3,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Lactate", 
                                                                    title_gp = gpar(fontsize = 16),
                                                                    labels_gp = gpar(fontsize = 14)))

annot_df4 <- data.frame(Butyrate = meta_heatmap$Butyrate)
col4 = list(Butyrate= colorRamp2(c(0, 15, 30, 50), 
                                       c("#fec682","#E18727","#BC3C29", "#B22222")))
sidebar_annotation4 <- rowAnnotation(df = annot_df4,
                                     col=col4,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Butyrate", 
                                                                    title_gp = gpar(fontsize = 16),
                                                                    labels_gp = gpar(fontsize = 14)))
annot_df5 <- data.frame(Acetate = meta_heatmap$Acetate)
col5 = list(Acetate= colorRamp2(c(0, 50, 100, 150), 
                                      c("#fec682","#E18727","#BC3C29", "#B22222")))
sidebar_annotation5 <- rowAnnotation(df = annot_df5,
                                     col=col5,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Acetate", 
                                                                    title_gp = gpar(fontsize = 16),
                                                                    labels_gp = gpar(fontsize = 14)))

sorted_samples <- meta_heatmap$Sample

sort_order <- top[sorted_samples,]

heatmap <- Heatmap(as.matrix(sort_order),  # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = F, # Cluster the rows using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 14), # Row name font size
                   column_names_gp = gpar(fontsize = 14, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 14), # Column title font size
                   row_title = "Samples", # Set row title
                   row_labels = rownames(sort_order),
                   row_names_side = c("left"),
                   row_title_gp = gpar(fontsize = 14), # Set row title font size
                   column_title = "", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 14), # Set legend title font size
                                               labels_gp = gpar(fontsize = 14))) # Set legend label font size

p <- heatmap + sidebar_annotation2 +
  sidebar_annotation3 + sidebar_annotation4 + sidebar_annotation5
p

svglite("Genus/Heatmap.svg", width = 12, height = 10)
draw(p)
dev.off()

#### Supplementary figures: BAL vs. Negatives ####
phylum <- read_tsv("Phylum/Phylum-all-samples-relative-abundance.txt") %>%
  column_to_rownames(., "...1") 
phylum = as.data.frame(t(phylum)) %>%
  rownames_to_column(., "sample-id")
meta_table = read_csv("Original-files/Propionate_metadata.csv")
meta_table = meta_table %>%
  select(Sample, Diagnosis, Burden)
df_phylum <- inner_join(meta_table, phylum, by=c("Sample"="sample-id")) %>%
  mutate(SampleType = case_when(Diagnosis == "IPF" | Diagnosis == "Controls" ~ "BAL",
                                 .default = "Negative control"))
df_phylum_meta = df_phylum %>%
  dplyr::select(Sample:Diagnosis, SampleType)

PCA <- prcomp(
  df_phylum[c("Firmicutes", "Proteobacteria")],
  scale. = F)

ggbiplot(PCA, obs.scale = 0, var.scale = 0,
         groups=df_phylum$Diagnosis, ellipse = T, circle=F,
         ellipse.linewidth = 0.4, ellipse.fill = F, ellipse.alpha = 0.2, varname.size = 4) +
  scale_color_discrete(name = "") +
  theme(legend.direction = "horizontal", legend.position = "top") + theme_classic() + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12))+
  scale_color_manual(values=c("#BC3C29FF", "#0072B5FF", "snow4")) 
ggsave("Phylum/Supplementary-biplot.svg", height=10, width=8)

BC_phylum = vegdist(df_phylum[,c("Actinobacteria", "Bacteroidetes", "Firmicutes", 
                                 "Fusobacteria", "Proteobacteria", "Verrucomicrobia")], "bray")
BC_pcoa = cmdscale(BC_phylum, k=2, eig=T)  
BC_pcoa_eig <- BC_pcoa$eig
BC_total_variance <- sum(BC_pcoa_eig[BC_pcoa_eig > 0])

percentage_explained_pco1 <- (BC_pcoa_eig[1] / BC_total_variance) * 100
percentage_explained_pco2 <- (BC_pcoa_eig[2] / BC_total_variance) * 100

BC_pcoa_coord = BC_pcoa$points

colnames(BC_pcoa_coord) = c("PCoA1", "PCoA2")
plot.data <- cbind(df_phylum_meta, BC_pcoa_coord)
plot.data = dplyr::rename("Groups" = "Diagnosis", plot.data)

ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Groups, fill = Groups),
             size = 3) + 
  stat_ellipse(aes(fill = Groups, color = Groups),
               alpha=0.2, level = 0.95, geom = "polygon") +
  # scale_"" is used to design the plot
  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "snow4")) + 
  scale_colour_manual(values = c("#BC3C29FF", "#0072B5FF", "snow4")) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "",
       x = paste("PCA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face="italic"),
    legend.position.inside = c(0.9,0.10), legend.title = element_text(size=18, face = "bold"),
    legend.text = element_text(size=16, face = "italic"), legend.position = "inside",
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))
ggsave("Phylum/Supplementary-NMDS.svg", height=10, width=12)

### Stacked bar plot ###
df_long <- reshape2::melt(df_phylum, id.vars = c("Diagnosis", "SampleType"), 
                          measure.vars = c("Actinobacteria", "Bacteroidetes",
                                           "Firmicutes", "Fusobacteria", "Proteobacteria", 
                                           "Verrucomicrobia"), 
                          variable.name = "Phylum",
                          factorsAsStrings = TRUE, na.rm = TRUE)

df_long_summarised <- df_long %>%
  group_by(SampleType, Phylum) %>%
  dplyr::summarise(n_samples =n(),
                   mean_value=mean(value),
                   sd=sd(value),
                   median=median(value),
                   q1 = quantile(value, 0.25),  # 1st quartile
                   q3 = quantile(value, 0.75)) # 3rd quartile

ggplot(df_long_summarised, aes(x=SampleType, y=mean_value)) + 
  geom_bar(aes(y = mean_value, x = SampleType, fill = Phylum),
           stat="identity") +
  #geom_errorbar(aes(x=Genus, ymin=(q1), ymax=(q3)), width=0.3, color='black', linewidth=0.5) + 
  scale_fill_manual(values = Palette_6) +
  theme_classic() + 
  labs(x = "", y = "Mean relative abundance (%)",
       title = "",
       subtitle = "") + 
  theme(axis.text.x = element_text(face="bold", size=12),
        legend.text = element_text(),
        legend.position = "bottom",
        strip.text = element_text(size=16),
        plot.margin = unit(c(10,0,10,0), "pt"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
        axis.title = element_text(size=18), axis.text = element_text(size=16),
        panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
        axis.line = element_line(color="grey"),    
        guides(fill = guide_legend(title = "Genus")))

ggsave("Phylum/Top-ten-all.svg", width=8, height=10)

##### Boxplot: Burden ####
p <- ggplot(df_phylum, aes(x=Diagnosis, y=Burden, color=Diagnosis)) + 
  geom_boxplot()
p + theme_classic() + stat_compare_means(comparisons = c("Controls", "IPF", "Negative control")) + 
  scale_color_manual(values=c("#BC3C29FF", "#0072B5FF", "snow4")) + 
  labs(y="BAL bacterial burden \n Log10(16S rRNA gene copies/mL of BAL)", x="",
       subtitle = "Kruskal-Wallis with post-hoc Dunn's") + 
  scale_y_continuous(trans='log10') + 
  geom_jitter(position=position_jitter(0.1), size = 1.5) + theme(legend.position = "none") + 
  annotate("text", x = 1.5, y=45000000, label = "p=4.32e-10", size = 3) + 
  annotate("segment", x=1, xend=2, y=35000000, yend = 35000000)+
  annotate("text", x = 2.5, y=65000000, label = "p=7.78e-11", size = 3) + 
  annotate("segment", x=2, xend=3, y=50000000, yend = 50000000)

ggsave("Genus/Burden.svg", width=4, height=6)

kruskal.test(Burden ~ Diagnosis, data = df_phylum)
dunnTest(Burden ~ Diagnosis, data = df_phylum)

