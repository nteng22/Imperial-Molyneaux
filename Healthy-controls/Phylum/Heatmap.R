#install.packages("phyloseq")
library("phyloseq")
#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("vegan")
library("vegan")
#install.packages("reshape2")
library("reshape2")

# Make a complex heatmap clustered by overall similarity of the genus microbiota using a Bray-Curtis matrix.
setwd("~/GitHub/Imperial-Molyneaux/Healthy-controls/Phylum/")

#install.packages("ComplexHeatmap")
# if (!requireNamespace("BiocManager", quietly = TRUE))

library("ComplexHeatmap")

##Select data
normalised_data <- as.data.frame(read_csv("Phylum-normalised.csv")) %>%
  filter(!Diagnosis=="Negative")
normalised_data <- dplyr::rename("Thermi" = "[Thermi]", normalised_data)
abund_table <- normalised_data %>%
  select(5:ncol(normalised_data))

metadata <- normalised_data %>%
  select("sample-id":Diagnosis)

# Make the data into percentages
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
#Extract list of top N Taxa
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(abund_table)[2] # This needs to be the total number of colums in the dataframe.
# dim(normalised_data)[2], calls the second element of the vector dim(normalised_data)
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
top_other <- cbind(metadata$`sample-id`, top_other)
top_other <- column_to_rownames(top_other, var = "metadata$`sample-id`")

top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

# Use the data to calculate the clustering for the heatmap
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows
# You can change the method e.g. Jacard, see manual for more information
data_dist_rows <- vegdist(top_other, method = "bray")
# Cluster the rows by hierarchical clustering
row_clustering <- hclust(data_dist_rows, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns
data_dist_columns <- vegdist(t(top_other), method = "bray")
# Cluster the columns by hierarchical clustering
col_clustering <- hclust(data_dist_columns, "average")

# Define the color palette for the heatmap colours
colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)

# This is to add annotation to the heatmap
# Make a dataframe of the variable you want to have 'annotated' with the heatmap
annot_df1 <- data.frame(Diagnosis = metadata$Diagnosis)

# Make a list, where the levels of the dataframe are allocated a colour
# colour list can be found here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
col1 = list(Diagnosis = c("COVID"= "hotpink3",
                          "healthy" = "darkseagreen",
                          "Bronchiectasis" = "deepskyblue2",
                          "Asthma" = "plum",
                          "COPD" = "gold",
                          "Post infective" = "red",
                          "Eosinophilc Pneumonia" = "lightblue",
                          "Empysema" = "lightgreen",
                          "Adenocarcinoma" = "slategrey",
                          "Aspiration no fibrois" = "orchid4",
                          "Antiphospholipid Syndrome" = "orange2"))

# Add everything together, change the df, col and title accordingly. 
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Diagnosis", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

# Create the heatmap
heatmap <- Heatmap(as.matrix(top_other), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   #cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 6), # Row name font size
                   column_names_gp = gpar(fontsize = 8, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 8), # Column title font size
                   row_title = "Samples", # Set row title
                   row_title_gp = gpar(fontsize = 8), # Set row title font size
                   column_title = "Phylum", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 8), # Set legend title font size
                                               labels_gp = gpar(fontsize = 8))) # Set legend label font size

# Combine the heatmap and the annotation together in the order in which they are to appear
p <- heatmap# + sidebar_annotation1
p

pdf(file = "Heatmap-unfiltered-Phylum", width = 11, height = 8)

# print(p) saves the figure into a file
print(p)
dev.off()
