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

#install.packages("ComplexHeatmap")
# if (!requireNamespace("BiocManager", quietly = TRUE))

library("ComplexHeatmap")

setwd("/GitHub/Imperial-Molyneaux/COVID-BAL-OralRinse/Filtered data/")
data <- read_csv("../Genus-normalised-metadata.csv") %>%
  filter(!Status == "control") # filter out control

# Select abundance data and put it into a new variable
abund_table <- data %>%
  select(5:ncol(data))
meta_table <- data %>%
  select(index:Status)

# Make the data into percentages
abund_table <- abund_table/rowSums(abund_table) * 100

# Just check all the rows sum to 100.
rowSums(abund_table)

# Use the data to calculate the clustering for the heatmap
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows
# You can change the method e.g. Jacard, see manual for more information
data_dist_rows <- vegdist(abund_table, method = "bray")
# Cluster the rows by hierarchical clustering
row_clustering <- hclust(data_dist_rows, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns
data_dist_columns <- vegdist(t(abund_table), method = "bray")
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
annot_df1 <- data.frame(Status = meta_table$Status)
annot_df2 <- data.frame(Sample.type = meta_table$`Sample-type`)

# Make a list, where the levels of the dataframe are allocated a colour
# colour list can be found here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
col1 = list(Status = c("disease"= "purple",
                           "healthy" = "goldenrod1"))
col2 = list(Sample.type = c("BAL" = "maroon",
                              "Oral" = "orange3"))

# Add everything together, change the df, col and title accordingly. 
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Status", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation2 <- rowAnnotation(df = annot_df2, # Dataframe containing treatment groups
                                     col = col2, # The list of treatment groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Sample type", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
# Create the heatmap
heatmap <- Heatmap(as.matrix(abund_table), # The dataframe containing the heatmap data
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
                   column_title = "Genus", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 8), # Set legend title font size
                                               labels_gp = gpar(fontsize = 8))) # Set legend label font size

# Combine the heatmap and the annotation together in the order in which they are to appear
p <- heatmap + sidebar_annotation1 + sidebar_annotation2
p

pdf(file = "Heatmap-proportionalASV-Class.pdf", width = 11, height = 8)

# print(p) saves the figure into a file
print(p)
dev.off()
