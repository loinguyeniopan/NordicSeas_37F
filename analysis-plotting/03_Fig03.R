# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(gridExtra)
library(vegan)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(ggrepel)
library(webr)

# Load required datasets
# Read and TRANSPOSE the ASV table, then convert it to a matrix
ASV_table <- read.csv("ASV_table_37F.csv", row.names = 1, header = TRUE)
com <- as.matrix(t(ASV_table))
# Read the metadata CSV file into a dataframe 
metadata <- read.csv("nordic_metadata.csv", row.names = 1, header = TRUE)
# Read the taxon assignment table CSV file into a dataframe 
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)
# Add the "sample_type" as the first column of the ASV_table dataframe
metadata_id <- metadata %>% rownames_to_column("sample_id")

# Filter the data based on the benthic group (including unsigned) and extract row names
benthic_ASV <- rownames(taxon[taxon$group %in% c("benthic", "unassigned"),])

# Subset the 'com' table for the benthic group
benthic_com <- com[rownames(com) %in% benthic_ASV , ]

# Create a new column 'merge_taxon' by combining 'taxon1' and 'taxon2'
taxon <- taxon %>%
  mutate(merge_taxon = paste(taxon1, taxon2, sep = "_"))

# Extract 'merge_taxon' based on the row names of benthic_com
merge_taxon_column <- taxon[rownames(benthic_com), "merge_taxon"]

# Add 'merge_taxon' as the first column of benthic_com
benthic_com <- cbind(merge_taxon = merge_taxon_column, benthic_com)
benthic_com <- as.data.frame(benthic_com)


# Convert columns (except 'merge_taxon') to numeric
benthic_com <- benthic_com %>%
  mutate(across(-merge_taxon, ~ as.numeric(as.character(.))))

# Group by 'merge_taxon', sum values, and set row names
benthic_summarized <- benthic_com %>%
  group_by(merge_taxon) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("merge_taxon")

# Merge all other Monothalamea 
merged_row <- benthic_summarized %>%
  filter(rownames(benthic_summarized) %in% c("Monothalamea_Monothalamea_X", "Monothalamea_Xenophyophoroidea", "Monothalamea_", "Nodosariata_Polymorphinida")) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

rownames(merged_row) <- "Monothalamea_Monothalamea_X"

# Remove original rows and add the merged row back
benthic_summarized <- benthic_summarized[!rownames(benthic_summarized) %in% c("Monothalamea_Monothalamea_X",
                                                                              "Monothalamea_Xenophyophoroidea",
                                                                              "Monothalamea_",
                                                                              "Nodosariata_Polymorphinida"), ] %>%
  bind_rows(merged_row)


# Calculate percentage and replace values < 2 with 0
benthic_percentage <- benthic_summarized / colSums(benthic_summarized)[col(benthic_summarized)] * 100
benthic_percentage[benthic_percentage < 2] <- 0

# Remove rows with all zeros
benthic_percentage <- benthic_percentage[rowSums(benthic_percentage) != 0, ]


# Add a new row for remaining percentages
new_row <- data.frame(t(100 - colSums(benthic_percentage)))
rownames(new_row) <- "Other < 2%"
benthic_percentage <- rbind(benthic_percentage, new_row)

# Set read of missing samples to 0
benthic_percentage[c("NS02_01", "NS38_01", "NS47_01", "KH08_01", "KH08_02", "KH08_03", 
                     "KH04_02", "KH04_03", "KH05_04", "KH07_03")] <- 0

# Convert row names to a new first column
benthic_percentage <- rownames_to_column(benthic_percentage, var = "taxon")


# reshape the data to long format
taxon_data_long <- pivot_longer(benthic_percentage , cols = -taxon, names_to = "sample", values_to = "percentage")

# Define taxon names and color palette
taxon_name <- c("Globothalamea_Robertinida",   "Globothalamea_Rotaliida",     "Globothalamea_Textulariida", 
                "Monothalamea_CladeA",         "Monothalamea_CladeB",         "Monothalamea_CladeBM",       
                "Monothalamea_CladeC",         "Monothalamea_CladeD",         "Monothalamea_CladeE",        
                "Monothalamea_CladeG",         "Monothalamea_CladeL_O",       "Monothalamea_CladeV",        
                "Monothalamea_CladeY",         "Monothalamea_ENFOR2",         "Monothalamea_ENFOR3",        
                "Monothalamea_ENFOR5",         "Monothalamea_ENFOR8",
                "Monothalamea_Monothalamea_X", "Tubothalamea_Miliolida",  "Other < 2%",
                "unassigned_unassigned" )       
colors <- c("#1F78B4", "#289BE9", "#31BFFF",
            "#FF7F00", "#FF850C", "#FF8C19", "#FF9326", "#FF9A33",
            "#FFA13F", "#FFA84C", "#FFAF59", "#FFB566", "#FFBC73",
            "#FFC37F", "#FFCA8C", "#FFD199", "#FFD8A6", "#FFDFB3",
            "#ae2012", "#0a9396", "gray80")
color_mapping <- setNames(colors, taxon_name)

# Fig. 3A - List of SW samples for filtering
lists1 <- c("NS06_01", "NS02_01", "NS04_01", "NS01_01", "NS10_01", "NS47_01", "NS12_01", "NS45_01", 
            "NS43_01", "NS38_01", "HR07_01", "HR04_01", "IS01_01", "KV01_01", "KH01_01", "IS02_01", 
            "KH10_01", "NC01_01", "KH05_01", "KH02_01", "KH03_01", "KH04_ICE1", "KH04_ICE2", "KH09_01", 
            "KH08_01", "KH07_01")
           

# Filtering rows based on 'sample' column matching elements in lists1
filtered_data <- taxon_data_long %>% filter(sample %in% lists1)
lists1 <- rev(lists1)

# create the bar chart
p1 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists1) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()

p1


# Fig. 3B - List of 100mW samples for filtering
lists2 <- c("NS06_02", "NS02_02", "NS04_02", "NS01_02", "NS10_02", "NS47_02", "NS12_02", "NS45_02", 
            "NS43_02", "NS38_02", "HR07_02", "HR04_02", "IS01_02", "KV01_02", "KH01_02", "IS02_02", 
            "KH10_02", "NC01_02", "KH05_02", "KH02_02", "KH03_02", "KH04_02", "KH04_03", "KH09_02", 
            "KH08_02", "KH07_02")

# Filtering rows based on 'sample' column matching elements in lists2
filtered_data <- taxon_data_long %>% filter(sample %in% lists2)
lists2 <- rev(lists2)

# create the bar chart
p2 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists2) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping ) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()
p2

# Fig. 3C - List of BW samples for filtering
lists3 <- c("NS06_03", "NS02_03", "NS04_03", "NS01_03", "NS10_03", "NS47_03", "NS12_03", "NS45_03", 
            "NS43_03", "NS38_03", "HR07_03", "HR04_03", "IS01_03", "KV01_03", "KH01_03", "IS02_03", 
            "KH10_03", "NC01_03", "KH05_03", "KH02_03", "KH03_03", "KH04_02", "KH04_03", "KH09_03", 
            "KH08_03", "KH07_03")


# Filtering rows based on 'sample' column matching elements in lists3
filtered_data <- taxon_data_long %>% filter(sample %in% lists3)
lists3 <- rev(lists3)

# create the bar chart
p3 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists3) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()

p3

# Fig. 3D - List of SED samples for filtering
lists4 <- c("NS06_04", "NS02_04", "NS04_04", "NS01_04", "NS10_04", "NS47_04", "NS12_04", "NS45_04", 
            "NS43_04", "NS38_04", "HR07_04", "HR04_04", "IS01_04", "KV01_04", "KH01_04", "IS02_04", 
            "KH10_04", "NC01_04", "KH05_04", "KH02_04", "KH03_04", "KH04_02", "KH04_03", "KH09_04", 
            "KH08_04", "KH07_04")

# Filtering rows based on 'sample' column matching elements in lists2
filtered_data <- taxon_data_long %>% filter(sample %in% lists4)
lists4 <- rev(lists4)

# create the bar chart
p4 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists4) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()
p4


# Sum total reads across columns in benthic_summarized and convert to a data frame
total_read <- as.data.frame(t(colSums(benthic_summarized)))
rownames(total_read) <- "total"

# Set read of missing samples to 1
total_read[c("NS02_01", "NS38_01", "NS47_01", "KH08_01", "KH08_02", "KH08_03", 
             "KH04_02", "KH04_03", "KH05_04", "KH07_03")] <- 1
total_read <- as.data.frame(t(total_read))

# Convert row names to a new first column
total_read <- rownames_to_column(total_read, var = "station")


# Fig. 3A - Filtering rows based on lists1
filtered_data <- total_read %>% filter(station %in% lists1)

# Basic bar plot using ggplot2
p5 <- ggplot(filtered_data, aes(x = factor(station, levels = lists1), y = log10(total))) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() + ylim (0,6) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p5


# Fig. 3B - Filtering rows based on lists2
filtered_data <- total_read %>% filter(station %in% lists2)

# Basic bar plot using ggplot2
p6 <- ggplot(filtered_data, aes(x = factor(station, levels = lists2), y = log10(total))) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() + ylim (0,6) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p6


# Fig. 3C - Filtering rows based on lists3
filtered_data <- total_read %>%  filter(station %in% lists3)

# Basic bar plot using ggplot2
p7 <- ggplot(filtered_data, aes(x = factor(station, levels = lists3), y = log10(total))) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() + ylim (0,6) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p7



# Fig. 3D - Filtering rows based on lists4
filtered_data <- total_read %>% filter(station %in% lists4)

# Basic bar plot using ggplot2
p8 <- ggplot(filtered_data, aes(x = factor(station, levels = lists4), y = log10(total))) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() + ylim (0,6) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p8



###################################################################################

# Fig S2. supplement - ASV richness

# Convert benthic abundance to benthic foram richness
# Apply the condition: if value > 1, set as 1; otherwise, set as 0
benthic_richness <- benthic_com
benthic_richness[] <- lapply(benthic_richness, function(x) {
  if(is.numeric(x)) {
    x[x >= 1] <- 1    # Set values >= 1 to 1 (present)
    x[x < 1] <- 0     # Set values < 1 to 0 (absent)
  }
  return(x)
})

# Convert columns (except 'merge_taxon') to numeric
benthic_richness <- benthic_richness %>%
  mutate(across(-merge_taxon, ~ as.numeric(as.character(.))))

# Group by 'merge_taxon', sum values, and set row names
benthic_summarized <- benthic_richness %>%
  group_by(merge_taxon) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("merge_taxon")

# Merge all other Monothalamea 
merged_row <- benthic_summarized %>%
  filter(rownames(benthic_summarized) %in% c("Monothalamea_Monothalamea_X", "Monothalamea_Xenophyophoroidea", "Monothalamea_", "Nodosariata_Polymorphinida")) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

rownames(merged_row) <- "Monothalamea_Monothalamea_X"

# Remove original rows and add the merged row back
benthic_summarized <- benthic_summarized[!rownames(benthic_summarized) %in% c("Monothalamea_Monothalamea_X",
                                                                              "Monothalamea_Xenophyophoroidea",
                                                                              "Monothalamea_",
                                                                              "Nodosariata_Polymorphinida"), ] %>%
  bind_rows(merged_row)


# Calculate percentage and replace values < 2 with 0
benthic_percentage <- benthic_summarized / colSums(benthic_summarized)[col(benthic_summarized)] * 100
benthic_percentage[benthic_percentage < 2] <- 0

# Remove rows with all zeros
benthic_percentage <- benthic_percentage[rowSums(benthic_percentage) != 0, ]


# Add a new row for remaining percentages
new_row <- data.frame(t(100 - colSums(benthic_percentage)))
rownames(new_row) <- "Other < 2%"
benthic_percentage <- rbind(benthic_percentage, new_row)

# Set read of missing samples to 0
benthic_percentage[c("NS02_01", "NS38_01", "NS47_01", "KH08_01", "KH08_02", "KH08_03", 
                     "KH04_02", "KH04_03", "KH05_04", "KH07_03")] <- 0

# Convert row names to a new first column
benthic_percentage <- rownames_to_column(benthic_percentage, var = "taxon")


# Fig. SA, filtering rows based on lists1
filtered_data <- taxon_data_long %>% filter(sample %in% lists1)

# create the bar chart
p1 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists1) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()

p1


# Fig. SB, filtering rows based on lists2
filtered_data <- taxon_data_long %>% filter(sample %in% lists2)

# create the bar chart
p2 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists2) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping ) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()
p2


# Fig. SC, filtering rows based on lists3
filtered_data <- taxon_data_long %>% filter(sample %in% lists3)

# create the bar chart
p3 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists3) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()

p3

# Fig. SD, filtering rows based on lists3
filtered_data <- taxon_data_long %>% filter(sample %in% lists4)


# create the bar chart
p4 <- ggplot(filtered_data, aes(x = factor(sample, levels = lists4) ,  y = percentage, fill = taxon)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = scales::percent) +
  labs( x = "", y = "Read proportion",
        fill = "Taxon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none") + coord_flip()
p4



# Sum total ASV across columns in benthic_summarized and convert to a data frame
total_ASV<- as.data.frame(t(colSums(benthic_summarized)))
rownames(total_ASV) <- "total"

# Set read of missing samples to 1
total_ASV[c("NS02_01", "NS38_01", "NS47_01", "KH08_01", "KH08_02", "KH08_03", 
             "KH04_02", "KH04_03", "KH05_04", "KH07_03")] <- 0
total_ASV <- as.data.frame(t(total_ASV))

# Convert row names to a new first column
total_ASV <- rownames_to_column(total_ASV, var = "station")


# Fig. SA - Filtering rows based on lists1
filtered_data <- total_ASV %>% filter(station %in% lists1)

# Basic bar plot using ggplot2
p5 <- ggplot(filtered_data, aes(x = factor(station, levels = lists1), y = total)) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p5

# Fig. SB - Filtering rows based on lists2
filtered_data <- total_ASV %>% filter(station %in% lists2)


# Basic bar plot using ggplot2
p6 <- ggplot(filtered_data, aes(x = factor(station, levels = lists2), y = total)) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p6


# Fig. SC - Filtering rows based on lists3
filtered_data <- total_ASV %>% filter(station %in% lists3)

# Basic bar plot using ggplot2
p7 <- ggplot(filtered_data, aes(x = factor(station, levels = lists3), y = total)) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p7

# Fig. SD - Filtering rows based on lists4
filtered_data <- total_ASV %>% filter(station %in% lists4)

# Basic bar plot using ggplot2
p8 <- ggplot(filtered_data, aes(x = factor(station, levels = lists4), y = total)) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.5) +  # Create a bar plot
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip()

p8




