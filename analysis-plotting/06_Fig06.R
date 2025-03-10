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


# Read the ASV table CSV file into a dataframe with row names and headers
pc <- read.csv("ASV_table_37F.csv", row.names = 1, header = TRUE)
com <- as.matrix(pc)

# Remove samples with N < 
com <- com[rowSums(com) >= 10000, ]

set.seed(123)

#Normalize the dataset

N <- rowSums(com)
S <- specnumber(com)
min_seq <- min(N)

richness <- rarefy(com, min(N))
com <- rrarefy(com, min(N))


# Remove rows for KH04_ICE1 and KH04_ICE2
com <- com[!rownames(com) %in% c("KH04_ICE1", "KH04_ICE2"), ]


com <- as.matrix(t(com))


# Read the metadata CSV file into a dataframe with row names and headers
metadata <- read.csv("nordic_metadata.csv", row.names = 1, header = TRUE)


# Read the ASV table CSV file into a dataframe with row names and headers
taxon <- read.csv("ASV_assignment.csv",row.names = 1, header = TRUE)


# Add the "sample_type" as the first column of the pc dataframe
metadata_id <- metadata %>% rownames_to_column("sample_id")

# Filter the data based on the group and extract row names
benthic_ASV <- rownames(taxon[taxon$group %in% c("benthic", "unassigned"),])

# Subset the 'com' table to keep only rows with row names in benthic_planktonic_list
benthic_com <- com[rownames(com) %in% benthic_ASV, ]


# Calculate total reads for each ASV and prepare selected_ASV for merging
selected_ASV <- as.data.frame(rowSums(benthic_com))
colnames(selected_ASV) <- "total_read"
selected_ASV <- rownames_to_column(selected_ASV, "ASV_id")

# Merge with taxon data
taxon$ASV_id <- rownames(taxon)
selected_ASV <- merge(selected_ASV, taxon, by = "ASV_id", all.x = TRUE)

# Filter merged data based on criteria
selected_ASV <- selected_ASV[!is.na(selected_ASV$mean_similarity) &
                               selected_ASV$mean_similarity > 0.97, ]


benthic_com <- as.data.frame(benthic_com)

benthic_com <- rownames_to_column(benthic_com, "ASV_id")

# Merge selected_ASV2 with benthic_percentage
selected_ASV2 <- merge(selected_ASV, benthic_com, by = "ASV_id", all.x = TRUE)

selected_ASV2 <- selected_ASV2 %>%
  mutate(merged_taxon = paste(taxon1, taxon5, sep = ";")) %>%
  dplyr::select(ASV_id, merged_taxon, everything())

write.csv(selected_ASV2, "selected_ASV2.csv")

# Get the column names from the 10th column to the last column
remaining_columns <- colnames(selected_ASV2)[11:ncol(selected_ASV2)]

# Combine the "merged_taxon" column with the remaining columns
extracted_data <- selected_ASV2[, c("merged_taxon", "total_read", remaining_columns)]

# Group by merged_taxon and summarize the numeric columns by summing their values
merged_data <- extracted_data %>%
  group_by(merged_taxon) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) sum(.) else first(.)))

write.csv(merged_data, "merged_data.csv")


# Load required libraries
library(ComplexHeatmap)
library(circlize)



selected_taxon <- read.csv("merged_data2.csv", header = TRUE, row.names = 1)
selected_taxon <- log1p(selected_taxon)



# Define the desired column order
station_order <- c(
  "NS06_01",  "NS04_01",  "NS01_01",  "NS10_01",  "NS12_01",  "NS45_01",  "NS43_01",
  "NS38_01",  "HR07_01",  "HR04_01",  "IS01_01",  "KV01_01",  "KH01_01",  "IS02_01",  "KH10_01",  "NC01_01",
  "KH05_01",  "KH02_01",  "KH03_01",  "KH09_01",  "KH07_01", 
  "NS06_02",  "NS02_02",  "NS04_02",  "NS01_02",  "NS10_02",  "NS12_02",  "NS45_02",  "NS43_02",
  "NS38_02",  "HR07_02",  "HR04_02",  "IS01_02",  "KV01_02",  "KH01_02",  "IS02_02",  "KH10_02",  "NC01_02",
  "KH05_02",  "KH02_02",  "KH03_02",  "KH09_02",  "KH07_02",
  "NS06_03",  "NS02_03",  "NS04_03",  "NS01_03",  "NS10_03",  "NS12_03",  "NS45_03",  "NS43_03",
  "NS38_03",  "HR07_03",  "HR04_03",  "IS01_03",  "KV01_03",  "KH01_03",  "IS02_03",  "KH10_03",  "NC01_03",
  "KH05_03",  "KH02_03",  "KH03_03",  "KH09_03",  
  "NS06_04",  "NS02_04",  "NS04_04",  "NS01_04",  "NS10_04",  "NS47_04",  "NS12_04",  "NS45_04",  "NS43_04",
  "NS38_04",  "HR07_04",  "HR04_04",  "IS01_04",  "KV01_04",  "KH01_04",  "IS02_04",  "KH10_04",  "NC01_04",
  "KH02_04",  "KH03_04",  "KH09_04",  "KH08_04",  "KH07_04")



group_labels <- c(rep("01", 8), 
                  rep("02", 6),
                  rep("03", 2), 
                  rep("04", 3),
                  rep("05", 2),
                  rep("06", 9), 
                  rep("07", 6),
                  rep("08", 2),
                  rep("10", 3),
                  rep("11", 2),
                  rep("12", 9),
                  rep("13", 6), 
                  rep("14", 2), 
                  rep("15", 3),
                  rep("16", 1),
                  rep("17", 10),
                  rep("18", 6), 
                  rep("19", 2), 
                  rep("20", 2),
                  rep("21", 3))
                  
                  
colnames(selected_taxon)

# Reorder the columns in selected_taxon based on the desired order
selected_taxon <- selected_taxon[, station_order]

selected_taxon <- t(selected_taxon)

head (selected_taxon)

col_scale  = colorRamp2(c(12, 10, 4, 2, 0), c("#DF1423",
                                              "#DA5331",
                                              "#1C48A4",
                                              "#5BA0D7",
                                              "gray95"))


# Read the contents of the file
file_contents <- readLines("list.txt")

# Print the contents of the file
print(file_contents)


# Create HeatmapAnnotation object for station annotation
ha <- rowAnnotation(type_sample = file_contents, border = TRUE,
                    col = list(type_sample = c("surface_water" = "#19A9A9",
                                               "100m_water" = "#053F5C",
                                               "bottom_water" = "#F2BC52",
                                               "sediment" = "#F68E51")),
                    annotation_legend_param = list(
                      type_sample = list(
                        title = "Type of sample",
                        at = c("surface_water", "100m_water", "bottom_water",  "sediment"),
                        labels = c("Surface water", "100m water", "Bottom water", "Surface sediment"),
                        title_position = "lefttop-rot"
                      )
                    )
)

# Perform clustering on columns of selected_taxon
dend = as.dendrogram(hclust(dist(t(selected_taxon))))

# Color branches based on 4 clusters
dend = dendextend::color_branches(dend, k = 6)

p2 <- Heatmap(selected_taxon, 
              cluster_columns = dend,
              cluster_rows = FALSE,
              row_split = group_labels,  
              border_gp = gpar(col = "black", lty = 2),
              rect_gp = gpar(col = "white", lwd = 0.1),
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 10),
              col = col_scale,
              left_annotation = ha,
              heatmap_legend_param = list(
                title = "log of nornalized reads",
                legend_height = unit(4, "cm"),
                title_position = "lefttop-rot"),
)

p2



# Table S3
# Read and convert it to a matrix
ASV_table <- read.csv("ASV_table_37F.csv", row.names = 1, header = TRUE)

# Remove stations with sea-ice KH04
ASV_table <- ASV_table[!grepl("^KH04_", rownames(ASV_table)), ]

com <- as.matrix(ASV_table)

# Read the metadata CSV file into a dataframe 
metadata <- read.csv("nordic_metadata.csv", row.names = 1, header = TRUE)
# Read the taxon assignment table CSV file into a dataframe 
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)

# Add the "sample_type" as the first column of the ASV_table dataframe
metadata_id <- metadata %>% rownames_to_column("sample_id")


merged_asv_list <- read.csv("merged_asv_list.csv", header = TRUE)


# Extract benthic ASVs and match with water-sediment shared ASVs
benthic_ASV <- rownames(taxon[taxon$group %in% c("benthic", "unassigned"),])
t_com <- t(com)
benthic_com <- t_com[rownames(t_com) %in% benthic_ASV , ]

# Filter the benthic_com dataframe to include only rows where the ASV is in list_ASV
list_ASV <- merged_asv_list$ASV
filtered_benthic_com <- benthic_com[rownames(benthic_com) %in% list_ASV, ]
filtered_benthic_com <- as.data.frame(filtered_benthic_com)


# Convert to long format, with row names as 'ASV' and columns as 'station'
long_benthic_com <- filtered_benthic_com %>%
  tibble::rownames_to_column(var = "ASV") %>%
  pivot_longer(cols = -ASV, names_to = "sample_id", values_to = "abundance")

# Perform a left join to add columns from metadata_id to long_benthic_com based on sample_id
long_benthic_com <- long_benthic_com %>%
  left_join(metadata_id %>% select(sample_id, station, 
                                   sample_type, location,
                                   distance, water_depth, sed_depth), by = "sample_id")

# View the updated data frame
head(long_benthic_com)
write.csv(long_benthic_com, "long_benthic_com.csv")

long_benthic_com <- long_benthic_com %>%
  filter(abundance >= 100) %>%
  filter(sed_depth >= 200)


summary_distance <- long_benthic_com %>%
  group_by(ASV,sample_type) %>%
  summarise(
    min_distance = min(distance, na.rm = TRUE),
    max_distance = max(distance, na.rm = TRUE),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_distance)



# Read the taxon assignment table CSV file into a dataframe 
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)
taxon <- taxon  %>% rownames_to_column("ASV")

summary_distance <- summary_distance %>%
  left_join(taxon %>% select(ASV, taxon5, mean_similarity), by = "ASV")


summary_distance <- summary_distance %>%
  filter(mean_similarity >= 0.97)


write.csv(summary_distance, "summary_distance.csv")











