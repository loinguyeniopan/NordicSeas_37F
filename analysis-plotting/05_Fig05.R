# Load required packages
library(VennDiagram)
library(dplyr)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(ggrepel)
library(webr)
library(geosphere)

# Load required datasets
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

ASV_sample_type <- ASV_table %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata_id %>% select(sample_id, sample_type), by = "sample_id") %>%
  select(sample_id, sample_type, everything()) %>%
  column_to_rownames("sample_id")

# Group by the sample type and calculate the sum of each ASV
totals_by_ASV <- ASV_sample_type %>%
  group_by(sample_type) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ]


# Convert totals_by_ASV to numeric, excluding the row names
venn_ASV <- as.data.frame(apply(totals_by_ASV, 2, as.numeric))

# Keep the row names and reorder the columns
rownames(venn_ASV) <- rownames(totals_by_ASV)
venn_ASV <- venn_ASV[, c("surface_water", "100m_water", "bottom_water", "sediment")]


# Subset the data to keep only the ASVs with total greater than 0 in all samples
subset_data <- venn_ASV[rowSums(venn_ASV) > 0, ]


# Create a list of the ASV names in each circle
asv_lists <- list(
  A = rownames(subset_data)[subset_data[,"surface_water"] > 0], 
  B = rownames(subset_data)[subset_data[,"sediment"] > 0], 
  C = rownames(subset_data)[subset_data[,"100m_water"] > 0], 
  D = rownames(subset_data)[subset_data[,"bottom_water"] > 0] 
)

sample_color <- c("#19A9A9", "#053F5C", "#F2BC52", "#F68E51")

# Create the Venn diagram
venn.plot <- draw.quad.venn(
  area1 = length(asv_lists$A),
  area2 = length(asv_lists$B),
  area3 = length(asv_lists$C),
  area4 = length(asv_lists$D),
  n12 = length(intersect(asv_lists$A, asv_lists$B)),
  n13 = length(intersect(asv_lists$A, asv_lists$C)),
  n14 = length(intersect(asv_lists$A, asv_lists$D)),
  n23 = length(intersect(asv_lists$B, asv_lists$C)),
  n24 = length(intersect(asv_lists$B, asv_lists$D)),
  n34 = length(intersect(asv_lists$C, asv_lists$D)),
  n123 = length(intersect(asv_lists$A, intersect(asv_lists$B, asv_lists$C))),
  n124 = length(intersect(asv_lists$A, intersect(asv_lists$B, asv_lists$D))),
  n134 = length(intersect(asv_lists$A, intersect(asv_lists$C, asv_lists$D))),
  n234 = length(intersect(asv_lists$B, intersect(asv_lists$C, asv_lists$D))),
  n1234 = length(intersect(asv_lists$A, intersect(asv_lists$B, intersect(asv_lists$C, asv_lists$D)))),
  category = c("Surface water \n (SW)",
               "Sediment \n (SED)",
               "100m water \n (100mW)",
               "Bottom water \n (SED)"),
  fill = c("#19A9A9", "#F68E51", "#053F5C", "#F2BC52"),
  alpha = rep(0.4, 4),
  lty = "blank",
  cex = 2,
  fontface = "bold",
  cat.cex = 1.5
)


# Calculate read percentages for each area
read_perc <- round(100 * colSums(subset_data) / sum(subset_data), 1)

A = rownames(subset_data)[subset_data[,1] > 0]
B = rownames(subset_data)[subset_data[,4] > 0]
C = rownames(subset_data)[subset_data[,2] > 0]
D = rownames(subset_data)[subset_data[,3] > 0] 

# Get the ASV lists for each overlap area
AB <- intersect(asv_lists$A, asv_lists$B)
AC <- intersect(asv_lists$A, asv_lists$C)
AD <- intersect(asv_lists$A, asv_lists$D)
BC <- intersect(asv_lists$B, asv_lists$C)
BD <- intersect(asv_lists$B, asv_lists$D)
CD <- intersect(asv_lists$C, asv_lists$D)
ABC <- intersect(AB, asv_lists$C)
ABD <- intersect(AB, asv_lists$D)
ACD <- intersect(AC, asv_lists$D)
BCD <- intersect(BC, asv_lists$D)
ABCD <- intersect(ABC, asv_lists$D)

# Get the unique ASVs for each area
A_u <- setdiff(asv_lists$A, union(B, union(C, D)))
B_u <- setdiff(asv_lists$B, union(A, union(C, D)))
C_u <- setdiff(asv_lists$C, union(A, union(B, D)))
D_u <- setdiff(asv_lists$D, union(A, union(B, C)))

AB_u <- setdiff(setdiff(union(A, B), union(C, D)), union(A_u,B_u))
AC_u <- setdiff(setdiff(union(A, C), union(B, D)), union(A_u,C_u))
AD_u <- setdiff(setdiff(union(A, D), union(C, B)), union(A_u,D_u))
BC_u <- setdiff(setdiff(union(B, C), union(A, D)), union(B_u,C_u))
BD_u <- setdiff(setdiff(union(B, D), union(A, C)), union(B_u,D_u))
CD_u <- setdiff(setdiff(union(C, D), union(A, B)), union(C_u,D_u))

ABC_u <- setdiff(ABC, ABCD)
ABD_u <- setdiff(ABD, ABCD)
ACD_u <- setdiff(ACD, ABCD)
BCD_u <- setdiff(BCD, ABCD)
ABCD


# Get the read counts for the unique ASVs in each area
asv_table <- subset_data %>%
  as.data.frame() %>%
  rownames_to_column("ASV")

asv_list_table_A_u <- asv_table[asv_table$ASV %in% A_u,]
asv_list_reads_A_u <- sum(asv_list_table_A_u[,2:ncol(asv_list_table_A_u)])

asv_list_table_B_u <- asv_table[asv_table$ASV %in% B_u,]
asv_list_reads_B_u <- sum(asv_list_table_B_u[,2:ncol(asv_list_table_B_u)])

asv_list_table_C_u <- asv_table[asv_table$ASV %in% C_u,]
asv_list_reads_C_u <- sum(asv_list_table_C_u[,2:ncol(asv_list_table_C_u)])

asv_list_table_D_u <- asv_table[asv_table$ASV %in% D_u,]
asv_list_reads_D_u <- sum(asv_list_table_D_u[,2:ncol(asv_list_table_D_u)])

asv_list_table_AB_u <- asv_table[asv_table$ASV %in% AB_u,]
asv_list_reads_AB_u <- sum(asv_list_table_AB_u[,2:ncol(asv_list_table_AB_u)])

asv_list_table_AC_u <- asv_table[asv_table$ASV %in% AC_u,]
asv_list_reads_AC_u <- sum(asv_list_table_AC_u[,2:ncol(asv_list_table_AC_u)])

asv_list_table_AD_u <- asv_table[asv_table$ASV %in% AD_u,]
asv_list_reads_AD_u <- sum(asv_list_table_AD_u[,2:ncol(asv_list_table_AD_u)])

asv_list_table_BC_u <- asv_table[asv_table$ASV %in% BC_u,]
asv_list_reads_BC_u <- sum(asv_list_table_BC_u[,2:ncol(asv_list_table_BC_u)])

asv_list_table_BD_u <- asv_table[asv_table$ASV %in% BD_u,]
asv_list_reads_BD_u <- sum(asv_list_table_BD_u[,2:ncol(asv_list_table_BD_u)])

asv_list_table_CD_u <- asv_table[asv_table$ASV %in% CD_u,]
asv_list_reads_CD_u <- sum(asv_list_table_CD_u[,2:ncol(asv_list_table_CD_u)])

asv_list_table_ABC_u <- asv_table[asv_table$ASV %in% ABC_u,]
asv_list_reads_ABC_u <- sum(asv_list_table_ABC_u[,2:ncol(asv_list_table_ABC_u)])

asv_list_table_ABD_u <- asv_table[asv_table$ASV %in% ABD_u,]
asv_list_reads_ABD_u <- sum(asv_list_table_ABD_u[,2:ncol(asv_list_table_ABD_u)])

asv_list_table_ACD_u <- asv_table[asv_table$ASV %in% ACD_u,]
asv_list_reads_ACD_u <- sum(asv_list_table_ACD_u[,2:ncol(asv_list_table_ACD_u)])

asv_list_table_BCD_u <- asv_table[asv_table$ASV %in% BCD_u,]
asv_list_reads_BCD_u <- sum(asv_list_table_BCD_u[,2:ncol(asv_list_table_BCD_u)])

asv_list_table_ABCD <- asv_table[asv_table$ASV %in% ABCD,]
asv_list_reads_ABCD <- sum(asv_list_table_ABCD[,2:ncol(asv_list_table_ABCD)])

# Vector containing all asv_list_reads values
asv_list_reads_all <- c(asv_list_reads_A_u, asv_list_reads_B_u, asv_list_reads_C_u, asv_list_reads_D_u, 
                        asv_list_reads_AB_u, asv_list_reads_AC_u, asv_list_reads_AD_u, asv_list_reads_BC_u, 
                        asv_list_reads_BD_u, asv_list_reads_CD_u, asv_list_reads_ABC_u, asv_list_reads_ABD_u, 
                        asv_list_reads_ACD_u, asv_list_reads_BCD_u, asv_list_reads_ABCD)

# Total sum of all asv_list_reads values
total_asv_list_reads <- sum(asv_list_reads_all)
total_asv <- sum(rowSums(subset_data))

asv_list_reads_all_perc <- round(100 * asv_list_reads_all / sum(asv_list_reads_all), 5)

# Define names for the ASV read counts
asv_names <- c(
  "A_u", "B_u", "C_u", "D_u",
  "AB_u", "AC_u", "AD_u", "BC_u",
  "BD_u", "CD_u", "ABC_u", "ABD_u",
  "ACD_u", "BCD_u", "ABCD"
)

# Create the result data frame
result_read_percentage <- data.frame(
  Name = asv_names,
  ASV_List_Reads = asv_list_reads_all,
  ASV_List_Reads_Percentage = asv_list_reads_all_perc)
print (result_read_percentage)


# Merge the data frames vertically of water-sediment shared ASVs
merged_asv_list <- rbind(asv_list_table_BD_u,
                         asv_list_table_BCD_u,
                         asv_list_table_BC_u,
                         asv_list_table_ABCD,
                         asv_list_table_AB_u,
                         asv_list_table_ABD_u,
                         asv_list_table_ABC_u)

# export for distance to to coast in table S3
write.csv(merged_asv_list, "merged_asv_list.csv")

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
  left_join(metadata_id %>% select(sample_id, station, sample_type, location, lat, long, water_depth), by = "sample_id")

# View the updated data frame
head(long_benthic_com)
write.csv(long_benthic_com, "long_benthic_com.csv")





#import the shared ASVs and set threhold for abundance
shared_ASV <- long_benthic_com


# Initialize an empty data frame to store all results
distance_results <- data.frame()

# Define the thresholds you want to test
thresholds <- c(0, 10, 100, 1000, 10000)

# Loop through each threshold
for (abundance_threshold in thresholds) {
  cat("\nProcessing with abundance threshold:", abundance_threshold, "\n")
  
  # Update abundance values: set to 0 if less than the current threshold
  filtered_ASV <- shared_ASV %>%
    mutate(abundance = if_else(abundance < abundance_threshold, 0, abundance))
  
  # Loop through each ASV with progress notifications
  for (asv in unique(filtered_ASV$ASV)) {
    cat("\nProcessing ASV:", asv, "\n")  # Notify which ASV is being processed
    
    # Subset data for the current ASV
    asv_data <- filtered_ASV %>% filter(ASV == asv)
    
    # Filter sediment stations with the ASV present
    sediment_stations <- asv_data %>%
      filter(sample_type == "sediment" & abundance > 0) %>%
      filter(!is.na(station) & station != "")
    
    # Skip if no valid sediment stations
    if (nrow(sediment_stations) == 0) {
      cat("  - No valid sediment stations for ASV:", asv, "\n")
      next
    }
    
    # Loop through each sediment station
    for (i in 1:nrow(sediment_stations)) {
      first_station <- sediment_stations[i, ]
      
      # Check if first_station is valid
      if (nrow(first_station) == 0 || is.na(first_station$station) || first_station$station == "") {
        cat("  - Skipping invalid sediment station\n")
        next
      }
      
      cat("  - Found ASV in 'Sediment' at station", first_station$station, "\n")
      
      # Check all other stations where ASV is present in water but absent in sediment
      water_stations <- asv_data %>%
        filter(
          station != first_station$station,
          sample_type != "sediment",
          abundance > 0
        )
      
      sediment_absent_stations <- asv_data %>%
        filter(
          station %in% water_stations$station,
          sample_type == "sediment",
          abundance == 0
        )
      
      # Ensure that stations meet the criteria
      valid_stations <- inner_join(water_stations, sediment_absent_stations, by = "station") %>%
        filter(!is.na(station), !is.na(lat.x), !is.na(long.x), !is.na(sample_type.x))
      
      # Calculate distances for each valid pair
      if (nrow(valid_stations) > 0) {
        for (j in 1:nrow(valid_stations)) {
          water_station <- valid_stations[j, ]
          
          # Skip invalid stations
          if (is.na(water_station$station) || water_station$station == "" || 
              is.na(water_station$lat.x) || is.na(water_station$long.x)) {
            cat("      Skipping invalid station\n")
            next
          }
          
          cat("    -> Checking distance to station", water_station$station, "of type", water_station$sample_type.x, "\n")
          
          # Calculate distance using distHaversine
          distance <- distHaversine(
            cbind(first_station$long, first_station$lat),
            cbind(water_station$long.x, water_station$lat.x)
          )
          
          if (is.na(distance)) {
            cat("      Distance could not be calculated (missing data)\n")
            next
          }
          
          cat("      Distance calculated:", round(distance, 2), "meters\n")
          
          # Add unique results to distance_results, include the threshold as a column
          distance_results <- rbind(distance_results, data.frame(
            ASV = asv,
            station_1 = first_station$station,
            station_2 = water_station$station,
            sample_type = water_station$sample_type.x,
            d2s = distance, # Distance 2 stations
            threshold = abundance_threshold  # Add the threshold column
          ))
        }
      } else {
        cat("    -> No valid stations found for ASV:", asv, "\n")
      }
    }
  }
}


head(distance_results)

write.csv(distance_results, "distance_results.csv")

results_100 <- distance_results %>%
  filter(threshold >= 100)


distance_summary <- results_100 %>%
  group_by(ASV) %>%
  summarise(
    min_d2s = min(d2s, na.rm = TRUE),
    max_d2s = max(d2s, na.rm = TRUE),
    mean_d2s = mean(d2s, na.rm = TRUE),
    median_d2s = median(d2s, na.rm = TRUE)
  )

head(distance_summary)


# Read the taxon assignment table CSV file into a dataframe 
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)

head(taxon)


# Check the updated distance_results
head(distance_results)


# Convert row names of taxon to a column for merging
taxon_df <- taxon %>%
  rownames_to_column(var = "ASV")

# Merge taxon data with distance_summary
distance_summary_taxon <- distance_summary %>%
  left_join(taxon_df, by = "ASV")

write.csv(distance_summary_taxon, "distance_summary_taxon.csv")


# Create a mapping table from `metadata_id` with only unique station and location pairs
location_mapping <- metadata_id %>%
  select(station, location) %>%
  distinct()

# Join `location_mapping` twice to get `region_1` and `region_2` in `distance_results`
distance_results <- distance_results %>%
  # Add region_1 based on station_1
  left_join(location_mapping, by = c("station_1" = "station")) %>%
  dplyr::rename(region_1 = location) %>%
  # Add region_2 based on station_2
  left_join(location_mapping, by = c("station_2" = "station")) %>%
  dplyr::rename(region_2 = location)



# Perform the merge to add the abundance column
distance_results <- merge(distance_results, shared_ASV[, c("ASV", "station", "sample_type", "abundance")], 
                          by.x = c("ASV", "station_2", "sample_type"), 
                          by.y = c("ASV", "station", "sample_type"), 
                          all.x = TRUE)

# Display the results with distances for ASVs that match criteria
head(distance_results)

write.csv(distance_results, "distance_results_00.csv")

distance_results <- read.csv("distance_results_00.csv", header = TRUE, row.names = 1 )

head(distance_results)



# Assuming your dataframe is named 'distance_results'
min_distance_summary <- distance_results %>%
  group_by(ASV) %>%
  slice_min(order_by = d2s, n = 1) %>%
  ungroup()

min_distance_summary

write.csv(min_distance_summary,"min_distance_summary.csv")

# Figure S 

# Ensure sample_type is a factor with the specified order
distance_results$sample_type <- factor(distance_results$sample_type, 
                                     levels = c("bottom_water",
                                                "100m_water",
                                                "surface_water"))

# Ensure sample_type is a factor with the specified order
distance_results$threshold <- factor(distance_results$threshold, 
                                     levels = c(10000,
                                                1000,
                                                100,
                                                10,
                                                0))



sample_color = rev(c("#19A9A9", "#053F5C", "#F2BC52"))
taxon_color =c("#2078B4", "#F47E1F", "#AD2423", "gray")
pair_comparisons <- list( c("surface_water", "100m_water"),
                          c("surface_water", "bottom_water"),
                          c("100m_water","bottom_water"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

p1 <- ggplot(distance_results, aes(x = sample_type, y = d2s / 1000, fill = sample_type)) +
  geom_violin(aes(color = interaction(sample_type, threshold)), trim = T) +
  scale_fill_manual(values = sample_color) + 
  stat_summary(fun.data = mean_sd, geom = "pointrange", color = "red", size = 0.5) +
  stat_compare_means (method = "wilcox.test",
                      label = "p.signif", tip.length = 0.01,
                      symnum.args = symnum.args, 
                      show.legend = F,  hide.ns = F, 
                      comparisons = pair_comparisons) + 
  labs(
    title = "Distance distribution between stations for dispersal ASV",
    x = "Water Type",
    y = "Distance (km)") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.border = element_rect(linetype = "solid", fill = NA )) + 
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold

p1 


# Read the ASV table CSV file into a dataframe with row names and headers
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)

distance_results <- merge(distance_results, taxon[, c("taxon1", "taxon2")], 
                          by.x = "ASV", by.y = "row.names", all.x = TRUE)


# Calculate the total unique count of taxon1 for each threshold and sample_type
total_counts <- distance_results %>%
  distinct(ASV, sample_type, threshold, .keep_all = TRUE) %>%
  group_by(threshold, sample_type) %>%
  dplyr::summarize(total_taxon1_count = n(), .groups = "drop")


# Calculate the percentage of taxon1 within each threshold and sample_type
distance_results_percent <- distance_results %>%
  distinct(ASV, sample_type, threshold, .keep_all = TRUE) %>%
  group_by(threshold, sample_type, taxon1) %>%
  tally() %>%
  ungroup() %>%
  left_join(total_counts, by = c("threshold", "sample_type")) 

# Create a bar plot showing the percentage of taxon1 for each sample_type, faceted by threshold
p2<-ggplot(distance_results_percent, aes(x = sample_type, y = n, fill = taxon1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual (values = taxon_color) + 
  labs(
    title = "Percentage of Taxon1 for Each Water Type, Faceted by Threshold",
    x = "Water Type",
    y = "Percentage") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.border = element_rect(linetype = "solid", fill = NA )) +  
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold

p2


# Create a bar plot showing the percentage of taxon1 for each sample_type, faceted by threshold
ggplot(distance_results_percent, aes(x = sample_type, y = n, fill = taxon1)) +
  geom_bar(stat = "identity") +
  geom_text(data = total_counts, aes(x = sample_type, y = 950, label = total_taxon1_count), 
            inherit.aes = FALSE, color = "black", size = 4, vjust = -0.5) +  # Add total count as text labels
  scale_fill_manual (values = taxon_color) + 
  labs(
    title = "Percentage of Taxon1 for Each Water Type, Faceted by Threshold",
    x = "Water Type",
    y = "ASV richness") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed") # Add gridlines for outliers
  ) + 
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold


# Calculate the percentage of taxon1 within each threshold and sample_type
distance_results_abundance <- distance_results %>%
  distinct( ASV, station_2, sample_type, threshold, abundance, .keep_all = TRUE) %>%
  group_by( ASV, station_2, sample_type, threshold, taxon1, abundance) %>%
  tally() %>%
  ungroup() 

# Create the bar plot
p3 <- ggplot(distance_results_abundance, aes(x = sample_type, y = abundance, fill = taxon1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual (values = taxon_color) + 
  labs(title = "Abundance by Sample Type", 
       x = "Sample Type", 
       y = "Abundance", 
       fill = "Taxon1") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.border = element_rect(linetype = "solid", fill = NA )) + 
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold


p3




library(gridExtra)

# Fig. 5A
# Filter distance_results for Threshold = 0
threshold_100_results <- distance_results %>%
  filter(threshold == 100)


p1 <- ggplot(threshold_100_results, aes(x = sample_type, y = d2s / 1000, fill = sample_type)) +
  geom_violin(aes(color = sample_type), trim = T) +
  scale_fill_manual(values = sample_color) + 
  stat_summary(fun.data = mean_sd, geom = "pointrange", color = "red", size = 0.5) +
  stat_compare_means (method = "wilcox.test",
                      label = "p.signif", tip.length = 0.01,
                      symnum.args = symnum.args, 
                      show.legend = F,  hide.ns = F, 
                      comparisons = pair_comparisons) + 
  labs(
    title = "Distance distribution between stations for dispersal ASV",
    x = "Water Type",
    y = "Distance (km)") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.border = element_rect(linetype = "solid", fill = NA )) + 
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold
p1 

# Calculate the percentage of taxon1 within each threshold and sample_type
threshold_100_percent <- threshold_100_results %>%
  distinct(ASV, sample_type, threshold, .keep_all = TRUE) %>%
  group_by(threshold, sample_type, abundance, taxon1) %>%
  tally() %>%
  ungroup() 
 
head(threshold_100_percent)

# Calculate total abundance and percentage
threshold_100_summary <- threshold_100_percent %>%
  group_by(taxon1) %>%  # Group by taxon1
  dplyr::summarize(
    total_abundance = sum(n),  # Sum abundance for each taxon
    .groups = "drop"
  ) %>%
  mutate(
    percentage = (total_abundance / sum(total_abundance)) * 100  # Calculate percentage
  )


# Calculate total abundance and percentage for each taxon1 within each sample_type
threshold_100_summary <- threshold_100_percent %>%
  group_by(sample_type, taxon1) %>%  # Group by sample_type and taxon1
  dplyr::summarize(
    total_abundance = sum(n),  # Total abundance for each group
    .groups = "drop"  # Drop grouping for cleaner output
  ) %>%
  group_by(sample_type) %>%  # Regroup by sample_type to calculate percentage
  mutate(
    percentage = (total_abundance / sum(total_abundance)) * 100  # Calculate percentage
  ) %>%
  ungroup()

# Calculate total abundance and percentage for each taxon1 within each sample_type
threshold_100_summary <- threshold_100_percent %>%
  group_by(sample_type, taxon1) %>%  # Group by sample_type and taxon1
  dplyr::summarize(
    total_abundance = sum(n),  # Total abundance for each group
    .groups = "drop"  # Drop grouping for cleaner output
  ) %>%
  group_by(sample_type) %>%  # Regroup by sample_type to calculate percentage
  mutate(
    percentage = (total_abundance / sum(total_abundance)) * 100  # Calculate percentage
  ) %>%
  ungroup()

# Fig. 5B
# Create a bar plot showing the percentage of taxon1 for each sample_type, faceted by threshold
p2<-ggplot(threshold_100_percent, aes(x = sample_type, y = n, fill = taxon1)) +
  scale_fill_manual (values = taxon_color) + 
  geom_bar(stat = "identity", position = "fill") +  # Bar plot with percentage values
  labs(
    title = "Percentage of Taxon1 for Each Water Type, Faceted by Threshold",
    x = "Water Type",
    y = "Percentage"
  ) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.border = element_rect(linetype = "solid", fill = NA )) +  
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold
p2

# Calculate the percentage of taxon1 within each threshold and sample_type
threshold_100_abundance  <- threshold_100_results %>%
  distinct(ASV, station_2, sample_type, threshold,  abundance, .keep_all = TRUE) %>%
  group_by( ASV, station_2, sample_type, threshold, taxon1, abundance) %>%
  tally() %>%
  ungroup() 


# Calculate total abundance and percentage for each taxon1 within each sample_type
threshold_100_summary <- threshold_100_abundance %>%
  group_by(sample_type, taxon1) %>%  # Group by sample_type and taxon1
  dplyr::summarize(
    total_abundance = sum(abundance),  # Total abundance for each group
    .groups = "drop"  # Drop grouping for cleaner output
  ) %>%
  group_by(sample_type) %>%  # Regroup by sample_type to calculate percentage
   mutate(
    percentage = (total_abundance / sum(total_abundance)) * 100  # Calculate percentage
  ) %>%
  ungroup()



# Calculate total abundance and percentage
threshold_100_summary <- threshold_100_abundance %>%
  group_by(taxon1) %>%  # Group by taxon1
  dplyr::summarize(
    total_abundance = sum(abundance),  # Sum abundance for each taxon
    .groups = "drop"
  ) %>%
  mutate(
    percentage = (total_abundance / sum(total_abundance)) * 100  # Calculate percentage
  )


# Fig. 5C
# Create a bar plot showing the percentage of taxon1 for each sample_type, faceted by threshold
p3<-ggplot(threshold_100_abundance, aes(x = sample_type, y = abundance, fill = taxon1)) +
  scale_fill_manual (values = taxon_color) + 
  geom_bar(stat = "identity", position ="fill") +  # Bar plot with percentage values
  labs(
    title = "Percentage of Taxon1 for Each Water Type, Faceted by Threshold",
    x = "Water Type",
    y = "Percentage"
  ) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), # Add gridlines for outliers
    panel.border = element_rect(linetype = "solid", fill = NA )) + 
  facet_wrap(~ threshold, ncol =1)  # Facet by Threshold
p3

## Fig. 5A, B, C
grid.arrange(p1, p2, p3, ncol = 3)




lists1 = rev(c("NS06",
           "NS02",
           "NS04",
           "NS01",
           "NS10",
           "NS47",
           "NS12",
           "NS45",
           "NS43",
           "NS38",
           "HR07",
           "HR04",
           "IS01",
           "KV01",
           "KH01",
           "IS02",
           "KH10",
           "NC01",
           "KH05",
           "KH02",
           "KH03",
           "KH09",
           "KH08",
           "KH07"))




# Count unique ASVs for each combination of station_1 and region_2
unique_asvs <- threshold_100_results %>%
  distinct(ASV, sample_type, .keep_all = TRUE) %>%
  tidyr::complete(station_1, region_2, sample_type, fill = list(ASV = NA)) %>%  # Fill missing combinations
  group_by(station_1,  region_2, sample_type) %>%  # Group by both station_1 and sample_type
  dplyr::summarize(
    unique_ASV_count = n_distinct(ASV, na.rm = TRUE),  # Count unique ASVs, ignoring NAs
    .groups = "drop"  # Drop the grouping after summarizing
  )


# Summarize unique ASV count for station_1 based on region_2
asv_station_water <- unique_asvs %>%
  group_by(station_1, sample_type) %>%  # Group by station_1 and region_2
  dplyr::summarize(
    total_ASV_count = sum(unique_ASV_count, na.rm = TRUE),  # Sum ASV counts within each group
    .groups = "drop"  # Drop grouping after summarization
  )


# Create the heatmap using ggplot2

p1 <- ggplot(asv_station_water, aes(x = factor(station_1, levels = lists1), y = sample_type, fill = total_ASV_count)) +
  geom_tile(color = "white", size = 1) +  # Create tiles with white borders
   geom_text(aes(label = total_ASV_count), color = "black", size = 4) +  # Add ASV counts as text
  scale_fill_gradientn(
    colors = c("white", "#90e0ef", "#0077b6", "#03045e"),  # Include gray for 0
    values = scales::rescale(c(0, 1, 20, 40)),          # Map 0 to gray explicitly
    limits = c(0, 40),                                  # Adjust scale limits
    breaks = seq(0, 40, by = 10)                        # Breaks for the legend
  ) +
  labs(
    title = "Dispersal ASV Count for Each Station to others by Water type",
    x = "Station",
    y = "Water type",
    fill = "Unique ASVs"
  ) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = -90),  # Rotate x-axis labels
    axis.text.y = element_text(angle = 0),  # Rotate y-axis labels (optional)
    panel.grid.major = element_line(color = "white"),  # Add major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

p1


# Summarize unique ASV count for station_1 based on region_2
asv_station_region <- unique_asvs %>%
  group_by(station_1, region_2) %>%  # Group by station_1 and region_2
  dplyr::summarize(
    total_ASV_count = sum(unique_ASV_count, na.rm = TRUE),  # Sum ASV counts within each group
    .groups = "drop"  # Drop grouping after summarization
  )

# Calculate the sum of total_ASV_count
sum_total_ASV_count <- sum(asv_station_region$total_ASV_count, na.rm = TRUE)

# View the sum
sum_total_ASV_count



# dispersal to region 2
lists2 = rev(c("North Svalbard",
               "West Svalbard",
               "North Norway",
               "NE Greenland",
               "Greenland Ridge"))


# Create the heatmap using ggplot2

p2 <- ggplot(asv_station_region, aes(x = factor(station_1, levels = lists1), 
                               y = factor(region_2, levels = lists2),
                               fill = total_ASV_count)) +
  geom_tile(color = "white", size = 1) +  # Create tiles with white borders
  geom_text(aes(label = total_ASV_count), color = "black", size = 4) +  # Add ASV counts as text
  scale_fill_gradientn(
    colors = c("white", "#90e0ef", "#0077b6", "#03045e"),  # Include gray for 0
    values = scales::rescale(c(0, 1, 20, 40)),          # Map 0 to gray explicitly
    limits = c(0, 40),                                  # Adjust scale limits
    breaks = seq(0, 40, by = 10)                        # Breaks for the legend
  ) +
  labs(
    title = "Dispersal ASV Count for Each Station to others region",
    x = "Station",
    y = "Water type",
    fill = "Unique ASVs"
  ) +
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = -90),  # Rotate x-axis labels
    axis.text.y = element_text(angle = 0),  # Rotate y-axis labels (optional)
    panel.grid.major = element_line(color = "white"),  # Add major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

p2


grid.arrange(p1, p2, ncol = 1)



#check the where the ASVs from other stations via water layers and regions
# Count unique ASVs for each combination of station_1 and region_2
unique_asvs <- threshold_100_results %>%
  distinct(ASV, sample_type, .keep_all = TRUE) %>%
  tidyr::complete(station_2, region_1, sample_type, fill = list(ASV = NA)) %>%  # Fill missing combinations
  group_by(station_2,  region_1, sample_type) %>%  # Group by both station_1 and sample_type
  dplyr::summarize(
    unique_ASV_count = n_distinct(ASV, na.rm = TRUE),  # Count unique ASVs, ignoring NAs
    .groups = "drop"  # Drop the grouping after summarizing
  )



# Summarize unique ASV count for station_1 based on region_2
asv_station_water <- unique_asvs %>%
  group_by(station_2, sample_type) %>%  # Group by station_1 and region_2
  dplyr::summarize(
    total_ASV_count = sum(unique_ASV_count, na.rm = TRUE),  # Sum ASV counts within each group
    .groups = "drop"  # Drop grouping after summarization
  )


# Create the heatmap using ggplot2

p1 <- ggplot(asv_station_water, aes(x = factor(station_2, levels = lists1), y = sample_type, fill = total_ASV_count)) +
  geom_tile(color = "white", size = 1) +  # Create tiles with white borders
  geom_text(aes(label = total_ASV_count), color = "black", size = 4) +  # Add ASV counts as text
  scale_fill_gradientn(
    colors = c("white", "#90e0ef", "#0077b6", "#03045e"),  # Include gray for 0
    values = scales::rescale(c(0, 1, 25, 50)),          # Map 0 to gray explicitly
    limits = c(0, 50),                                  # Adjust scale limits
    breaks = seq(0, 50, by = 10)                        # Breaks for the legend
  ) +
  labs(
    title = "Dispersal ASV Count for Each Station to others by Water type",
    x = "Station",
    y = "Water type",
    fill = "Unique ASVs"
  ) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = -90),  # Rotate x-axis labels
    axis.text.y = element_text(angle = 0),  # Rotate y-axis labels (optional)
    panel.grid.major = element_line(color = "white"),  # Add major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

p1


# Summarize unique ASV count for station_1 based on region_2
asv_station_region <- unique_asvs %>%
  group_by(station_2, region_1) %>%  # Group by station_1 and region_2
  dplyr::summarize(
    total_ASV_count = sum(unique_ASV_count, na.rm = TRUE),  # Sum ASV counts within each group
    .groups = "drop"  # Drop grouping after summarization
  )

# Calculate the sum of total_ASV_count
sum_total_ASV_count <- sum(asv_station_region$total_ASV_count, na.rm = TRUE)

# View the sum
sum_total_ASV_count



# dispersal to region 2
lists2 = rev(c("North Svalbard",
               "West Svalbard",
               "North Norway",
               "NE Greenland",
               "Greenland Ridge"))


# Create the heatmap using ggplot2

p2 <- ggplot(asv_station_region, aes(x = factor(station_2, levels = lists1), 
                                     y = factor(region_1, levels = lists2),
                                     fill = total_ASV_count)) +
  geom_tile(color = "white", size = 1) +  # Create tiles with white borders
  geom_text(aes(label = total_ASV_count), color = "black", size = 4) +  # Add ASV counts as text
  scale_fill_gradientn(
    colors = c("white", "#90e0ef", "#0077b6", "#03045e"),  # Include gray for 0
    values = scales::rescale(c(0, 1, 25, 50)),          # Map 0 to gray explicitly
    limits = c(0, 50),                                  # Adjust scale limits
    breaks = seq(0, 50, by = 10)                        # Breaks for the legend
  ) +
  labs(
    title = "Dispersal ASV Count for Each Station from others by Water type",
    x = "Station",
    y = "Water type",
    fill = "Unique ASVs"
  ) +
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = -90),  # Rotate x-axis labels
    axis.text.y = element_text(angle = 0),  # Rotate y-axis labels (optional)
    panel.grid.major = element_line(color = "white"),  # Add major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

p2


grid.arrange(p1, p2, ncol = 1)








