# Load required libraries
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(ggrepel)
library(dplyr)
library(webr)

# Read the ASV table CSV file into a dataframe
ASV_table <- read.csv("ASV_table_37F.csv", row.names = 1, header = TRUE)
com <- as.matrix(ASV_table)

# Read the metadata CSV file into a dataframe 
metadata <- read.csv("nordic_metadata.csv", row.names = 1, header = TRUE)


#Fig. 2A Pie_chart
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

# Convert totals_by_ASV to numeric and preserve row names
pie_ASV <- totals_by_ASV %>%
  apply(2, as.numeric) %>%
  as.data.frame()

# Keep row names and reorder columns
rownames(pie_ASV) <- rownames(totals_by_ASV)
pie_ASV <- pie_ASV[, c("surface_water", "100m_water", "bottom_water", "sediment")]


#Pie-Donut chart based on read abundance
# Add group of foram information to ASV table
group_foram <- taxon %>% rownames_to_column("ASV_id")

# Merge ASV abundance data with taxon group information
pie_taxon <- pie_ASV %>%
  rownames_to_column("ASV_id") %>%
  left_join(group_foram %>% select(ASV_id, group), by = "ASV_id") %>%
  select(ASV_id, group, everything()) %>%
  column_to_rownames("ASV_id")


# Group by the sample type and calculate the sum of each ASV
totals_by_abun <- pie_taxon %>%
  group_by(group) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  rownames_to_column("sample_type") %>%
  pivot_longer(cols = -sample_type, names_to = "group", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>%
  group_by(sample_type, group) %>%
  summarise(n = sum(n))

# Reorder factors for plotting
totals_by_abun$sample_type <- factor(totals_by_abun$sample_type, levels=c("surface_water",  "100m_water", "bottom_water", "sediment"))
totals_by_abun$group <- factor(totals_by_abun$group, levels=c("planktonic",  "benthic", "unassigned"))


# Create Pie-Donut chart
p1 <- PieDonut(totals_by_abun , aes( sample_type, group, count = n), r0 = 0.45, r1 = 0.9)


#Pie-Donut chart based on ASV Richness
# Modify the dataframe so values > 0 are set to 1 and values = 0 remain 0
totals_by_rich  <- pie_taxon %>%
  mutate(across(surface_water:sediment, ~ ifelse(. > 0, 1, 0))) %>%
  group_by(group) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  rownames_to_column("sample_type") %>%
  pivot_longer(cols = -sample_type, names_to = "group", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>%
  group_by(sample_type, group) %>%
  summarise(n = sum(n))

# Reorder factors for plotting
totals_by_rich$sample_type <- factor(totals_by_rich$sample_type, levels=c("surface_water",  "100m_water", "bottom_water", "sediment"))
totals_by_rich$group <- factor(totals_by_rich$group, levels=c("planktonic",  "benthic", "unassigned"))


# Create Pie-Donut chart
p2 <- PieDonut(totals_by_rich , aes( sample_type, group, count = n), r0 = 0.45, r1 = 0.9)


#Fig. 2B Venn diagram
# Load required packages
library(VennDiagram)

# Convert totals_by_ASV to numeric, excluding the row names
venn_ASV <- as.data.frame(apply(totals_by_ASV, 2, as.numeric))

# Keep the row names andReorder the columns
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
               "Bottom water \n (BW)"),
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
result_read_percenta <- data.frame(
  Name = asv_names,
  ASV_List_Reads = asv_list_reads_all,
  ASV_List_Reads_Percentage = asv_list_reads_all_perc)
print (result_read_percenta)


# Fig. 2C Accumulation curve
# Filter rows for "surface_water" sample type and remove columns with sum equal to 0
SW <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type == "surface_water"))) %>%
  select_if(~ sum(.) != 0)

# Filter rows for "100m_water" sample type and remove columns with sum equal to 0
M100 <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type == "100m_water"))) %>%
  select_if(~ sum(.) != 0)

# Filter rows for "bottom_water" sample type and remove columns with sum equal to 0
BW <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type == "bottom_water"))) %>%
  select_if(~ sum(.) != 0)

# Filter rows for "sediment" sample type and remove columns with sum equal to 0
SED <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type == "sediment"))) %>%
  select_if(~ sum(.) != 0)


# List of dataframes and colors
curves <- list(SW, M100, BW, SED)

# Colors for curves and confidence intervals
sample_color <- c("#19A9A9", "#053F5C", "#F2BC52", "#F68E51")
ci_colors <- c(rgb(0, 0.8, 1, 0.2), rgb(0, 0.2, 0.1, 0.2), rgb(1, 1, 0.2, 0.2), rgb(1, 0.8, 0.4, 0.2))

# Plot each curve
for (i in seq_along(curves)) {
  curve <- specaccum(curves[[i]], method = "random", permutations = 1000)
  plot(curve, ci.type = "poly", col = sample_color[i], lwd = 3, lty = 1,
       ci.lty = 0, ci.col = ci_colors[i], ylim = c(0, 3800), xlim = c(1, 25),
       xlab = "Number of samples", ylab = "Number of ASVs", add = i > 1)
}

# Legend
sample_type <- c("Surface water", "100m water", "Bottom water", "Sediment")
legend(x = 1, y = 3700, legend = sample_type, 
       fill = sample_color, box.lty = 0)


# Fig. 2D Shannon
# Calculate the Shannon diversity index for each sample
shannon <- diversity(com, index = "shannon")

# Add the Shannon index to the metadata and set sample type as a factor
metadata$shannon <- shannon
metadata$sample_type <- factor(metadata$sample_type, levels=c("surface_water",  "100m_water", "bottom_water", "sediment"))

# Perform pairwise Wilcoxon test between sample types
pairwise.wilcox.test(metadata$shannon, metadata$sample_type,
                     p.adjust.method = "none")

# Create a violin plot of the Shannon diversity index by sample type
pair_comparisons <- list( c("surface_water", "sediment"),
                          c("surface_water", "bottom_water"),
                          c("100m_water","sediment"),
                          c("100m_water","bottom_water"),
                          c("bottom_water","sediment"))

# Define significance symbols for the Wilcoxon test results
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                   symbols = c("****", "***", "**", "*", "ns"))

# Sample type labels for the plot
sample_type <- c("SW", "100mW", "BW", "SED")

# Create a violin plot of the Shannon diversity index by sample type
ggplot(metadata, aes(x = sample_type, y = shannon), color = sample_type) +
  geom_boxplot(aes(color = sample_type), alpha = 0.5) +
  geom_jitter(aes(color = sample_type), size= 2.5, alpha=0.9) +
  scale_colour_manual(values = sample_color) +
  labs(x = "Sample Type", y = "Shannon diversity index") +
  stat_compare_means (method = "wilcox.test",
                      label = "p.signif", tip.length = 0.01,
                      symnum.args = symnum.args, 
                      show.legend = F,  hide.ns = TRUE, 
                      comparisons = pair_comparisons) + 
  annotate("text", x = 3.5, y = 0, label = "Wilcoxon tests, ***P < 0.001 and ****P < 0.0001", size = 4) +
  theme_classic() + guides(color = FALSE) +
  theme (axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
         axis.text.x = element_text(colour = "black", face = "bold", size = 10),
         axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
         axis.title.y = element_text(face = "bold", size = 10, colour = "black"))+
  scale_x_discrete(labels = sample_type)


