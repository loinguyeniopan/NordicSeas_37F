# Load required libraries
library(vegan)
library(metagenomeSeq)
library(devtools)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(gridExtra)

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

# Filter the data based on the group and extract benthic ASVs
benthic_ASV <- rownames(taxon[taxon$group %in% c("benthic", "unassigned"),])
benthic_com <- com[rownames(com) %in% benthic_ASV , ]
benthic_com <- t(benthic_com)

# Remove rows with sum less than 100
N <- rowSums(benthic_com)
benthic_com <- benthic_com[N >= 100, ]


#Extract the corresponding meta_table for the samples in abund_table
metadata <- metadata[rownames(benthic_com),]

#Import for CSS normalization
OTU_read_count <- as.matrix(t(benthic_com))


# Normalize dataset by CSS
# convert OTU table into package format
metaSeqObject = newMRexperiment(OTU_read_count) 

# CSS normalization
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )

# Convert CSS normalized data into data.frame-formatted OTU table (log transformed data)
OTU_read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))

# Turn abundance data frame into a matrix
m_com = as.matrix(t(OTU_read_count_CSS))


# Permutational multivariate analysis of variance
set.seed(123)
ado = adonis2(m_com ~ time*sample_type*location*long*distance*water_depth,
              data=metadata, method = "bray", permutations=999, na.action = na.omit)
ado
summary(ado)

# Mantel test
library(vegan)
library(geosphere)

commdist<-vegdist(m_com, method="bray")

# Extract the rows from `env` that have the same row names as in `com`
metadata <- metadata[rownames(metadata) %in% rownames(m_com), ]
geo <- metadata[ ,7:8]
coast = metadata$distance
env <- metadata[ ,10:16]


# Geographic data frame - haversine distance 
geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(geo)
dist.coast = dist(coast, method = "euclidean")

# Scale data 
scale.env = scale(env, center = TRUE, scale = TRUE)

# Create distance matrix of scaled data
dist.env = dist(scale.env, method = "euclidean")

# Mantel test for geography
mantel(commdist, dist.geo, method="pearson", permutations=999)

# Mantel test for distance from the coast
mantel(commdist, dist.coast, method="pearson", permutations=999)

# Mantel test for selected environmental parameters
mantel(commdist, dist.env, method="pearson", permutations=999)

aa = as.vector(commdist)
gg = as.vector(dist.geo)
dd = as.vector(dist.coast)

#Fig. S3
# New data frame with vectorized distance matrices
mat = data.frame(aa,gg,dd)

mm = ggplot(mat, aes(y = aa, x = gg/1000)) + 
  geom_point(size = 2, alpha = 0.75, colour = "black",shape = 21, aes(fill = dd/1000)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.5) + 
  labs(x = " Geographic separation (Km)", y = "Bray-Curtis Dissimilarity", fill = "Offshore distance separation (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
   scale_fill_continuous(high = "navy", low = "skyblue")

mm


# Fig. 4 inset
# Calculates beta diversity dispersion within each sample type, when comparing 2 or more sample types

# Calculate Bray-Curtis dissimilarity matrix
m_com.bray <- vegdist(m_com, "bray", na.rm = TRUE)

# Convert 'sample_type' column to a factor with custom order
metadata$sample_type <- factor(metadata$sample_type, levels = c("surface_water",
                                                                "100m_water",
                                                                "bottom_water",
                                                                "sediment"))

# Perform beta diversity dispersion analysis
disp <- betadisper(m_com.bray, metadata$sample_type, type = "centroid")

# Test if centroid distances are significantly different from each other using ANOVA
anova_result <- anova(disp)
print(anova_result)

# Conduct a permutation-based test
adonis_result <- adonis2(dist(disp$distances) ~ metadata$sample_type)
print(adonis_result)

# Test significance between each group using Tukey's Honest Significant Difference test
disp.TukeyHSD <- TukeyHSD(disp)
print(disp.TukeyHSD)

# Plot beta diversity dispersion with ellipses
plot(disp, hull = FALSE, ellipse = TRUE)

# Perform pairwise Wilcoxon tests
pairwise_wilcox_result <- pairwise.wilcox.test(disp$distances, metadata$sample_type,
                                               p.adjust.method = "none")
print(pairwise_wilcox_result)

# Define comparisons for violin plot
pair_comparisons <- list(c("surface_water", "100m_water"),
                         c("surface_water", "sediment"),
                         c("surface_water", "bottom_water"),
                         c("100m_water", "bottom_water"),
                         c("100m_water", "sediment"),
                         c("bottom_water", "sediment"))

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
sample_color = c("#19A9A9", "#053F5C", "#F2BC52", "#F68E51")


# Add distances to metadata
metadata$distances <- disp$distances

# Create a violin plot of the distances to the group centroid by sample type
ggplot(metadata, aes(x = sample_type, y = distances, fill = sample_type)) +
  geom_violin(trim = FALSE, color = "black") +
  scale_fill_manual(values = sample_color) +
  stat_summary(fun.data = mean_sd, geom = "pointrange", color = "red", size = 1) +
  labs(x = "Sample Type", y = "Distances to the Group Centroid") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", tip.length = 0.01, 
                     symnum.args = symnum.args, 
                     show.legend = FALSE, 
                     comparisons = pair_comparisons) +
  annotate("text", x = 2, y = 0.4, label = "Wilcoxon tests, *P < 0.05, ****P < 0.0001", size = 4) +
  theme_gray() + 
  guides(color = FALSE) +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 10, colour = "black"))


# Plot Fig. 4 nMDS
set.seed(123)

# Perform nMDS using Bray-Curtis distance
nmds = metaMDS(m_com, autotransform = FALSE, distance = "bray")
nmds

plot(nmds, display = c("sites", "text"))  # Set type = "n" to prevent automatic plotting
text(nmds, display = "sites", col = "red")  # Adjust color and other parameters as needed


# Select relevant environmental variables for fitting
selected_env = metadata[,7:12]

# Fit environmental vectors to nMDS ordination
envfit_result = envfit(nmds, selected_env, permutations = 999, na.rm = TRUE)
envfit_result

# Extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds)$sites)

# Add metadata to the NMDS scores data frame
data.scores$name = metadata$name
names(data.scores)[c(1, 2)] <- c("x", "y")

data.scores$sample_type = metadata$sample_type
data.scores$location = metadata$location
data.scores$water_depth = metadata$water_depth
head(data.scores)

# Convert 'sample_type' column to a factor with a custom order
data.scores$sample_type <- factor(data.scores$sample_type, levels = c("surface_water",
                                                                      "100m_water",
                                                                      "bottom_water",
                                                                      "sediment"))

# Convert 'location' column to a factor with a custom order
data.scores$location <- factor(data.scores$location, levels = c("North Norway",
                                                                "West Svalbard",
                                                                "North Svalbard",
                                                                "NE Greenland",
                                                                "Greenland Ridge"))

head(data.scores)



# Extract environmental vector coordinates and scale appropriately
en_coord_cont = as.data.frame(scores(envfit_result, "vectors")) * ordiArrowMul(envfit_result)

# Fit a smooth surface to the ordination based on water depth
ordi <- ordisurf(nmds, metadata$water_depth, plot = FALSE, bs = "ds")

# Extract the ordisurf object grid
ordi.grid <- ordi$grid
str(ordi.grid) # Check structure of the grid object

# Convert the grid object to a format suitable for plotting
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
ordi.mite$z <- as.vector(ordi.grid$z) # Unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) # Remove NA values

# Create the plot
p <- ggplot() +
  geom_point(data = data.scores, aes(x, y, color = sample_type, fill = sample_type, shape = location), size = 4) +
  scale_shape_manual(values = c(22, 24, 25, 21, 23)) +
  scale_color_manual(values = sample_color) +
  scale_fill_manual(values = sample_color) +
  labs(colour = "Sample layers", shape = "Location", x = "NMDS1", y = "NMDS2") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face = "bold", colour = "black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  geom_text_repel() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = en_coord_cont, size = 1, colour = "red") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "red", 
            fontface = "bold", label = row.names(en_coord_cont)) +
  guides(fill = FALSE)  # Remove legend for fill

p



# Fig. S4
# Filter rows for "surface_water" sample type and remove columns with sum equal to 0
SW <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type2 == "surface_water"))) %>%
  select_if(~ sum(.) != 0)

# Filter rows for "100m_water" sample type and remove columns with sum equal to 0
M100 <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type2 == "100m_water"))) %>%
  select_if(~ sum(.) != 0)

# Filter rows for "bottom_water" sample type and remove columns with sum equal to 0
BW <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type2 == "bottom_water"))) %>%
  select_if(~ sum(.) != 0)

# Filter rows for "sediment" sample type and remove columns with sum equal to 0
SED <- ASV_table %>%
  filter(row.names(.) %in% rownames(subset(metadata, sample_type2 == "sediment"))) %>%
  select_if(~ sum(.) != 0)


# Process 'SW' data

# Environmental fitting
env_vars <- c("lat", "salinity", "temperature", "long", "distance")

metadata_SW <- metadata %>%
  filter(row.names(metadata) %in% row.names(SW))


OTU_read_count <- as.data.frame(t(SW))
metaSeqObject <- newMRexperiment(OTU_read_count)

metaSeqObject_CSS <- cumNorm(metaSeqObject, p = cumNormStatFast(metaSeqObject))
OTU_read_count_CSS <- data.frame(MRcounts(metaSeqObject_CSS, norm = TRUE, log = TRUE))

m_com <- as.matrix(t(OTU_read_count_CSS))

# nMDS Analysis
set.seed(123)
nmds <- metaMDS(m_com, autotransform = FALSE, distance = "bray")

# Extract NMDS scores
data.scores <- as.data.frame(scores(nmds)$sites)
names(data.scores) <- c("x", "y")
data.scores$location <- factor(metadata_SW$location, levels = c("North Norway", "West Svalbard", "North Svalbard", "NE Greenland", "Greenland Ridge"))


ordi_results <- lapply(env_vars, function(var) ordisurf(nmds, metadata_SW[[var]], plot = FALSE, bs = "ds"))

ordi_grid_list <- lapply(ordi_results, function(ordi) {
  ordi.grid <- ordi$grid
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
  ordi.mite$z <- as.vector(ordi.grid$z)
  na.omit(ordi.mite)
})

region_color <- c("#0196C1", "#B964A6", "#d62728", "#66a61e", "#FF7F00")

#ANOSIM
ano = anosim(m_com, metadata_SW$time, distance = "bray", permutations = 9999)
ano

ano = anosim(m_com, metadata_SW$location, distance = "bray", permutations = 9999)
ano


# Process 'M100' data
metadata_M100 <- metadata %>%
  filter(row.names(metadata) %in% row.names(M100))

OTU_read_count_M100 <- as.data.frame(t(M100))
metaSeqObject_M100 <- newMRexperiment(OTU_read_count_M100)

metaSeqObject_CSS_M100 <- cumNorm(metaSeqObject_M100, p = cumNormStatFast(metaSeqObject_M100))
OTU_read_count_CSS_M100 <- data.frame(MRcounts(metaSeqObject_CSS_M100, norm = TRUE, log = TRUE))

m_com_M100 <- as.matrix(t(OTU_read_count_CSS_M100))

# nMDS Analysis for M100
set.seed(123)
nmds_M100 <- metaMDS(m_com_M100, autotransform = FALSE, distance = "bray")

# Extract NMDS scores for M100
data.scores_M100 <- as.data.frame(scores(nmds_M100)$sites)
names(data.scores_M100) <- c("x", "y")
data.scores_M100$location <- factor(metadata_M100$location, levels = c("North Norway", "West Svalbard", "North Svalbard", "NE Greenland", "Greenland Ridge"))

# Environmental fitting for M100
ordi_results_M100 <- lapply(env_vars, function(var) ordisurf(nmds_M100, metadata_M100[[var]], plot = FALSE, bs = "ds"))

ordi_grid_list_M100 <- lapply(ordi_results_M100, function(ordi) {
  ordi.grid <- ordi$grid
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
  ordi.mite$z <- as.vector(ordi.grid$z)
  na.omit(ordi.mite)
})


#ANOSIM
ano = anosim(m_com_M100, metadata_M100$time, distance = "bray", permutations = 9999)
ano

ano = anosim(m_com_M100, metadata_M100$location, distance = "bray", permutations = 9999)
ano


# Process 'BW' data
metadata_BW <- metadata %>%
  filter(row.names(metadata) %in% row.names(BW))

OTU_read_count_BW <- as.data.frame(t(BW))
metaSeqObject_BW <- newMRexperiment(OTU_read_count_BW)

metaSeqObject_CSS_BW <- cumNorm(metaSeqObject_BW, p = cumNormStatFast(metaSeqObject_BW))
OTU_read_count_CSS_BW <- data.frame(MRcounts(metaSeqObject_CSS_BW, norm = TRUE, log = TRUE))

m_com_BW <- as.matrix(t(OTU_read_count_CSS_BW))

# nMDS Analysis for BW
set.seed(123)
nmds_BW <- metaMDS(m_com_BW, autotransform = FALSE, distance = "bray")

# Extract NMDS scores for BW
data.scores_BW <- as.data.frame(scores(nmds_BW)$sites)
names(data.scores_BW) <- c("x", "y")
data.scores_BW$location <- factor(metadata_BW$location, levels = c("North Norway", "West Svalbard", "North Svalbard", "NE Greenland", "Greenland Ridge"))

# Environmental fitting for BW
ordi_results_BW <- lapply(env_vars, function(var) ordisurf(nmds_BW, metadata_BW[[var]], plot = FALSE, bs = "ds"))

ordi_grid_list_BW <- lapply(ordi_results_BW, function(ordi) {
  ordi.grid <- ordi$grid
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
  ordi.mite$z <- as.vector(ordi.grid$z)
  na.omit(ordi.mite)
})


#ANOSIM
ano = anosim(m_com_BW, metadata_BW$time, distance = "bray", permutations = 9999)
ano

ano = anosim(m_com_BW, metadata_BW$location, distance = "bray", permutations = 9999)
ano


# Process 'SED' data
metadata_SED <- metadata %>%
  filter(row.names(metadata) %in% row.names(SED))

OTU_read_count_SED <- as.data.frame(t(SED))
metaSeqObject_SED <- newMRexperiment(OTU_read_count_SED)

metaSeqObject_CSS_SED <- cumNorm(metaSeqObject_SED, p = cumNormStatFast(metaSeqObject_SED))
OTU_read_count_CSS_SED <- data.frame(MRcounts(metaSeqObject_CSS_SED, norm = TRUE, log = TRUE))

m_com_SED <- as.matrix(t(OTU_read_count_CSS_SED))

# nMDS Analysis for SED
set.seed(123)
nmds_SED <- metaMDS(m_com_SED, autotransform = FALSE, distance = "bray")

# Extract NMDS scores for SED
data.scores_SED <- as.data.frame(scores(nmds_SED)$sites)
names(data.scores_SED) <- c("x", "y")
data.scores_SED$location <- factor(metadata_SED$location, levels = c("North Norway", "West Svalbard", "North Svalbard", "NE Greenland", "Greenland Ridge"))

# Environmental fitting for SED
ordi_results_SED <- lapply(env_vars, function(var) ordisurf(nmds_SED, metadata_SED[[var]], plot = FALSE, bs = "ds"))

ordi_grid_list_SED <- lapply(ordi_results_SED, function(ordi) {
  ordi.grid <- ordi$grid
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
  ordi.mite$z <- as.vector(ordi.grid$z)
  na.omit(ordi.mite)
})


#ANOSIM
ano = anosim(m_com_SED, metadata_SED$time, distance = "bray", permutations = 9999)
ano

ano = anosim(m_com_SED, metadata_SED$location, distance = "bray", permutations = 9999)
ano



# Plot Fig. S4 nMDS

plot_nmds <- function(data.scores, ordi_mite_na, contour_col, title_label) {
  ggplot() +
    geom_point(data = data.scores, aes(x, y, fill = location, shape = location), size = 4) +
    scale_fill_manual(values = region_color, name = "Location") +  # Specify legend name and values for fill
    scale_shape_manual(values = c(22, 24, 25, 21, 23), name = "Location") +  # Specify legend name and values for shape
    stat_contour(data = ordi_mite_na, aes(x = x, y = y, z = z, colour = "contour1"), size = 0.5, position = "identity") +
    scale_colour_manual(name = "Contour", values = contour_col, labels = title_label) +  # Specify legend name and values for colour
    labs(color = "Contour", x = "NMDS1", y = "NMDS2") +  # Update color legend title
    theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10),
          legend.text = element_text(size = 10, face = "bold", colour = "black"),
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key = element_blank()) 
}

# Define function to plot nMDS with contours
plot_nmds_combined <- function(data_scores, ordi_mite_na2, ordi_mite_na3, contour_col2, contour_col3, title_label2, title_label3) {
  ggplot() +
    geom_point(data = data_scores, aes(x, y, fill = location, shape = location), size = 4) +
    scale_shape_manual(values = c(22, 24, 25, 21, 23),name = "Location") +
    scale_fill_manual(values = region_color, name = "Location") +
    stat_contour(data = ordi_mite_na2, aes(x = x, y = y, z = z, colour = "contour2"), size = 0.5, position = "identity") +
    stat_contour(data = ordi_mite_na3, aes(x = x, y = y, z = z, colour = "contour3"), size = 0.5, position = "identity") +
    scale_colour_manual(name = " ", values = c("contour2" = contour_col2, "contour3" = contour_col3),
                        labels = c(title_label2, title_label3)) +
    labs(color = " ", x = "NMDS1", y = "NMDS2") +
    theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10),
          legend.text = element_text(size = 10, face = "bold", colour = "black"),
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key = element_blank()) 
  
}

# Plot combined contours for ordi_grid_list 2 (sal) and 3 (tem)

p1 <- plot_nmds(data.scores, ordi_grid_list[[1]], "#0000FF", "Latitude")

# Plot combined contours for ordi_grid_list 2 (sal) and 3 (tem)
p2 <- plot_nmds_combined(data.scores, ordi_grid_list[[2]], ordi_grid_list[[3]], "darkblue", "red", "Salinity", "Temperature")

# Plot combined contours for ordi_grid_list 4 (long) and 5 (distance)
p3 <- plot_nmds_combined(data.scores, ordi_grid_list[[4]], ordi_grid_list[[5]],  "#f45d01", "#62b6cb", "Longitude", "Offshore Distacne")

# Arrange plots in a grid
grid.arrange(p1, p2, p3, ncol = 1)

# Plot nMDS for M100
p4 <- plot_nmds(data.scores_M100, ordi_grid_list_M100[[1]], "#0000FF", "Latitude")

# Plot combined contours for ordi_grid_list 2 (sal) and 3 (tem)
p5 <- plot_nmds_combined(data.scores_M100, ordi_grid_list_M100[[2]], ordi_grid_list_M100[[3]],"darkblue", "red",  "Salinity", "Temperature")

# Plot combined contours for ordi_grid_list 4 (long) and 5 (distance)
p6 <- plot_nmds_combined(data.scores_M100, ordi_grid_list_M100[[4]],ordi_grid_list_M100[[5]], "#f45d01", "#62b6cb",  "Longitude", "Offshore Distacne")

# Arrange plots in a grid
grid.arrange(p4, p5, p6, ncol = 1)


# Plot nMDS for BW
p7 <- plot_nmds(data.scores_BW, ordi_grid_list_BW[[1]], "#0000FF", "Latitude")

# Plot combined contours for ordi_grid_list 2 (sal) and 3 (tem)
p8 <- plot_nmds_combined(data.scores_BW, ordi_grid_list_BW[[2]], ordi_grid_list_BW[[3]], "darkblue", "red",  "Salinity", "Temperature")

# Plot combined contours for ordi_grid_list 4 (long) and 5 (distance)
p9 <- plot_nmds_combined(data.scores_BW, ordi_grid_list_BW[[4]],ordi_grid_list_BW[[5]], "#f45d01", "#62b6cb",  "Longitude", "Offshore Distacne")

# Arrange plots in a grid
grid.arrange(p7, p8, p9, ncol = 1)


# Plot nMDS for SED
p10 <- plot_nmds(data.scores_SED, ordi_grid_list_SED[[1]], "#0000FF", "Latitude")

# Plot combined contours for ordi_grid_list 2 (sal) and 3 (tem)
p11 <- plot_nmds_combined(data.scores_SED, ordi_grid_list_SED[[2]], ordi_grid_list_SED[[3]], "darkblue", "red",  "Salinity", "Temperature")

# Plot combined contours for ordi_grid_list 4 (long) and 5 (distance)
p12 <- plot_nmds_combined(data.scores_SED, ordi_grid_list_SED[[4]],ordi_grid_list_SED[[5]], "#f45d01", "#62b6cb", "Longitude", "Offshore Distacne")

# Arrange plots in a grid
grid.arrange(p10, p11, p12, ncol = 1)


######################
# Figure S4 - S7

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(ggrepel)
library(plotly)
library(dplyr)
library(vegan)


# Fig. S5
# ABUNDANCE
# Read the ASV table CSV file into a dataframe with row names and headers
pc <- read.csv("ASV_table_37F.csv", row.names = 1, header = TRUE)
com <- as.matrix(pc)

# Read the metadata CSV file into a dataframe with row names and headers
metadata <- read.csv("nordic_metadata.csv", row.names = 1, header = TRUE)

# Read the ASV table CSV file into a dataframe with row names and headers
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)

# Add the "sample_type" as the first column of the pc dataframe
metadata_id <- metadata %>% rownames_to_column("sample_id")

# Assuming metadata_id is your dataframe
extracted_data <- metadata_id %>%
  select(sample_id, station, sample_type, time, distance, location)

extracted_data$total <- rowSums(com)

# Create a new dataset with only 'station' and 'water_depth' for sediment samples
station_depth <- metadata_id %>%
  filter(sample_type == "sediment") %>%    # Filter for sediment samples
  select(station, water_depth)              # Select only the station and water_depth columns

# Join extracted_data with station_depth based on the station column
extracted_data <- extracted_data %>%
  left_join(station_depth, by = "station")

# Update the water_depth for specific stations
extracted_data <- extracted_data %>%
  mutate(water_depth = case_when(
    station == "KH04" ~ 0,         # Set water_depth to 0 for station KH04
    station == "KH05" ~ 426,       # Set water_depth to 426 for station KH05
    TRUE ~ water_depth )) %>%  
  mutate(sample_type = case_when(
    sample_id == "KH09_03" ~ "100m_water",         
    TRUE ~ sample_type ))            


# Read the ASV table CSV file into a dataframe with row names and headers
taxon <- read.csv("ASV_assignment.csv", row.names = 1, header = TRUE)

# Filter the data based on the group and extract row names
benthic_ASV <- rownames(taxon[taxon$group %in% c("benthic", "unassigned"),])

# Subset the 'com' table to keep only rows with row names in benthic_planktonic_list
t_com <- t(com)

benthic_com <- t_com[rownames(t_com) %in% benthic_ASV, ]
benthic_com <- t(benthic_com)

extracted_data$benthic_read <- rowSums(benthic_com)

# Calculate the percentage and add it as a new column
extracted_data$percentage <- (extracted_data$benthic_read * 100) / extracted_data$total
extracted_data$richness <- specnumber(benthic_com)
extracted_data$distance <- extracted_data$distance / 1000
extracted_data$water_depth <- extracted_data$water_depth / 1000

region_color <- c("#d62728", "#B964A6", "#0196C1",  "#66a61e", "#FF7F00")

sample_color = c("#2a9d8f", "#e9c46a", "#e76f51")

extracted_data$location <- factor(extracted_data$location, 
                                  levels = c("North Svalbard", 
                                             "West Svalbard", 
                                             "North Norway", 
                                             "NE Greenland", 
                                             "Greenland Ridge"))

# Filter rows where type is "type sample"

SW <- extracted_data %>%
  filter(sample_type == "surface_water")


# Show a bubbleplot
p1 <- SW %>%
  mutate(percentage = percentage) %>%
  arrange(desc(percentage)) %>%
  mutate(sample_id = factor(sample_id, sample_id)) %>%
  ggplot(aes(x=distance, y= water_depth, size = percentage,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,19),
                        name= "Relative abundance", breaks = c(10,50,75,100)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "Time (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p1





M100 <- extracted_data %>%
  filter(sample_type == "100m_water")


# Show a bubbleplot
p2 <- M100 %>%
  mutate(percentage = percentage) %>%
  arrange(desc(percentage)) %>%
  mutate(sample_id = factor(sample_id, sample_id)) %>%
  ggplot(aes(x=distance, y= water_depth, size = percentage,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,19),
                        name= "Relative abundance", breaks = c(10,50,75,100)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p2


BW <- extracted_data %>%
  filter(sample_type == "bottom_water")


# Show a bubbleplot
p3 <- BW %>%
  mutate(percentage = percentage) %>%
  arrange(desc(percentage)) %>%
  mutate(sample_id = factor(sample_id, sample_id)) %>%
  ggplot(aes(x=distance, y= water_depth, size = percentage,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,19),
                        name= "Relative abundance", breaks = c(10,50,75,100)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p3



SED <- extracted_data %>%
  filter(sample_type == "sediment")


# Show a bubbleplot
p4 <- SED %>%
  mutate(percentage = percentage) %>%
  arrange(desc(percentage)) %>%
  mutate(sample_id = factor(sample_id, sample_id)) %>%
  ggplot(aes(x=distance, y= water_depth, size = percentage,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,19),
                        name= "Relative abundance", breaks = c(10,50,75,100)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p4


# Fig. S6
# RICHNESS
# Show a bubbleplot
p1 <- SW %>%
  mutate(richness = richness) %>%
  arrange(desc(richness)) %>%
  mutate(sample_id  = factor(sample_id , sample_id )) %>%
  ggplot(aes(x=distance, y= water_depth, size = richness,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(1, 1200), range = c(1,19),
                        name= "Relative abundance", breaks = c(100,500,750,1200)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p1


max(extracted_data$richness)


# Show a bubbleplot
p2 <- M100 %>%
  mutate(richness = richness) %>%
  arrange(desc(richness)) %>%
  mutate(sample_id  = factor(sample_id , sample_id )) %>%
  ggplot(aes(x=distance, y= water_depth, size = richness,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(1, 1200), range = c(1,19),
                        name= "Relative abundance", breaks = c(100,500,750,1200)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p2



# Show a bubbleplot
p3 <- BW %>%
  mutate(richness = richness) %>%
  arrange(desc(richness)) %>%
  mutate(sample_id  = factor(sample_id , sample_id )) %>%
  ggplot(aes(x=distance, y= water_depth, size = richness,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(1, 1200), range = c(1,19),
                        name= "Relative abundance", breaks = c(100,500,750,1200)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p3


# Show a bubbleplot
p4 <- SED %>%
  mutate(richness = richness) %>%
  arrange(desc(richness)) %>%
  mutate(sample_id  = factor(sample_id , sample_id )) %>%
  ggplot(aes(x=distance, y= water_depth, size = richness,  color = location)) +
  geom_point(aes(fill = location), alpha=0.7, pch = 21, color = "black") + 
  ylim(c(3.5,-0.2)) +
  xlim(c(-20, 400)) +
  scale_size_continuous(limits = c(1, 1200), range = c(1,19),
                        name= "Relative abundance", breaks = c(100,500,750,1200)) + 
  scale_color_manual(values = region_color) +
  scale_fill_manual(values = region_color) +
  theme_minimal()+
  labs(color = "location (month)", x = "Distance from coast (km)", y = "Water depth of station (km)") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

p4

# Fig. S7
# Load packages and colors
library(gtools)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(data.table)
library(ape)
library(psych)



merged_data <- bind_rows(SW, M100)


### Relative abundance of sinking plankton in the sediment

gam_mod1 <- mgcv::gam(SW$percentage ~ s(SW$distance, k=3))
summ_gam1 <- summary(gam_mod1)
ex_dev1 <- round(summ_gam1$dev.expl * 100, 1)
p1 <- summ_gam1$s.table[,'p-value']
if(p1 < 0.05) {
  if(p1 >= 0.01) sig1 <- "*"
  if(p1 >= 0.001) sig1 <- "**"
  if(p1 < 0.001) sig1 <- "***"
} else {sig1 <- "ns" }




gam_mod2 <- mgcv::gam(M100$percentage ~ s(M100$distance, k=3))
summ_gam2 <- summary(gam_mod2)
ex_dev2 <- round(summ_gam2$dev.expl * 100, 1)
p2 <- summ_gam1$s.table[,'p-value']
if(p2 < 0.05) {
  if(p2 >= 0.01) sig2 <- "*"
  if(p2 >= 0.001) sig2 <- "**"
  if(p2 < 0.001) sig2 <- "***"
} else {sig2 <- "ns" }

summ_gam <- paste0("Expl. deviance", 
                   "\n Surface: ", ex_dev1, "%", sig1, 
                   "\n 100m water: ", ex_dev2, "%", sig2)



# Create the plot
pp1 <- ggplot(merged_data, aes(x = percentage, y = distance, colour = sample_type)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = c("surface_water" = "#19A9A9", "100m_water" = "#053F5C")) + 
  stat_smooth(aes(fill=sample_type), method = "gam", formula = y ~ s(x, k= 2), size = 2, alpha = 0.1, se = F,  level = 0.95) +
  theme_classic() + 
  theme(legend.position="top") +
  ylab("Distance to the coast (km)") +
  xlab("Proportion of benthic reads (%)") + xlim(0,100) + ylim(0,400)+
  annotate("text", x = 80, y = 350, label = summ_gam) 

pp1 


### Richness of sinking plankton in the sediment

gam_mod3 <- mgcv::gam(SW$richness ~ s(SW$distance, k=3))
summ_gam3 <- summary(gam_mod3)
ex_dev3 <- round(summ_gam3$dev.expl * 100, 1)
p3 <- summ_gam3$s.table[,'p-value']
if(p3 < 0.05) {
  if(p3 >= 0.01) sig3 <- "*"
  if(p3 >= 0.001) sig3 <- "**"
  if(p3 < 0.001) sig3 <- "***"
} else {sig3 <- "ns" }




gam_mod4 <- mgcv::gam(M100$richness ~ s(M100$distance, k=3))
summ_gam4 <- summary(gam_mod4)
ex_dev4 <- round(summ_gam4$dev.expl * 100, 1)
p4 <- summ_gam4$s.table[,'p-value']
if(p4 < 0.05) {
  if(p4 >= 0.01) sig4 <- "*"
  if(p4 >= 0.001) sig4 <- "**"
  if(p4 < 0.001) sig4 <- "***"
} else {sig4 <- "ns" }

summ_gam <- paste0("Expl. deviance", 
                   "\n Surface: ", ex_dev3, "%", sig3, 
                   "\n 100m water: ", ex_dev4, "%", sig4)


# Create the plot
pp2 <- ggplot(merged_data, aes(x = richness, y = distance, colour = sample_type)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = c("surface_water" = "#19A9A9", "100m_water" = "#053F5C")) + 
  stat_smooth(aes(fill=sample_type), method = "gam", formula = y ~ s(x, k= 2), size = 2, alpha = 0.1, se = F,  level = 0.95) +
  theme_classic() + 
  theme(legend.position="top") +
  ylab("Distance to the coast (km)") +
  xlab("Number of benthic ASV") + xlim(0,280) + ylim(0,400)+
  annotate("text", x = 250, y = 350, label = summ_gam) 

pp2 

grid.arrange(pp1, pp2, ncol = 1)




