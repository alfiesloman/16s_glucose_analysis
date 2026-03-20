# Microbiome Analysis Workflow in R - GENUS LEVEL
# This script performs a complete analysis workflow on genus-level abundance data.

# --- 1. Load Libraries ---
library(phyloseq)
library(tidyverse) # Includes ggplot2, dplyr, etc.
library(vegan)

# --- Define File Paths ---
abundance_file <- "genus_abundance_table.csv"
metadata_file <- "metadata.csv"
output_prefix <- "genus_analysis_R"

# --- Load Data ---
# row.names=1 sets the first column (taxa names) as the row names
abundance_table_df <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
# row.names=1 sets the first column (SampleID) as the row names
metadata_table <- read.csv(metadata_file, row.names = 1)


# --- FIX: Robustly convert the abundance data to a numeric matrix ---
# This multi-step process is more robust than direct conversion. It handles
# cases where R might silently keep data as non-numeric characters or factors.

# 1. Unlist the data frame into a single vector and force it to be numeric.
#    This will introduce NA for any values that can't be converted.
numeric_vector <- as.numeric(unlist(abundance_table_df))

# 2. Reshape this vector back into a matrix with the original dimensions.
abundance_matrix <- matrix(numeric_vector,
                           nrow = nrow(abundance_table_df),
                           ncol = ncol(abundance_table_df))

# 3. Restore the original row (taxa) and column (sample) names.
rownames(abundance_matrix) <- rownames(abundance_table_df)
colnames(abundance_matrix) <- colnames(abundance_table_df)

# 4. Replace any NAs that were introduced during coercion with 0.
abundance_matrix[is.na(abundance_matrix)] <- 0

# Ensure metadata is a data.frame for phyloseq
metadata_table <- as.data.frame(metadata_table)

# FIX: Convert 'Day' column to a factor for plotting aesthetics
metadata_table$Day <- as.factor(metadata_table$Day)


# --- Create Phyloseq Object ---
# This object conveniently bundles your abundance and metadata tables.
OTU <- otu_table(abundance_matrix, taxa_are_rows = TRUE)
META <- sample_data(metadata_table)

# FIX: Create and add a taxonomy table to resolve heatmap error
taxa_names <- rownames(abundance_matrix)
tax_table_df <- data.frame(Genus = taxa_names) # Changed to Genus
rownames(tax_table_df) <- taxa_names
TAX <- tax_table(as.matrix(tax_table_df))

physeq <- phyloseq(OTU, META, TAX)

# --- FIX: Create and apply descriptive sample names AFTER phyloseq object creation ---
# This is more robust and avoids errors from metadata/abundance table mismatches.
meta_df <- as(sample_data(physeq), "data.frame")
descriptive_sample_names <- paste(meta_df$Glucose_Concentration, "_D", meta_df$Day, sep = "")
sample_names(physeq) <- descriptive_sample_names


# Normalize to relative abundance (for beta diversity and heatmaps)
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

print("Data loaded into phyloseq object and normalized successfully.")

# --- 2. Alpha Diversity Analysis ---
print("--- Step 2: Alpha Diversity Analysis ---")

# --- FIX: Changed visualization from boxplot to line plot ---
# A line plot is more appropriate for visualizing time-series data with single replicates.

# --- FIX: Create a clean data frame for plotting to resolve aesthetic error ---
# Combine diversity estimates and metadata into a single, clean data frame.
alpha_div_estimates <- estimate_richness(physeq, measures = "Shannon")
metadata_df_for_plot <- as(sample_data(physeq), "data.frame")
plot_data <- cbind(alpha_div_estimates, metadata_df_for_plot)


# Create the line plot using the new clean data frame
alpha_plot <- ggplot(plot_data, aes(x = Day, y = Shannon, group = Glucose_Concentration, color = Glucose_Concentration)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Alpha Diversity (Shannon Index) Over Time - Genus Level",
       x = "Day",
       y = "Shannon Diversity Index",
       color = "Glucose") +
  theme_bw()

# Save the plot
ggsave(paste0(output_prefix, "_alpha_diversity_lineplot.png"), plot = alpha_plot, width = 8, height = 6)
print(paste("Alpha diversity line plot saved to", paste0(output_prefix, "_alpha_diversity_lineplot.png")))


# Statistical test (Kruskal-Wallis) for Glucose effect on Day 4
print("Running Kruskal-Wallis test for Glucose effect on Day 4...")
day4_data_for_test <- subset(plot_data, Day == 4)
kw_test <- kruskal.test(Shannon ~ Glucose_Concentration, data = day4_data_for_test)
print(kw_test)

# --- Pairwise Alpha Diversity Tests by Day (Wilcoxon Rank-Sum Test) ---
print("--- Running Pairwise Alpha Diversity Tests between Days ---")

# Compare Day 1 vs Day 2
print("--- Day 1 vs Day 2 ---")
data_d1_d2 <- subset(plot_data, Day %in% c("1", "2"))
wilcox_d1_d2 <- wilcox.test(Shannon ~ Day, data = data_d1_d2)
print(wilcox_d1_d2)

# Compare Day 1 vs Day 4
print("--- Day 1 vs Day 4 ---")
data_d1_d4 <- subset(plot_data, Day %in% c("1", "4"))
wilcox_d1_d4 <- wilcox.test(Shannon ~ Day, data = data_d1_d4)
print(wilcox_d1_d4)

# Compare Day 2 vs Day 4
print("--- Day 2 vs Day 4 ---")
data_d2_d4 <- subset(plot_data, Day %in% c("2", "4"))
wilcox_d2_d4 <- wilcox.test(Shannon ~ Day, data = data_d2_d4)
print(wilcox_d2_d4)


# --- 3. Beta Diversity Analysis ---
print("--- Step 3: Beta Diversity Analysis ---")

# Calculate Bray-Curtis dissimilarity and perform PCoA using the relative abundance data
pcoa_bray <- ordinate(physeq_rel, method = "PCoA", distance = "bray")

# Plot the PCoA
beta_plot <- plot_ordination(physeq_rel, pcoa_bray, color = "Glucose_Concentration", shape = "Day") +
  geom_point(size = 5) +
  labs(title = "Beta Diversity (Bray-Curtis PCoA) - Genus Level",
       color = "Glucose",
       shape = "Day") +
  theme_bw()

# Save the plot
ggsave(paste0(output_prefix, "_beta_diversity_pcoa.png"), plot = beta_plot, width = 8, height = 7)
print(paste("Beta diversity PCoA plot saved to", paste0(output_prefix, "_beta_diversity_pcoa.png")))

# To test for significance of clustering, we use PERMANOVA from the vegan package
print("Running PERMANOVA test on combined model...")
bray_dist <- phyloseq::distance(physeq_rel, method = "bray")
permanova_result <- adonis2(bray_dist ~ Glucose_Concentration + Day, data = as(sample_data(physeq), "data.frame"))
print(permanova_result)

# --- PERMANOVA by individual factors to determine a hierarchy of effect sizes ---
print("--- Running PERMANOVA for each factor to compare R2 values ---")

# Test the effect of Day alone
print("--- Testing the effect of Day alone ---")
permanova_day <- adonis2(bray_dist ~ Day, data = as(sample_data(physeq), "data.frame"))
print(permanova_day)

# Test the effect of Glucose Concentration alone
print("--- Testing the effect of Glucose Concentration alone ---")
permanova_glucose <- adonis2(bray_dist ~ Glucose_Concentration, data = as(sample_data(physeq), "data.frame"))
print(permanova_glucose)

# --- NEW: Pairwise PERMANOVA to test specific hypotheses ---
print("--- Running Pairwise PERMANOVA to test Water vs 20 Glucose ---")
# Create a subset of the data containing only the groups we want to compare
physeq_subset_wv20 <- subset_samples(physeq_rel, Glucose_Concentration %in% c("Water", "20"))
# Calculate Bray-Curtis distance for this subset
bray_dist_subset <- phyloseq::distance(physeq_subset_wv20, method = "bray")
# Run PERMANOVA on the subset
permanova_pairwise <- adonis2(bray_dist_subset ~ Glucose_Concentration, data = as(sample_data(physeq_subset_wv20), "data.frame"))
print(permanova_pairwise)


# --- 4. Community Composition Barplots (Genus) ---
print("--- Step 4: Community Composition Barplots (Genus) ---")

# Aggregate the relative abundance data (already normalized in physeq_rel)
physeq_genus_glom <- tax_glom(physeq_rel, taxrank="Genus")

# Melt the phyloseq object for ggplot2
plot_data_melt <- psmelt(physeq_genus_glom)

# Identify the top 10 most abundant genera for the plot, grouping the rest as 'Other'
top_n <- 10
top_taxa <- names(sort(tapply(plot_data_melt$Abundance, plot_data_melt$Genus, sum), decreasing=TRUE)[1:top_n])
plot_data_melt$Taxa_Group <- factor(ifelse(plot_data_melt$Genus %in% top_taxa,
                                           as.character(plot_data_melt$Genus),
                                           "Other"),
                                    levels = c(top_taxa, "Other")) # REVERTED: Removed rev() for standard stacking

# Filter out the "Other" group from the plot data
plot_data_filtered <- plot_data_melt %>%
  filter(Taxa_Group != "Other")

# Order the samples by Day and Glucose concentration for clean plotting
meta_df_for_ordering_bp <- as(sample_data(physeq_rel), "data.frame")
meta_df_for_ordering_bp$Day <- factor(meta_df_for_ordering_bp$Day, levels = c("1", "2", "4"), ordered = TRUE)
meta_df_for_ordering_bp$Glucose_Concentration <- factor(meta_df_for_ordering_bp$Glucose_Concentration, levels = c("Water", "5", "10", "20"), ordered = TRUE)
meta_df_for_ordering_bp <- meta_df_for_ordering_bp[order(meta_df_for_ordering_bp$Day, meta_df_for_ordering_bp$Glucose_Concentration), ]
ordered_samples_names <- rownames(meta_df_for_ordering_bp)

# Apply the new sample order to the melted data frame
plot_data_filtered$Sample <- factor(plot_data_filtered$Sample, levels = ordered_samples_names)

# Define a custom color palette with a dark to light green gradient for the top 10 genera
# We will use the RColorBrewer 'Greens' palette which has 9 colors, so we need 10.
# We'll use 10 distinct hex codes from a green gradient.
green_gradient <- c(
  "#00441B", # Darkest
  "#006D2C",
  "#238B45",
  "#41AB5D",
  "#74C476",
  "#A1D99B",
  "#C7E9C0",
  "#E5F5E0",
  "#F7FCF5", # Lightest from Greens palette
  "#33a02c"  # Added a distinct tenth green for visual separation
)

# Identify all unique genus names present in the filtered data
taxa_to_color <- levels(droplevels(plot_data_filtered$Taxa_Group))

# Assign the green gradient colors to the genera, repeating or taking the first N
# The first 10 colors from the gradient will be used for the 10 taxa
color_map <- setNames(green_gradient[1:length(taxa_to_color)], taxa_to_color)


barplot_genus <- ggplot(plot_data_filtered, aes(x = Sample, y = Abundance * 100, fill = Taxa_Group)) + # Multiply Abundance by 100
  geom_bar(stat = "identity", position = "stack") +
  # Use the descriptive sample names for the x-axis labels
  scale_x_discrete(labels = descriptive_sample_names) +
  scale_fill_manual(values = color_map) + # CHANGED to custom green gradient
  # Title changed to match reference plot aesthetic
  labs(title = paste("Community Composition of Top", top_n, "Genus Level"),
       x = "Sample (Glucose Concentration and Day)",
       y = "Relative Abundance (%)", # Y-axis label updated
       fill = "Genus") +
  # Faceting removed to match reference plot aesthetic
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        strip.background = element_rect(fill = "grey80"))

# Save the plot
ggsave(paste0(output_prefix, "_composition_barplot_genus.png"), plot = barplot_genus, width = 12, height = 7)
print(paste("Genus composition barplot saved to", paste0(output_prefix, "_composition_barplot_genus.png")))

# --- 5. Differential Abundance Analysis (Heatmap) ---
print("--- Step 5: Differential Abundance Analysis (Heatmap) ---")

# Get the top 10 most abundant genera to make the plot cleaner
top_taxa_n <- 10
top_taxa <- names(sort(taxa_sums(physeq), decreasing = TRUE)[1:top_taxa_n])
physeq_top_taxa <- prune_taxa(top_taxa, physeq_rel)

# --- FIX: Create a custom sample order for the heatmap x-axis ---
# Get metadata to determine the order
meta_df_for_ordering <- as(sample_data(physeq_top_taxa), "data.frame")
# Create ordered factors to sort by
meta_df_for_ordering$Day <- factor(meta_df_for_ordering$Day, levels = c("1", "2", "4"), ordered = TRUE)
meta_df_for_ordering$Glucose_Concentration <- factor(meta_df_for_ordering$Glucose_Concentration, levels = c("Water", "5", "10", "20"), ordered = TRUE)
# Order the data frame based on Day, then Glucose
meta_df_for_ordering <- meta_df_for_ordering[order(meta_df_for_ordering$Day, meta_df_for_ordering$Glucose_Concentration), ]
# Get the correctly ordered sample names from the row names of the ordered metadata
ordered_samples <- rownames(meta_df_for_ordering)


# FIX: Create the heatmap, providing the new custom order
heatmap_plot <- plot_heatmap(physeq_top_taxa,
                             method = NULL, # Disable automatic sample reordering
                             sample.order = ordered_samples, # Provide our custom order
                             taxa.label = "Genus", # Changed to Genus
                             low = "#FFFFFF", high = "#008000", # Color gradient
                             title = paste("Top", top_taxa_n, "Most Abundant Genera")) + # Changed Title
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) # Rotate labels

# Save the plot with dimensions suitable for the single heatmap format
ggsave(paste0(output_prefix, "_abundance_heatmap.png"), plot = heatmap_plot, width = 10, height = 8)
print(paste("Abundance heatmap saved to", paste0(output_prefix, "_abundance_heatmap.png")))

print("Genus Workflow complete.")







#############################################################################
#         SPECIES ###########################################################
#############################################################################
# Microbiome Analysis Workflow in R
# This script performs the same analysis as the Python equivalent, using
# the phyloseq, vegan, and ggplot2 packages.

# --- 1. Load Libraries ---
library(phyloseq)
library(tidyverse) # Includes ggplot2, dplyr, etc.
library(vegan)

print("--- Step 1: Loading and Preparing Data ---")

# --- Define File Paths ---
abundance_file <- "species_abundance_table.csv"
metadata_file <- "metadata.csv"
output_prefix <- "species_analysis_R"

# --- Load Data ---
# row.names=1 sets the first column (taxa names) as the row names
abundance_table_df <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
# row.names=1 sets the first column (SampleID) as the row names
metadata_table <- read.csv(metadata_file, row.names = 1)


# --- FIX: Robustly convert the abundance data to a numeric matrix ---
# This multi-step process is more robust than direct conversion. It handles
# cases where R might silently keep data as non-numeric characters or factors.

# 1. Unlist the data frame into a single vector and force it to be numeric.
#    This will introduce NA for any values that can't be converted.
numeric_vector <- as.numeric(unlist(abundance_table_df))

# 2. Reshape this vector back into a matrix with the original dimensions.
abundance_matrix <- matrix(numeric_vector,
                           nrow = nrow(abundance_table_df),
                           ncol = ncol(abundance_table_df))

# 3. Restore the original row (taxa) and column (sample) names.
rownames(abundance_matrix) <- rownames(abundance_table_df)
colnames(abundance_matrix) <- colnames(abundance_table_df)

# 4. Replace any NAs that were introduced during coercion with 0.
abundance_matrix[is.na(abundance_matrix)] <- 0

# Ensure metadata is a data.frame for phyloseq
metadata_table <- as.data.frame(metadata_table)

# FIX: Convert 'Day' column to a factor for plotting aesthetics
metadata_table$Day <- as.factor(metadata_table$Day)


# --- Create Phyloseq Object ---
# This object conveniently bundles your abundance and metadata tables.
OTU <- otu_table(abundance_matrix, taxa_are_rows = TRUE)
META <- sample_data(metadata_table)

# FIX: Create and add a taxonomy table to resolve heatmap error
taxa_names <- rownames(abundance_matrix)
tax_table_df <- data.frame(Species = taxa_names) # Changed from Genus to Species
rownames(tax_table_df) <- taxa_names
TAX <- tax_table(as.matrix(tax_table_df))

physeq <- phyloseq(OTU, META, TAX)

# --- FIX: Create and apply descriptive sample names AFTER phyloseq object creation ---
# This is more robust and avoids errors from metadata/abundance table mismatches.
meta_df <- as(sample_data(physeq), "data.frame")
descriptive_sample_names <- paste(meta_df$Glucose_Concentration, "_D", meta_df$Day, sep = "")
sample_names(physeq) <- descriptive_sample_names


# Normalize to relative abundance (for beta diversity and heatmaps)
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

print("Data loaded into phyloseq object and normalized successfully.")

# --- 2. Alpha Diversity Analysis ---
print("--- Step 2: Alpha Diversity Analysis ---")

# --- FIX: Changed visualization from boxplot to line plot ---
# A line plot is more appropriate for visualizing time-series data with single replicates.

# Get all alpha diversity estimates at once
alpha_div_estimates <- estimate_richness(physeq, measures = "Shannon")
# Get metadata as a standard data frame
metadata_df_for_plot <- as(sample_data(physeq), "data.frame")

# --- FIX: Create a clean data frame for plotting to resolve aesthetic error ---
# Combine diversity estimates and metadata into a single, clean data frame.
plot_data <- cbind(alpha_div_estimates, metadata_df_for_plot)


# Create the line plot using the new clean data frame
alpha_plot <- ggplot(plot_data, aes(x = Day, y = Shannon, group = Glucose_Concentration, color = Glucose_Concentration)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Alpha Diversity (Shannon Index) Over Time - Species Level",
       x = "Day",
       y = "Shannon Diversity Index",
       color = "Glucose") +
  theme_bw()

# Save the plot
ggsave(paste0(output_prefix, "_alpha_diversity_lineplot.png"), plot = alpha_plot, width = 8, height = 6)
print(paste("Alpha diversity line plot saved to", paste0(output_prefix, "_alpha_diversity_lineplot.png")))


# Statistical test (Kruskal-Wallis) for Glucose effect on Day 4
print("Running Kruskal-Wallis test for Glucose effect on Day 4...")
day4_data_for_test <- subset(plot_data, Day == 4)
kw_test <- kruskal.test(Shannon ~ Glucose_Concentration, data = day4_data_for_test)
print(kw_test)

# --- Pairwise Alpha Diversity Tests by Day (Wilcoxon Rank-Sum Test) ---
print("--- Running Pairwise Alpha Diversity Tests between Days ---")

# Compare Day 1 vs Day 2
print("--- Day 1 vs Day 2 ---")
data_d1_d2 <- subset(plot_data, Day %in% c("1", "2"))
wilcox_d1_d2 <- wilcox.test(Shannon ~ Day, data = data_d1_d2)
print(wilcox_d1_d2)

# Compare Day 1 vs Day 4
print("--- Day 1 vs Day 4 ---")
data_d1_d4 <- subset(plot_data, Day %in% c("1", "4"))
wilcox_d1_d4 <- wilcox.test(Shannon ~ Day, data = data_d1_d4)
print(wilcox_d1_d4)

# Compare Day 2 vs Day 4
print("--- Day 2 vs Day 4 ---")
data_d2_d4 <- subset(plot_data, Day %in% c("2", "4"))
wilcox_d2_d4 <- wilcox.test(Shannon ~ Day, data = data_d2_d4)
print(wilcox_d2_d4)


# --- 3. Beta Diversity Analysis ---
print("--- Step 3: Beta Diversity Analysis ---")

# Calculate Bray-Curtis dissimilarity and perform PCoA using the relative abundance data
pcoa_bray <- ordinate(physeq_rel, method = "PCoA", distance = "bray")

# Plot the PCoA
beta_plot <- plot_ordination(physeq_rel, pcoa_bray, color = "Glucose_Concentration", shape = "Day") +
  geom_point(size = 5) +
  labs(title = "Beta Diversity (Bray-Curtis PCoA) - Species Level",
       color = "Glucose",
       shape = "Day") +
  theme_bw()

# Save the plot
ggsave(paste0(output_prefix, "_beta_diversity_pcoa.png"), plot = beta_plot, width = 8, height = 7)
print(paste("Beta diversity PCoA plot saved to", paste0(output_prefix, "_beta_diversity_pcoa.png")))

# To test for significance of clustering, we use PERMANOVA from the vegan package
print("Running PERMANOVA test on combined model...")
bray_dist <- phyloseq::distance(physeq_rel, method = "bray")
# FIX: Simplified the model to test main effects, avoiding the "no residual" error
permanova_result <- adonis2(bray_dist ~ Glucose_Concentration + Day, data = as(sample_data(physeq), "data.frame"))
print(permanova_result)

# --- PERMANOVA by individual factors to determine a hierarchy of effect sizes ---
print("--- Running PERMANOVA for each factor to compare R2 values ---")

# Test the effect of Day alone
print("--- Testing the effect of Day alone ---")
permanova_day <- adonis2(bray_dist ~ Day, data = as(sample_data(physeq), "data.frame"))
print(permanova_day)

# Test the effect of Glucose Concentration alone
print("--- Testing the effect of Glucose Concentration alone ---")
permanova_glucose <- adonis2(bray_dist ~ Glucose_Concentration, data = as(sample_data(physeq), "data.frame"))
print(permanova_glucose)

# --- NEW: Pairwise PERMANOVA to test specific hypotheses ---
print("--- Running Pairwise PERMANOVA to test Water vs 20 Glucose ---")
# Create a subset of the data containing only the groups we want to compare
physeq_subset_wv20 <- subset_samples(physeq_rel, Glucose_Concentration %in% c("Water", "20"))
# Calculate Bray-Curtis distance for this subset
bray_dist_subset <- phyloseq::distance(physeq_subset_wv20, method = "bray")
# Run PERMANOVA on the subset
permanova_pairwise <- adonis2(bray_dist_subset ~ Glucose_Concentration, data = as(sample_data(physeq_subset_wv20), "data.frame"))
print(permanova_pairwise)


# --- 4. Community Composition Barplots (Species) ---
print("--- Step 4: Community Composition Barplots (Species) ---")

# Aggregate the relative abundance data (already normalized in physeq_rel)
physeq_species_glom <- tax_glom(physeq_rel, taxrank="Species")

# Melt the phyloseq object for ggplot2
plot_data_melt_species <- psmelt(physeq_species_glom)

# Identify the top 10 species for the plot, grouping the rest as 'Other'
top_n <- 10
top_taxa_species <- names(sort(tapply(plot_data_melt_species$Abundance, plot_data_melt_species$Species, sum), decreasing=TRUE)[1:top_n])
plot_data_melt_species$Taxa_Group <- factor(ifelse(plot_data_melt_species$Species %in% top_taxa_species,
                                                   as.character(plot_data_melt_species$Species),
                                                   "Other"),
                                            # We will re-order the levels below to ensure grouping and stacking consistency
                                            levels = c(top_taxa_species, "Other"))

# Filter out the "Other" group from the plot data (as done for Genus)
plot_data_filtered_species <- plot_data_melt_species %>%
  filter(Taxa_Group != "Other")

# --- FIX: Custom Ordering for Species Stacking ---
# 1. Extract the Genus from the Species names (assuming "Genus species" format)
#    We use a placeholder Genus for single-name species if needed, but here we assume the Genus is the first word.
#    Since we don't have a formal taxonomy table for species-level data, we extract the first word as the Genus.
genus_from_species <- sapply(strsplit(taxa_to_color_species, " "), `[`, 1)

# 2. Create a temporary data frame to calculate the overall abundance for ordering
temp_abundance_df <- plot_data_filtered_species %>%
  group_by(Taxa_Group) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Genus = sapply(strsplit(as.character(Taxa_Group), " "), `[`, 1))

# 3. Calculate the total abundance for each *Genus* group, then the mean species abundance
#    This allows us to order by major Genus, and then by the abundance of the species within that Genus.
ordered_taxa_species_df <- temp_abundance_df %>%
  group_by(Genus) %>%
  mutate(GenusTotalAbundance = sum(TotalAbundance)) %>%
  ungroup() %>%
  # Order first by Genus Total Abundance (descending), then by Species Total Abundance (descending)
  arrange(desc(GenusTotalAbundance), desc(TotalAbundance))

# 4. Get the final, custom-ordered list of species
custom_ordered_species <- ordered_taxa_species_df$Taxa_Group
custom_ordered_species_rev <- rev(custom_ordered_species) # Reverse for bottom-up stacking

# 5. Apply the new custom order (reversed) to the Taxa_Group factor levels
plot_data_filtered_species$Taxa_Group <- factor(plot_data_filtered_species$Taxa_Group, levels = custom_ordered_species_rev)


# Order the samples by Day and Glucose concentration for clean plotting
meta_df_for_ordering_bp_s <- as(sample_data(physeq_rel), "data.frame")
meta_df_for_ordering_bp_s$Day <- factor(meta_df_for_ordering_bp_s$Day, levels = c("1", "2", "4"), ordered = TRUE)
meta_df_for_ordering_bp_s$Glucose_Concentration <- factor(meta_df_for_ordering_bp_s$Glucose_Concentration, levels = c("Water", "5", "10", "20"), ordered = TRUE)
meta_df_for_ordering_bp_s <- meta_df_for_ordering_bp_s[order(meta_df_for_ordering_bp_s$Day, meta_df_for_ordering_bp_s$Glucose_Concentration), ]
ordered_samples_names_species <- rownames(meta_df_for_ordering_bp_s)

# Apply the new sample order to the melted data frame
plot_data_filtered_species$Sample <- factor(plot_data_filtered_species$Sample, levels = ordered_samples_names_species)

# Define the custom color palette with a dark to light green gradient (reused from Genus)
# We need to make sure we have enough colors for the ordered list
green_gradient <- c(
  "#00441B", # Darkest (bottom of the stack)
  "#006D2C",
  "#238B45",
  "#41AB5D",
  "#74C476",
  "#A1D99B",
  "#C7E9C0",
  "#E5F5E0",
  "#F7FCF5",
  "#33a02c"  # Tenth green (top of the stack)
)
# Ensure the color map is applied in the order of the new factor levels
# The factor levels are REV so the first color of the gradient is at the bottom
taxa_to_color_species <- levels(plot_data_filtered_species$Taxa_Group)
color_map_species <- setNames(green_gradient[1:length(taxa_to_color_species)], taxa_to_color_species)


barplot_species <- ggplot(plot_data_filtered_species, aes(x = Sample, y = Abundance * 100, fill = Taxa_Group)) + # Multiply Abundance by 100
  geom_bar(stat = "identity", position = "stack") +
  # Use the descriptive sample names for the x-axis labels
  scale_x_discrete(labels = descriptive_sample_names) +
  scale_fill_manual(values = color_map_species) + # Using custom green gradient and order
  # Title changed to match reference plot aesthetic, updated for Species
  labs(title = paste("Community Composition of Top", top_n, "Species Level (Ordered by Genus)"),
       x = "Sample (Glucose Concentration and Day)",
       y = "Relative Abundance (%)", # Y-axis label updated
       fill = "Species") +
  # Faceting removed to match reference plot aesthetic
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        strip.background = element_rect(fill = "grey80"))

# Save the plot
ggsave(paste0(output_prefix, "_composition_barplot_species.png"), plot = barplot_species, width = 12, height = 7)
print(paste("Species composition barplot saved to", paste0(output_prefix, "_composition_barplot_species.png")))


# --- 5. Differential Abundance Analysis (Heatmap) ---
print("--- Step 5: Differential Abundance Analysis (Heatmap) ---")

# This is a simplified approach, focusing on the most abundant taxa for visualization.
# For rigorous stats, DESeq2 or ANCOM-BC are recommended.

# Get the top 10 most abundant species to make the plot cleaner
top_taxa_n <- 10
top_taxa <- names(sort(taxa_sums(physeq), decreasing = TRUE)[1:top_taxa_n])
physeq_top_taxa <- prune_taxa(top_taxa, physeq_rel)

# --- FIX: Create a custom sample order for the heatmap x-axis ---
# Get metadata to determine the order
meta_df_for_ordering <- as(sample_data(physeq_top_taxa), "data.frame")
# Create ordered factors to sort by
meta_df_for_ordering$Day <- factor(meta_df_for_ordering$Day, levels = c("1", "2", "4"), ordered = TRUE)
meta_df_for_ordering$Glucose_Concentration <- factor(meta_df_for_ordering$Glucose_Concentration, levels = c("Water", "5", "10", "20"), ordered = TRUE)
# Order the data frame based on Day, then Glucose
meta_df_for_ordering <- meta_df_for_ordering[order(meta_df_for_ordering$Day, meta_df_for_ordering$Glucose_Concentration), ]
# Get the correctly ordered sample names from the row names of the ordered metadata
ordered_samples <- rownames(meta_df_for_ordering)


# FIX: Create the heatmap, providing the new custom order
heatmap_plot <- plot_heatmap(physeq_top_taxa,
                             method = NULL, # Disable automatic sample reordering
                             sample.order = ordered_samples, # Provide our custom order
                             taxa.label = "Species", # Changed from Genus to Species
                             low = "#FFFFFF", high = "#008000", # Color gradient
                             title = paste("Top", top_taxa_n, "Most Abundant Species")) + # Changed Title
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) # Rotate labels

# Save the plot with dimensions suitable for the single heatmap format
ggsave(paste0(output_prefix, "_abundance_heatmap.png"), plot = heatmap_plot, width = 10, height = 8)
print(paste("Abundance heatmap saved to", paste0(output_prefix, "_abundance_heatmap.png")))

print("Species Workflow complete.")
