################################################################################
# R Script for Trebouxia Photobiont Community Analysis and Figure Generation
################################################################################
# Project: SPRUCE Evernia mesomorpha photobiont analysis
# Purpose: Analyze photobiont OTU data and generate Figure 8 for:
#          Meyer, Valentin et al. "Climate warming causes photobiont 
#          degradation and C starvation in a boreal climate sentinel lichen"
# Author: Abigail Meyer
# Updated: July 6 2022
#
# This script processes Trebouxia OTU data from QIIME2 output to:
# 1. Rename and aggregate OTUs based on phylogenetic placement
# 2. Visualize photobiont community composition across treatments
# 3. Create phylogenetic tree visualizations
# 4. Analyze alpha diversity differences between core and tip samples
# 5. Test for warming effects on photobiont community composition
################################################################################

################################################################################
# SETUP: Set Working Directory
################################################################################
# Update this path to match your local directory structure
setwd(~/Desktop/PhD/Projects/Evernia_SPRUCE/Evernia_R_code/Manuscript_Files)

################################################################################
# Load Required Libraries
################################################################################
library(dplyr)      # Data manipulation and aggregation
library(tidyr)      # Data reshaping (wide to long format)
library(ggplot2)    # Data visualization and figure generation

################################################################################
# PART 1: DATA IMPORT
################################################################################
# Read in the three main data files from QIIME2 processing

# Metadata: Sample information including treatment conditions
# Columns: sample_id, individual, thallus_type, temp, CO2
treb_meta <- read.csv("Trebouxia_Metadata.csv")

# OTU table: Read counts for each OTU in each sample
# Rows = samples, Columns = OTUs (named by reference match from Muggia et al. 2020)
treb_otu <- read.csv("Trebouxia_OTU_table.csv")

# Convert factors to character for easier manipulation
treb_meta$sample_id <- as.character(treb_meta$sample_id)
treb_meta$temp <- as.character(treb_meta$temp)

################################################################################
# PART 2: RENAME UNKNOWN OTUs
################################################################################
# Unknown OTUs (UNK_01 through UNK_07) were sequences that didn't match the
# reference database at 97% similarity. These were placed in the Muggia et al.
# 2020 multilocus phylogeny using RAxML to determine their clade identity.
#
# Naming convention: 
# - I = Trebouxia clade 'I' (impressa/gelatinosa group)
# - A = Trebouxia clade 'A' (arboricola/gigantea group)  
# - S = Trebouxia clade 'S' (simplex/jamesii group)
# Number after indicates placement within that clade

treb_otu <- treb_otu %>%
  rename(I02 = UNK_01,      # Placed in clade I, group 02
         A13_1 = UNK_03,    # Placed in clade A, group 13 (first duplicate)
         S01 = UNK_04,      # Placed in clade S, group 01
         I18 = UNK_05,      # Placed in clade I, group 18
         A13_2 = UNK_06,    # Placed in clade A, group 13 (second duplicate)
         A04 = UNK_07       # Placed in clade A, group 04
  )

# Remove UNK_02 (19th column) - failed to align properly in phylogenetic analysis
treb_otu <- treb_otu[, -19]

################################################################################
# PART 3: AGGREGATE DUPLICATE OTUs
################################################################################
# Some OTUs have multiple representative sequences but belong to the same
# phylogenetic group (e.g., A13_1 and A13_2). Sum read counts for OTUs that
# share the same first 3 characters of their name.
#
# split.default: Splits columns by first 3 characters of column name
# rowSums: Sums read counts across duplicate columns

treb_otu <- lapply(split.default(treb_otu, substr(names(treb_otu), 1, 3)), 
                   rowSums, na.rm = TRUE)
treb_otu <- as.data.frame(treb_otu)

# Restore sample_id column and move to first position
treb_otu <- treb_otu %>%
  rename(sample_id = sam)
treb_otu <- treb_otu %>%
  select(sample_id, everything())

################################################################################
# PART 4: MERGE METADATA WITH OTU TABLE
################################################################################
# Combine treatment information with OTU counts for analysis
treb <- merge(treb_meta, treb_otu, by = "sample_id")

################################################################################
# PART 5: RESHAPE DATA FOR VISUALIZATION
################################################################################
# Convert from wide format (one column per OTU) to long format 
# (one row per sample-OTU combination). This format is required for ggplot2.
#
# gather: Converts columns 6-17 (the OTU columns) into two columns:
#   - OTU: The OTU name
#   - sequences: The read count

treb_long <- treb %>%
  gather(OTU, sequences, 6:17)

################################################################################
# PART 6: CREATE FIGURE 8A - Stacked Bar Plot of OTU Composition
################################################################################
# Visualize the proportion of each OTU in core vs. tip samples across 
# temperature treatments
#
# Key features:
# - Each bar represents one thallus type (core or tips) at one temperature
# - Colors represent different OTUs (organized by phylogenetic clade)
# - Y-axis shows proportion (0-100%)

library(scales)  # For percent formatting on y-axis

ggplot(treb_long, 
       aes(fill = factor(OTU, levels = c("A04", "A33", "A13", 
                                         "I08", "I11", "I01", "I02", "I18", 
                                         "S01", "S02", "S06", "S10")), 
           y = sequences, 
           x = thallus_type)) +
  geom_bar(stat = "identity", position = "fill") +
  # Facet by temperature treatment (in order of increasing warming)
  facet_grid(~ factor(temp, levels = c("0", "2", "4", "6", "8", "10"))) +
  # Clean theme without grid lines
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  # Manually set colors for each OTU
  # Color scheme: Similar clades have similar hues
  # A clade = browns/yellows, I clade = blues, S clade = greens
  scale_fill_manual(values = c("A04" = "#6e5c06",   # Dark brown
                                "A13" = "#a5994a",   # Medium brown
                                "A33" = "#ddda8c",   # Light yellow-brown
                                "I01" = "#5864AB",   # Medium blue
                                "I02" = "#323281",   # Dark blue
                                "I08" = "#B5CCFB",   # Light blue
                                "I11" = "#8397D3",   # Medium-light blue
                                "I18" = "#0C0056",   # Very dark blue
                                "S01" = "#9CF7C7",   # Light green
                                "S02" = "#6abc90",   # Medium-light green
                                "S06" = "#39835e",   # Medium-dark green
                                "S10" = "#014f2f")) + # Dark green
  guides(fill = guide_legend(title = "OTU")) +
  xlab("Thallus Type by Temperature Differential °C") +
  ylab("Percent") +
  scale_y_continuous(labels = percent)

################################################################################
# PART 7: CREATE PHYLOGENETIC TREE FOR FIGURE 8A LEGEND
################################################################################
# Load phylogenetic tree packages
library(ggtree)   # Grammar of graphics for phylogenetic trees
library(treeio)   # Tree input/output
library(ape)      # Phylogenetic analysis tools

# Read in RAxML best tree with representative sequences from Muggia et al. 2020
# This tree includes our new OTUs placed within the published Trebouxia phylogeny
treb_tree <- read.tree(file = "RAxML_bestTree.legend.tree")

# View tip labels to identify what to remove
treb_tree$tip.label

# Remove two tips that aren't relevant for the legend (14th and 15th tips)
drop <- c(treb_tree$tip.label[14], treb_tree$tip.label[15])
final_tree <- drop.tip(treb_tree, drop)
plot(final_tree)

################################################################################
# PART 8: SIMPLIFY TIP LABELS FOR CLEANER VISUALIZATION
################################################################################
# Extract just the OTU codes (first 3 characters) from full sequence names
# Special case: Label the outgroup as "Clade 'C'"

label2 <- c(final_tree$tip.label)
label2 <- substr(label2, 2, 4)  # Extract characters 2-4 from each label
label2 <- replace(label2, 13, "Clade 'C'")  # Replace 13th label with clade name

# Create data frame for renaming
d <- data.frame(label = final_tree$tip.label, label2 = label2)

# Apply new labels to tree
final_tree <- rename_taxa(final_tree, d, key = 1, value = 2)
plot(final_tree)

################################################################################
# PART 9: CREATE FIGURE 8A PHYLOGRAM (Tree Legend)
################################################################################
# Create phylogenetic tree visualization without branch lengths
# This serves as a legend showing relationships between OTUs
#
# Note: The flip and rotate operations organize the tree to match
# the color scheme in the stacked bar plot

# Initial tree plot
ggtree(final_tree, branch.length = "none") + 
  geom_tiplab() + 
  xlim(-10, 8) + 
  geom_text(aes(label = node))

# Store tree plot and manipulate topology for better visual arrangement
t <- ggtree(final_tree, branch.length = "none") + 
  geom_tiplab() + 
  xlim(-10, 8)

# Flip branches to match desired orientation
t <- flip(t, 20, 22)

# Rotate clade to improve readability
t <- ggtree::rotate(t, 22)

# Display final phylogram
t

################################################################################
# PART 10: ADD CLADE ASSIGNMENTS FOR STATISTICAL ANALYSIS
################################################################################
# Add a column indicating which major clade (I, S, or A) each OTU belongs to
# This allows calculation of clade-level summaries
#
# Clade determination based on first letter of OTU name:
# - I = impressa/gelatinosa clade
# - S = simplex/jamesii clade  
# - A = arboricola/gigantea clade

treb_long <- treb_long %>%
  mutate(clade = case_when(
    startsWith(OTU, "I") ~ "I",
    startsWith(OTU, "S") ~ "S",
    startsWith(OTU, "A") ~ "A"
  ))

################################################################################
# PART 11: CALCULATE CLADE-LEVEL PROPORTIONS
################################################################################
# Aggregate read counts to clade level and calculate proportions
# This addresses the question: Does the relative abundance of different
# Trebouxia clades change with warming?

# Sum sequences for each clade within each sample
treb_clade_sum <- aggregate((treb_long$sequences), 
                            by = list(treb_long$sample_id, 
                                     treb_long$individual, 
                                     treb_long$thallus_type, 
                                     treb_long$temp, 
                                     treb_long$CO2, 
                                     treb_long$clade), 
                            FUN = sum)

# Assign meaningful column names
colnames(treb_clade_sum) <- c("sample_id", "individual", "thallus_type", 
                              "temp", "CO2", "clade", "clade_sum")

# Calculate total sequences per sample (across all clades)
treb_clade_sum <- treb_clade_sum %>%
  group_by(sample_id) %>%
  mutate(ind_sum = sum(clade_sum))

# Calculate proportion of each clade within each sample
treb_clade_sum$clade_prop <- treb_clade_sum$clade_sum / treb_clade_sum$ind_sum

################################################################################
# PART 12: STATISTICAL ANALYSIS - Clade I vs. Temperature
################################################################################
# Test whether the proportion of clade 'I' increases with temperature
# Convert temperature to numeric for regression analysis

treb_clade_sum$temp <- as.numeric(treb_clade_sum$temp)

# Linear regression: proportion of clade I as a function of temperature
# Subset to only clade I samples
lm_prop_I <- lm(clade_prop ~ temp, 
                data = subset(treb_clade_sum, clade == "I"))

# View regression results (not printed in script, but check R console)
# summary(lm_prop_I)  # Uncomment to see detailed statistics

################################################################################
# PART 13: CREATE FIGURE 8B - Clade I Proportion vs. Temperature
################################################################################
# Scatter plot showing increase in clade 'I' proportion with warming
# Includes regression line from the model above
#
# Key result: Significant positive relationship (p = 0.0011, R² = 0.159)

ggplot(subset(treb_clade_sum, clade == "I"), 
       aes(x = factor(temp, levels = c("0", "2", "4", "6", "8", "10")), 
           y = clade_prop)) +
  geom_point() +
  theme_bw() +
  # Add regression line with coefficients from lm_prop_I model
  # intercept = 0.0007021, slope = 0.0196269
  geom_abline(intercept = 0.0007021, 
              slope = 0.0196269, 
              color = "#8397D3",  # Blue color matching clade I
              size = 2) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Temperature Differential °C") +
  ylab("Percent Clade 'I'")

################################################################################
# PART 14: ALPHA DIVERSITY ANALYSIS - Core vs. Tips
################################################################################
# Research question: Do core and tip regions of thalli differ in photobiont
# diversity? Hypothesis: Core regions maintain higher diversity as tips
# experience more environmental stress.

################################################################################
# PART 14A: Prepare Data for Metacoder Package
################################################################################
# Metacoder requires specific data structure with OTUs as rows

# Transpose OTU table so OTUs are rows and samples are columns
meta_otu <- as.data.frame(t(treb_otu))
colnames(meta_otu) <- meta_otu[1, ]  # First row becomes column names
meta_otu <- meta_otu[-1, ]           # Remove first row (now redundant)

# Add OTU_id column
meta_otu$OTU_id <- rownames(meta_otu)
meta_otu <- meta_otu %>%
  select(OTU_id, everything())  # Move OTU_id to first column
meta_otu <- as_tibble(meta_otu)

# Read and prepare taxonomy assignments
meta_tax <- read.csv("Trebouxia_tax_assignments.csv")
meta_tax <- as_tibble(meta_tax)

# Prepare sample metadata
meta_sample <- as_tibble(treb_meta)
meta_sample$sample_id <- as.character(meta_sample$sample_id)

# Ensure consistent data types for merging
meta_tax$OTU_id <- as.character(meta_tax$OTU_id)
meta_otu$OTU_id <- as.character(meta_otu$OTU_id)

################################################################################
# PART 14B: Merge OTU Table with Taxonomy
################################################################################
# Combine OTU counts with taxonomic classifications

meta_otu <- left_join(meta_otu, meta_tax,
                     by = c("OTU_id" = "OTU_id"))

################################################################################
# PART 14C: Create Taxmap Object for Rarefaction
################################################################################
# Taxmap is a specialized data structure for taxonomic data
# Allows rarefaction while maintaining taxonomic information

library(metacoder)

obj <- parse_tax_data(meta_otu,
                     class_cols = "taxonomy",  # Column with taxonomic string
                     class_sep = ";"           # Separator in taxonomy string
)

# Rename the data table for clarity
names(obj$data) <- "otu_counts"

################################################################################
# PART 14D: Rarefy to Even Sequencing Depth
################################################################################
# Rarefaction normalizes samples to the same number of reads
# Necessary because alpha diversity is sensitive to sequencing depth
#
# Rarefaction depth: 1790 reads (lowest sample total across dataset)
# This ensures fair comparison while maximizing sample retention

obj$data$otu_rarefied <- rarefy_obs(obj, "otu_counts", other_cols = TRUE)

################################################################################
# PART 14E: Calculate Shannon Diversity Index
################################################################################
# Shannon diversity index accounts for both richness (number of OTUs) and
# evenness (how evenly reads are distributed among OTUs)
# Higher values = more diverse photobiont communities

library(vegan)  # Ecological diversity analysis package

meta_sample$alpha.shannon <- diversity(obj$data$otu_rarefied[, meta_sample$sample_id],
                                      MARGIN = 2,        # Calculate by column (sample)
                                      index = "shannon")

################################################################################
# PART 14F: Prepare Paired Sample Data
################################################################################
# For Wilcoxon signed-rank test, need samples with BOTH core AND tip
# measurements from the same individual lichen
#
# This paired design accounts for individual-level variation

# Extract only relevant columns
paired <- meta_sample[, c(2, 3, 6)]  # individual, thallus_type, alpha.shannon

# Reshape to wide format: one row per individual, separate columns for core/tips
paired_wide <- pivot_wider(paired, 
                          names_from = thallus_type, 
                          values_from = alpha.shannon)

# Remove individuals without both measurements
paired_wide <- na.omit(paired_wide)

# Convert back to long format for statistical test
paired <- pivot_longer(paired_wide, 
                      cols = c(2, 3),           # Core and tips columns
                      names_to = "thallus_type", 
                      values_to = "alpha.shannon")

################################################################################
# PART 14G: Statistical Test - Core vs. Tips
################################################################################
# Wilcoxon signed-rank test: Non-parametric paired test
# Tests whether median diversity differs between core and tips
# Paired design controls for individual lichen variation
#
# Result: p-value = 0.01508 (significant difference)

wilcox <- wilcox.test(alpha.shannon ~ thallus_type, 
                     data = paired, 
                     paired = TRUE)

# Print test results
print(wilcox)

################################################################################
# PART 14H: Calculate Median Alpha Diversity
################################################################################
# Descriptive statistics for core and tip diversity

median(paired_wide$tips)  # Result: 0.1253272
median(paired_wide$core)  # Result: 0.3336968

# Interpretation: Core regions have ~2.7x higher median Shannon diversity
# than tips, suggesting photobiont communities in tips are more stressed
# or less diverse under the experimental conditions

################################################################################
# Script Complete
################################################################################
# Key outputs:
# - Figure 8A: Stacked bar plot of OTU composition + phylogram
# - Figure 8B: Regression plot of clade I proportion vs. temperature
# - Statistical results: 
#   * Clade I increases significantly with temperature
#   * Core regions have significantly higher diversity than tips
################################################################################
