#Load libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(lme4)
library(vegan)
library(microbiome)
library(phyloseq)
library(ANCOMBC)
library(RColorBrewer)
library(tidyr)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(ggforce)
library(patchwork)
library(ggrepel)

### Import, clean and structure data ### 
df<- read.delim("data/adhoc_data2.txt", header = TRUE, sep = "\t") %>%
  clean_names()

df$primer <- NA  # Initialize the primer column with NA
df$primer[df$paper_id == 1] <- "V5-V6"
df$primer[df$paper_id == 3] <- "V4"
df$primer[df$paper_id == 4] <- "V6-V8"
df$primer[df$paper_id == 6] <- "V1-V3"
df$primer[df$paper_id == 7] <- "V5-V6"
df$primer[df$paper_id == 8] <- "V4"
df$primer[df$paper_id == 9] <- "V5-V6"
df$primer[df$paper_id == 10] <- "V4"
df$primer[df$paper_id == 11] <- "V5-V6"
df$primer[df$paper_id == 12] <- "V5-V6"

df$controlled_for_contamination <- as.factor(df$controlled_for_contamination)
df$paper_id <- as.character(df$paper_id)
df$bee_id <- as.character(df$bee_id)

### Creating a dummy phyloseq plot to observe relative abundances between groups ###
# Prepare the OTU table
community_matrix <- df %>% # Create a community matrix
  select(bee_id, genus, average_relative_abundance) %>%
  spread(key = genus, value = average_relative_abundance)

community_matrix[is.na(community_matrix)] <- 0
rownames(community_matrix) <- community_matrix$bee_id # Set bee_id as row names
community_matrix <- community_matrix[, -1]

# Convert OTU_table to phyloseq object
OTU <- otu_table(as.matrix(community_matrix), taxa_are_rows = FALSE)

# Prepare the metadata
metadata <- df %>%
  select(bee_id, controlled_for_contamination, paper_id, bee_genus, bee_species, primer) %>%
  distinct() %>%
  mutate(combined_label = paste(paper_id, bee_genus, bee_species, sep = "_"))

rownames(metadata) <- metadata$bee_id # Set bee_id as row names

#Convert metadata to phyloseq object
sample_data <- sample_data(metadata)
rownames(sample_data) <- sample_data$bee_id

df$combined_label <- paste(df$paper_id, 
                           df$bee_genus, 
                           df$bee_species, 
                           sep = "_")

# Prepare the tax table
tax_data <- data.frame(
  OTU = colnames(community_matrix),
  Kingdom = "Bacteria",
  Phylum = "NA",
  Class = "NA",
  Order = "NA",
  Family = "NA",
  Genus = colnames(community_matrix)
)
# Convert to matrix and set row names
tax_matrix <- as.matrix(tax_data[, -1])
rownames(tax_matrix) <- tax_data$OTU # Set genera as row names

# Convert tax data into a phyloseq object
TAX <- tax_table(tax_matrix)# Turn into the tax_table object

# Combine together to create the dummy phyloseq object
ps <- phyloseq(OTU, sample_data, TAX)

###########################################################################
### Plot the relative abundance of bees from the two contamination control groups ###
# Prepare genera
ps_Genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
ps_melted <- psmelt(ps_Genus_ra)

#ps_melted <- ps_melted %>%
#  mutate(Genus = ifelse(Genus %in% PCTs, "PCT", Genus))

# Identify the Top 20 Most Abundant Genera (excluding "PCT")
top_20_genera_filtered <- ps_melted %>%
  filter(Genus != "Other") %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  top_n(20, total_abundance) %>%
  pull(Genus)

ps_melted <- ps_melted %>%
  mutate(Genus = ifelse(Genus %in% top_20_genera_filtered, Genus, "Other"))

ps_melted_top_20_filtered <- ps_melted %>%
  filter(Genus %in% c(top_20_genera_filtered, "Other"))

ps_melted_top_20_filtered <- ps_melted_top_20_filtered %>%
  mutate(Genus = factor(Genus, levels = Abundance$Genus[order(-genus_abundance_filtered$total_abundance)]))

genus_abundance_filtered <- ps_melted_top_20_filtered %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance))

ps_melted_top_20_filtered <- ps_melted_top_20_filtered %>%
  mutate(Genus = factor(Genus, levels = rev(genus_abundance_filtered$Genus[order(-genus_abundance_filtered$total_abundance)])))

#Colour palette specifically for PCTs 
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter4/Chapter4_BlanksSLR/data",
  "colour_list_ch4.csv") %>% read_csv
my_palette <- colours_df$colours
names(my_palette) <- colours_df$Genus
my_palette

(plot_top_20_filtered <- ps_melted_top_20_filtered %>%
    ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~controlled_for_contamination + paper_id + primer, scales = "free") +
    theme(axis.text.x = element_text(angle = 70, size = 12, face= 'bold', hjust=1),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.justification = "center") +
    labs(x = "", y = "Relative Abundance"))

ggsave("figures/adhoc_RelativeAbundances.png", height=9, width=12)


### Plot the relative abundance of bees from the two contamination control groups ###
# Prepare genera
ps_Genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
ps_melted <- psmelt(ps_Genus_ra)

#Vector for PCTs
PCT_vec <- read_delim("data/PCT_vec.txt", 
                      delim = "\n", 
                      col_names = "PCT")
PCTs <- PCT_vec$PCT

#Make all taxa that are not a PCT "other"
ps_melted <- ps_melted %>%
  mutate(Genus = ifelse(Genus %in% PCTs, Genus, "Other"))

#custom oder for the plot to have PCTs on top 
order <- read_delim("data/custom_order.txt", 
                    delim = "\n", 
                    col_names = "PCT")
custom_order <- order$PCT

ps_melted$Genus <- factor(ps_melted$Genus, levels = custom_order)

# Create the relative abundance plot
PCTS <- ps_melted %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~controlled_for_contamination + paper_id + primer, scales = "free") +
  theme(axis.text.x = element_text(angle = 70, size = 12, face= 'bold', hjust=1),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.justification = "center") +
  labs(x = "", y = "Relative Abundance")

ggsave("figures/adhoc_RelativeAbundances_PCTs.png", height=10, width=12)



### Alpha diversity ###
#Compare alpha diversity between groups
alpha_diversity <- alpha(ps, index = "Shannon")
metadata <- meta(ps)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

alpha_diversity_metadata %>%
  ggplot(aes(x = controlled_for_contamination, y = diversity_shannon , colour = controlled_for_contamination)) +
  geom_boxplot() +
  geom_point(size = 3)

### Beta diversity ###
# Calculate Bray-Curtis dissimilarity and perform PCoA
bray_curtis_dist <- distance(ps, method = "bray")
ordination <- ordinate(ps, method = "PCoA", distance = bray_curtis_dist)

# Create plot
ordination_plot <- plot_ordination(physeq = ps,
                                   ordination = ordination,
                                   color = "controlled_for_contamination",
                                   axes = c(1, 2)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_text_repel(aes(label = combined_label),  
                  max.overlaps = 15,
                  size = 3) +
  theme_minimal() +
  stat_ellipse(geom = "polygon", type="norm", alpha=0) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ordination_plot
ggsave("figures/ordination_plot.png", height=6, width=9)


# Adonis statistic to determine whether communities are significantly different
#Refactor metadata
md_combined_subset <- as(sample_data(ps), "data.frame")
#Run adonis with 9999 permutations
uw_adonis <- adonis2(bray_curtis_dist ~ controlled_for_contamination, 
                     data = as(sample_data(ps), "data.frame"), 
                     permutations = 9999)
uw_adonis


