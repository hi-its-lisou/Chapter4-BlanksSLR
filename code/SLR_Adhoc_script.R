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
df<- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()

df$controlled_for_contamination <- as.factor(df$controlled_for_contamination)
df$paper_id <- as.character(df$paper_id)
df$bee_id <- as.character(df$bee_id)

### Bees represented in the data ### 
unique(df$bee_genus)
unique(df$bee_species)

### Common contaminants in microbiome studies ###
# Create a vector for common contaminants taken from Eisenhofer et al. 2019
vec<- c("Actinomyces", "Corynebacterium", "Arthrobacter", "Rothia", 
        "Propionibacterium", "Atopobium", "Sediminibacterium", 
        "Porphyromonas", "Prevotella", "Chryseobacterium", 
        "Capnocytophaga", "Chryseobacterium", "Flavobacterium", 
        "Pedobacter", "UnclassifiedTM7", "Bacillus", "Geobacillus", 
        "Brevibacillus", "Paenibacillus", "Staphylococcus", 
        "Abiotrophia", "Granulicatella", "Enterococcus", 
        "Lactobacillus", "Streptococcus", "Clostridium", 
        "Coprococcus", "Anaerococcus", "Dialister", "Megasphaera", 
        "Veillonella", "Fusobacterium", "Leptotrichia", 
        "Brevundimonas", "Afipia", "Bradyrhizobium", "Devosia", 
        "Methylobacterium", "Mesorhizobium", "Phyllobacterium", 
        "Rhizobium", "Methylobacterium", "Phyllobacterium", 
        "Roseomonas", "Novosphingobium", "Sphingobium", 
        "Sphingomonas", "Achromobacter", "Burkholderia", 
        "Acidovorax", "Comamonas", "Curvibacter", "Pelomonas",
        "Cupriavidus", "Duganella", "Herbaspirillum", 
        "Janthinobacterium", "Massilia", "Oxalobacter", 
        "Ralstonia", "Leptothrix", "kingella", "Neisseria", 
        "Escherichia", "Haemophilus", "Acinetobacter", 
        "Enhydrobacter", "Pseudomonas", "Stenotrophomonas", 
        "Xanthomonas")

# Add a new column to relative abundance with vector
df$found_in_vector <- ifelse(df$genus %in% vec, "Yes", "No")
df$contaminant <- ifelse(df$genus %in% vec, 1, 0)

getcounts <- table(df$found_in_vector)
getcounts #2374 do not overlap and 559 overlap with common contaminants

### Most commonly reported bacteria across studies ###
genus_counts <- df %>%
  group_by(controlled_for_contamination, genus, found_in_vector) %>%
  reframe(Count = n()) %>%
  ungroup()

# The top 10 most commonly reported bacteria in controlled and uncontrolled groups
genus_counts_yes <- genus_counts %>%
  filter(controlled_for_contamination == "Yes") %>%
  arrange(desc(Count)) %>%
  head(10)
print(genus_counts_yes)

genus_counts_no <- genus_counts %>%
  filter(controlled_for_contamination == "No") %>%
  arrange(desc(Count)) %>%
  head(10)
print(genus_counts_no)

### Relative abundance of Potential Contaminant Taxa ### 
# Dataframe for relative abundance of PCTs per bee. Overlap is the variable term used to describe the relative abundance of PCTs.
overlap_per_bee_id <- df %>%
  group_by(bee_id, controlled_for_contamination, paper_id) %>%
  summarise(
    total_genus = n(),
    contaminant_genus = sum(contaminant),
    total_relative_abundance = sum(average_relative_abundance), #because it the averaged relative abundance per taxa, not always equal to 100%
    overlap = sum(average_relative_abundance * contaminant)
  ) %>%
  ungroup()

print(overlap_per_bee_id)

# Test the distribution of data
shapiro.test(overlap_per_bee_id$overlap) #Not normally distributed, 

# Test for significance 
wilcox.test(overlap ~ controlled_for_contamination,
            data = overlap_per_bee_id) # p-value = 0.7825

glmm <- glmer(overlap ~ controlled_for_contamination+(1|paper_id), 
              data = overlap_per_bee_id,
              family = binomial) 
summary(glmm) #p-value = 0.216


# Function to calculate standard error
standard_error <- function(x) sd(x) / sqrt(length(x))

# Average relative abundance of PCTs per bee in controlled and uncontrolled groups, including SE the 95% CIs for plotting.
average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_overlap = mean(overlap, na.rm = TRUE),
    standard_error_overlap = standard_error(overlap),
    n = n(),
    lower_ci = average_overlap - qt(0.975, df = n-1) * standard_error_overlap,
    upper_ci = average_overlap + qt(0.975, df = n-1) * standard_error_overlap
  )

print(average_overlap_per_group)
#A tibble: 2 x 6
# controlled_for_contamination    average_overlap   standard_error_overlap       n       lower_ci  upper_ci
# <fct>                              <dbl>                  <dbl>               <int>    <dbl>     <dbl>
#  1 No                              0.333                 0.0617                25      0.228     0.439
#  2 Yes                             0.243                 0.0605                15      0.136     0.349

# Plotting the boxplot with error bars for 95% confidence intervals
(figure3 <- (ggplot(overlap_per_bee_id, 
              aes(x=controlled_for_contamination, 
                  y=overlap, 
                  fill=controlled_for_contamination)) +
         geom_boxplot(width=0.8,
                      outlier.shape = 4, 
                      outlier.size = 3) +
         geom_errorbar(data = average_overlap_per_group, 
                       aes(x = controlled_for_contamination, 
                           ymin = lower_ci, 
                           ymax = upper_ci, 
                           y = average_overlap), 
                       width = 0.2,
                       color = "dark grey", 
                       size = 1.5) +
          geom_jitter(width=0.1, alpha = 0.9) +
          labs(x="",
              y = "Relative Abundance of PCTs per Bee") +
         scale_x_discrete(labels = c("No" = "Uncontrolled for contamination(n=25)", "Yes" = "Controlled for contamination (n=15)")) +
         scale_fill_manual(values = c("#A72525", "#009E73")) +  
         theme_light() +
         theme( legend.position = "none",
                axis.title = element_text(size = 12),
                axis.text = element_text(size = 12))))

ggsave("figures/Figure 3.png", height=5, width=7.5)

### Supplementary material ###
#Two outliers identified using the interquartile method are bees with high amounts of Acinetobacter, a common nectar associate
#Removing Acinetobacter as a potential contaminant

#Create second vector without Acinetobacter as a potential contaminant
vec2<- c("Actinomyces", "Corynebacterium", "Arthrobacter", "Rothia", 
         "Propionibacterium", "Atopobium", "Sediminibacterium", 
         "Porphyromonas", "Prevotella", "Chryseobacterium", 
         "Capnocytophaga", "Chryseobacterium", "Flavobacterium", 
         "Pedobacter", "UnclassifiedTM7", "Bacillus", "Geobacillus", 
         "Brevibacillus", "Paenibacillus", "Staphylococcus", 
         "Abiotrophia", "Granulicatella", "Enterococcus", 
         "Lactobacillus", "Streptococcus", "Clostridium", 
         "Coprococcus", "Anaerococcus", "Dialister", "Megasphaera", 
         "Veillonella", "Fusobacterium", "Leptotrichia", 
         "Brevundimonas", "Afipia", "Bradyrhizobium", "Devosia", 
         "Methylobacterium", "Mesorhizobium", "Phyllobacterium", 
         "Rhizobium", "Methylobacterium", "Phyllobacterium", 
         "Roseomonas", "Novosphingobium", "Sphingobium", 
         "Sphingomonas", "Achromobacter", "Burkholderia", 
         "Acidovorax", "Comamonas", "Curvibacter", "Pelomonas",
         "Cupriavidus", "Duganella", "Herbaspirillum", 
         "Janthinobacterium", "Massilia", "Oxalobacter", 
         "Ralstonia", "Leptothrix", "kingella", "Neisseria", 
         "Escherichia", "Haemophilus",  
         "Enhydrobacter", "Pseudomonas", "Stenotrophomonas", 
         "Xanthomonas")

#rerun code to add values from new vector
df$found_in_vector2 <- ifelse(df$genus %in% vec2, "Yes", "No")
df$contaminant2 <- ifelse(df$genus %in% vec2, 1, 0)

overlap_per_bee_id <- df %>%
  group_by(bee_id, controlled_for_contamination, paper_id) %>%
  summarise(
    total_genus = n(),
    contaminant_genus = sum(contaminant),
    total_relative_abundance = sum(average_relative_abundance),
    overlap = sum(average_relative_abundance * contaminant),
    contaminant_genus2 = sum(contaminant2),
    overlap_vec2 = sum(average_relative_abundance * contaminant2)
  ) %>%
  ungroup()

print(overlap_per_bee_id)

#Pivot longer for the two overlap vectors
tidy_overlap_per_bee_id <- overlap_per_bee_id %>%
  pivot_longer(cols = c(overlap, overlap_vec2), names_to = "vector", values_to = "average_PCT_RA") %>%
  mutate(dummy=CI)
  

#Function to summarise the data and obtain 95% confidence intervals
calc_ci <- function(x) {
  n <- length(x)
  mean_x <- mean(x, na.rm = TRUE)
  stderr <- sd(x, na.rm = TRUE) / sqrt(n)
  error <- qt(0.975, df = n-1) * stderr
  c(lower = mean_x - error, upper = mean_x + error)
}

#Function to summarise the data and obtain 95% confidence intervals
average_overlap_per_group <- tidy_overlap_per_bee_id %>%
  group_by(controlled_for_contamination, vector) %>%
  summarise(
    n = n(),
    average_PCTs_RA = mean(average_PCT_RA, na.rm = TRUE),
    lower_ci = calc_ci(average_PCT_RA)[1],
    upper_ci = calc_ci(average_PCT_RA)[2]) %>%
  mutate(dummy = "CI") # To add to the boxplot

print(average_overlap_per_group)

#Plot the relative abundance of PCTs for the two common contaminant vectors 
(ggplot(tidy_overlap_per_bee_id, 
        aes(x=interaction(controlled_for_contamination, vector), 
            y=average_PCT_RA, 
            fill=interaction(controlled_for_contamination, vector))) +
    geom_boxplot(position=position_dodge(width=1),
                 width=0.8,
                 outlier.shape = 4, 
                 outlier.size = 3) +
    geom_jitter(width=0.1) +
    geom_errorbar(data = average_overlap_per_group, 
                  aes(x=interaction(controlled_for_contamination, vector),
                      ymin = lower_ci, 
                      ymax = upper_ci, 
                      y = average_PCTs_RA), 
                  width = 0.2,
                  color = "dark grey", 
                  size = 1.5) +
    geom_line(data = average_overlap_per_group, 
              aes(x = interaction(controlled_for_contamination, vector),
                  y = average_PCTs_RA,
                  color = dummy),
              size = 1.5) +
    labs(x="",
         y = "Relative Abundance per Bee") +
    scale_x_discrete(labels = c("No.overlap" = "", 
                                "Yes.overlap" = "", 
                                "No.overlap_vec2" = "", 
                                "Yes.overlap_vec2" = "")) +
    scale_fill_manual(values = c("#A72525", "#009E73", "#FFA16E", "#B8F2BB"),
                      name = "Contamination control",
                      labels = c("No.overlap" = "All PCT, uncontrolled", 
                                 "Yes.overlap" = "All PCT, controlled", 
                                 "No.overlap_vec2" = "PCT -Acinetobacter, uncontrolled", 
                                 "Yes.overlap_vec2" = "PCT -Acinetobacter, controlled")) +
    scale_color_manual(values = c("dark grey"),
                       name = "Confidence Interval",
                       labels = c("CI" = "95% CI")) +
    theme_light() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12)))

#ggsave("figures/S5.png", height=6, width=9)


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
  select(bee_id, controlled_for_contamination, paper_id, bee_genus, bee_species) %>%
  distinct() %>%
  mutate(combined_label = paste(paper_id, bee_genus, bee_species, sep = "_"))

rownames(metadata) <- metadata$bee_id # Set bee_id as row names

ps$combined_label <- paste(sample_data_df$paper_id, 
                           sample_data_df$bee_genus, 
                           sample_data_df$bee_species, 
                           sep = "_")

#Convert metadata to phyloseq object
sample_data <- sample_data(metadata)
rownames(sample_data) <- sample_data$bee_id


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

#Colour palette specifically for PCTs 
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter4/Chapter4_BlanksSLR/data",
  "colour_list_ch4.csv") %>% read_csv
my_palette <- colours_df$colours
names(my_palette) <- colours_df$Genus
my_palette

#custom oder for the plot to have PCTs on top 
order <- read_delim("data/custom_order.txt", 
                      delim = "\n", 
                      col_names = "PCT")
custom_order <- order$PCT

ps_melted$Genus <- factor(ps_melted$Genus, levels = custom_order)

# Create the relative abundance plot
(plot <- ps_melted %>%
    ggplot(aes(x = combined_label, y = Abundance, fill = Genus)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_wrap(~ controlled_for_contamination, 
               scales = "free") +
    theme(axis.text.x = element_text(angle = 70, 
                                     size = 12, 
                                     face= 'bold',
                                     hjust=1),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.justification = "center") +
    labs(x = "",
         y = "Relative Abundance"))


ggsave("figures/FigureS7.png", height=15, width=10)

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

# Adonis statistic to determine whether communities are significantly different
#Refactor metadata
md_combined_subset <- as(sample_data(ps), "data.frame")
#Run adonis with 9999 permutations
uw_adonis <- adonis2(bray_curtis_dist ~ controlled_for_contamination, 
                     data = as(sample_data(ps), "data.frame"), 
                     permutations = 9999)
uw_adonis