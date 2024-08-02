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

ggsave("figures/Figure3.png", height=5, width=7.5)

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
  pivot_longer(cols = c(overlap, overlap_vec2), names_to = "vector", values_to = "average_PCT_RA")

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

ggsave("figures/S6.png", height=6, width=9)


