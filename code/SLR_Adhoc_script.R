#load the libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)

#import the data
df<- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()

#Load in function for standard error
standard_error <- function(x) sd(x)/sqrt(length(x)) 

# list the genus from the top 5 bacterial taxa of solitary bees
genus_counts <- table(df$genus)
genus_count_data <- data.frame(genus = names(genus_counts), Count = as.numeric(genus_counts))
genus_count_data
# Obtain top 10 reported genus
top_genus <- head(genus_count_data[order(-genus_count_data$Count), ], 10)
top_genus

#Create a vector for common contaminants taken from Eisenhofer et al. 2019
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

### Add a new column with vector
# 1) whether the bacteria taxa is "found_in_vector" with values "Yes" or "No"
df$found_in_vector <- ifelse(df$genus %in% vec, "Yes", "No")
# 2) whether it is a contaminant (1 if found in vector, 0 otherwise)
df$contaminant <- ifelse(df$genus %in% vec, 1, 0)
#write.table(df,"df.txt",sep= "\t",row.names=FALSE)
getcounts <- table(df$found_in_vector)
getcounts #2377 do not overlap and 559 overlap with common contaminants

#Create table with proportions and proportional_biomass of overlapping taxa per bee microbiome description (bee ID)
overlap_per_bee_id <- df %>%
  group_by(bee_id, controlled_for_contamination) %>%
  summarise(
    total_genus = n(),
    contaminant_genus = sum(contaminant),
    proportion_contaminant = mean(contaminant),
    overlap = sum(average_relative_abundance * contaminant),
    ) %>%
  ungroup()
print(overlap_per_bee_id)

# calculate the average relative abundance that overlaps with common contaminants per bee
average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_overlap = mean(overlap, na.rm = TRUE),
    standard_error_overlap = standard_error(overlap)
  )
print(average_overlap_per_group)

#test for normality
shapiro.test(overlap_per_bee_id$overlap)

# Not normally distributed, so perform Wilcoxon rank-sum test
# Is overlap significantly different between contamination control methods
wilcox_test <- wilcox.test(overlap ~ controlled_for_contamination, 
                           data = overlap_per_bee_id)
# Print the result of the Wilcoxon rank-sum test
print(wilcox_test) #not significant

# Create the box and whisker plot
ggplot(overlap_per_bee_id, 
       aes(x = controlled_for_contamination, 
           y = overlap)) +
  geom_boxplot(fill = "lightblue", 
               color = "black", 
               outlier.shape = 8, 
               outlier.size = 1.5) +
  geom_jitter(width = 0.2, 
              color = "red", 
              alpha = 0.5) +
  labs(title = "Contaminant Average Relative Abundance by Contamination Control",
       x = "Controlled for Contamination",
       y = "Average Relative Abundance Per Bee microbiome of Overlapping taxa") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Outliers identified - these are bees that have a high abudance of Acinetobacter
# Acinetobacter is a common nectar associate

# Remove the three bees that are outliers
DF <- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names() %>%
  filter(!bee_id %in% c("16", "30", "37"))

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

DF$found_in_vector <- ifelse(DF$genus %in% vec, "Yes", "No")
DF$contaminant <- ifelse(DF$genus %in% vec, 1, 0)

overlap_per_bee_id <- DF %>%
  group_by(bee_id, controlled_for_contamination) %>%
  summarise(
    total_genus = n(),
    contaminant_genus = sum(contaminant),
    proportion_contaminant = mean(contaminant),
    overlap = sum(average_relative_abundance * contaminant),
  ) %>%
  ungroup()
print(overlap_per_bee_id)

average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(average_overlap = mean(overlap),
            standard_error_overlap = standard_error(overlap))
print(average_overlap_per_group)

shapiro.test(overlap_per_bee_id$overlap) #now follows a normal distribution

# Is overlap significantly different between contamination control methods
t.test(overlap ~ controlled_for_contamination, data=overlap_per_bee_id)

# Create the box and whisker plot
ggplot(overlap_per_bee_id, aes(x = controlled_for_contamination, y = overlap)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.shape = 8, outlier.size = 1.5) +
  geom_jitter(width = 0.2, color = "red", alpha = 0.5) +
  labs(
    title = "Contaminant Average Relative Abundance by Contamination Control",
    x = "Controlled for Contamination",
    y = "Contaminant Average Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
