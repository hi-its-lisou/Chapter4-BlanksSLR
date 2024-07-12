#Load libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(lme4)

#Import and clean data
df<- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()

df$controlled_for_contamination <- as.factor(df$controlled_for_contamination)
df$bee_id <- as.character(df$bee_id)
df$paper_id <- as.character(df$paper_id)

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
df$found_in_vector <- ifelse(df$genus %in% vec, "Yes", "No")
df$contaminant <- ifelse(df$genus %in% vec, 1, 0)

#write.table(df,"df.txt",sep= "\t",row.names=FALSE)
getcounts <- table(df$found_in_vector)
getcounts #2377 do not overlap and 559 overlap with common contaminants

#Create table with prelative abundance of overlapping taxa per bee
overlap_per_bee_id <- df %>%
  group_by(bee_id, controlled_for_contamination, paper_id) %>%
  summarise(
    total_genus = n(),
    contaminant_genus = sum(contaminant),
    total_relative_abundance = sum(average_relative_abundance),
    overlap = sum(average_relative_abundance * contaminant),
    ) %>%
  ungroup()
print(overlap_per_bee_id)

# Create the box and whisker plot
ggplot(overlap_per_bee_id, 
       aes(x = controlled_for_contamination, 
           y = overlap)) +
  geom_boxplot(fill = "lightblue", 
               color = "black", 
               outlier.shape = 8, 
               outlier.size = 1.5) +
    labs(x = "Controlled for Contamination",
       y = "Average Relative Abundance of Overlapping taxa") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Three outliers identified - these are bees that have a high abudance of Acinetobacter
# Acinetobacter is a common nectar associate

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
    overlap2 = sum(average_relative_abundance * contaminant2)
  ) %>%
  ungroup()

#Plot new varibale data with old
#Pivot longer for the two overlap variables
tidy_overlap_per_bee_id <- overlap_per_bee_id %>%
  pivot_longer(cols = c(overlap, overlap2), names_to = "variable", values_to = "value")
#Plot it
(ggplot(tidy_overlap_per_bee_id, 
        aes(x=variable, 
            y=value, 
            fill=controlled_for_contamination)) +
 geom_boxplot(position=position_dodge(width=1),
              width=0.8,
              outlier.shape = 8, 
              outlier.size = 1.5) +
 labs(x="Controlled for Contamination",
      y = "Relative Abundance") +
 scale_x_discrete(labels = c("overlap" = "PPT", "overlap2" = "PPT -Acine")) +
 scale_fill_manual(values = c("#A72525", "355E3B"), name = "Contamination control") +  
 theme_light()+
 theme(axis.title = element_text(size = 12),
       axis.text = element_text(size = 12)))

ggsave("figures/boxplot.png", height=6, width=9)


### Is there a difference in the relative abundance of PCTs in papers that do and don't control contam? ###
# dataframe with the average relative abundance of PCT (including Acinetobacter) per bee 
standard_error <- function(x) sd(x)/sqrt(length(x)) 

average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_overlap = mean(overlap, na.rm = TRUE),
    standard_error_overlap = standard_error(overlap)
  )
print(average_overlap_per_group)

#Test the distribution of data
shapiro.test(overlap_per_bee_id$overlap)

# Not normally distributed, so using a generalised linear mixed model with random varibale set to paper_id
glmm <- glmer(overlap ~ controlled_for_contamination+(1|paper_id), 
              data = overlap_per_bee_id,
              family = binomial)
summary(glmm)

# dataframe with the average relative abundance of PCT (minus Acinetobacter) per bee 
average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_overlap2 = mean(overlap2, na.rm = TRUE),
    standard_error_overlap2 = standard_error(overlap2)
  )
print(average_overlap_per_group)

#Test the distribution of data
shapiro.test(overlap_per_bee_id$overlap) #now follows a normal distribution

# It is normally distributed, so using t-test
t.test(overlap2 ~ controlled_for_contamination, data=overlap_per_bee_id)


### Most reported bacteria across studies
# Group by controlled_for_contamination and genus, then count the occurrences
genus_counts <- df %>%
  group_by(controlled_for_contamination, genus, found_in_vector) %>%
  reframe(Count = n()) %>%
  ungroup()

# Separate the counts for "Yes" and "No" groups
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

