#load the libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)

#import the data
solitary_bee_bacteria_all<- read.delim("Data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()
print(solitary_bee_bacteria_all)

# list the genera from the top 5 bacterial taxa of solitary bees
genus_counts <- table(solitary_bee_bacteria_all$genera)
genus_count_data <- data.frame(genera = names(genus_counts), Count = as.numeric(genus_counts))
genus_count_data
# Obtain top 10 reported genera
top_genera <- head(genus_count_data[order(-genus_count_data$Count), ], 10)
top_genera

#Create a vector for common contaminants taken from Eisenhofer et al. 2019

#vec<- c("Actinomyces", "Corynebacterium", "Arthrobacter", "Rothia", 
#        "Propionibacterium", "Atopobium", "Sediminibacterium", 
#        "Porphyromonas", "Prevotella", "Chryseobacterium", 
#        "Capnocytophaga", "Chryseobacterium", "Flavobacterium", 
#        "Pedobacter", "UnclassifiedTM7", "Bacillus", "Geobacillus", 
#        "Brevibacillus", "Paenibacillus", "Staphylococcus", 
#        "Abiotrophia", "Granulicatella", "Enterococcus", 
#        "Lactobacillus", "Streptococcus", "Clostridium", 
#        "Coprococcus", "Anaerococcus", "Dialister", "Megasphaera", 
#        "Veillonella", "Fusobacterium", "Leptotrichia", 
#        "Brevundimonas", "Afipia", "Bradyrhizobium", "Devosia", 
#        "Methylobacterium", "Mesorhizobium", "Phyllobacterium", 
#        "Rhizobium", "Methylobacterium", "Phyllobacterium", 
#        "Roseomonas", "Novosphingobium", "Sphingobium", 
#        "Sphingomonas", "Achromobacter", "Burkholderia", 
#        "Acidovorax", "Comamonas", "Curvibacter", "Pelomonas",
#        "Cupriavidus", "Duganella", "Herbaspirillum", 
#        "Janthinobacterium", "Massilia", "Oxalobacter", 
#        "Ralstonia", "Leptothrix", "kingella", "Neisseria", 
#        "Escherichia", "Haemophilus", "Acinetobacter", 
#        "Enhydrobacter", "Pseudomonas", "Stenotrophomonas", 
#        "Xanthomonas")

#removing Acinetobacter and Lactobacillus as potential contaminants

vec<- c("Actinomyces", "Corynebacterium", "Arthrobacter", "Rothia", 
        "Propionibacterium", "Atopobium", "Sediminibacterium", 
        "Porphyromonas", "Prevotella", "Chryseobacterium", 
        "Capnocytophaga", "Chryseobacterium", "Flavobacterium", 
        "Pedobacter", "UnclassifiedTM7", "Bacillus", "Geobacillus", 
        "Brevibacillus", "Paenibacillus", "Staphylococcus", 
        "Abiotrophia", "Granulicatella", "Enterococcus", 
     "Streptococcus", "Clostridium", 
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

### Add a new column with vector
# 1) whether the bacteria taxa is "found_in_vector" with values "Yes" or "No"
solitary_bee_bacteria_all$found_in_vector <- ifelse(solitary_bee_bacteria_all$genera %in% vec, "Yes", "No")
# 2) whether it is a contaminant (1 if found in vector, 0 otherwise)
solitary_bee_bacteria_all$contaminant <- ifelse(solitary_bee_bacteria_all$genera %in% vec, 1, 0)
#write.table(solitary_bee_bacteria_all,"solitary_bee_bacteria_all.txt",sep= "\t",row.names=FALSE)
getcounts <- table(solitary_bee_bacteria_all$found_in_vector)
getcounts #172 do not overlap and 71 overlap with common contaminants

overlap_in_controlled_studies <- solitary_bee_bacteria_all %>%
  group_by(bee_id) %>%
  filter(controlled_for_contamination == "Yes") %>%
  pull(contaminant)

overlap_in_uncontrolled_studies <- solitary_bee_bacteria_all %>%
  group_by(bee_id) %>%
  filter(controlled_for_contamination == "No") %>%
  pull(contaminant)

# Perform the t-test
t_test_result <- t.test(overlap_in_controlled_studies, overlap_in_uncontrolled_studies)

print(t_test_result)

#Create table with proportions and proportional_biomass of overlapping taxa per bee microbiome description (bee ID)
proportion_contaminants_per_bee_id <- solitary_bee_bacteria_all %>%
  group_by(bee_id, controlled_for_contamination) %>%
  summarise(
    total_genera = n(),
    contaminant_genera = sum(contaminant),
    proportion_contaminant = mean(contaminant),
    total_proportional_biomass_of_top_genera = sum(proportional_biomass),
    contaminant_proportional_biomass = sum(proportional_biomass * contaminant),
    proportion_proportional_biomass_contaminant = sum(proportional_biomass * contaminant) / sum(proportional_biomass) * 100
    ) %>%
  ungroup()
print(proportion_contaminants_per_bee_id)

# calculate the average proportion for overlapping taxa for controlled and uncontrolled studies per bee
average_proportion_per_group <- proportion_contaminants_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(average_proportion_contaminant = mean(proportion_contaminant))
print(average_proportion_per_group)

# calculate the average proportional_biomass of overlapping taxa for controlled and uncontrolled studies per bee
average_contaminate_proportional_biomass_per_group <- proportion_contaminants_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(average_contaminant_proportional_biomass = mean(contaminant_proportional_biomass))
print(average_contaminate_proportional_biomass_per_group)



##############################################################################################################################################################################################
subsetted_df <- solitary_bee_bacteria_all %>%
 filter(proportional_biomass > 10)

overlap_in_controlled_studies <- subsetted_df %>%
  group_by(bee_id) %>%
  filter(controlled_for_contamination == "Yes") %>%
  pull(contaminant)

overlap_in_uncontrolled_studies <- subsetted_df %>%
  group_by(bee_id) %>%
  filter(controlled_for_contamination == "No") %>%
  pull(contaminant)

# Perform the t-test
t_test_result <- t.test(overlap_in_controlled_studies, overlap_in_uncontrolled_studies)

print(t_test_result)

#Create table with proportions and proportional_biomass of overlapping taxa per bee microbiome description (bee ID)
proportion_contaminants_per_bee_id <- subsetted_df %>%
  group_by(bee_id, controlled_for_contamination) %>%
  summarise(
    total_genera = n(),
    contaminant_genera = sum(contaminant),
    proportion_contaminant = mean(contaminant),
    total_proportional_biomass_of_top_genera = sum(proportional_biomass),
    contaminant_proportional_biomass = sum(proportional_biomass * contaminant),
    proportion_proportional_biomass_contaminant = sum(proportional_biomass * contaminant) / sum(proportional_biomass) * 100
  ) %>%
  ungroup()

print(proportion_contaminants_per_bee_id)

# calculate the average proportion for overlapping taxa for controlled and uncontrolled studies per bee
average_proportion_per_group <- proportion_contaminants_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(average_proportion_contaminant = mean(proportion_contaminant))
print(average_proportion_per_group)

# calculate the average proportional_biomass of overlapping taxa for controlled and uncontrolled studies per bee
average_contaminate_proportional_biomass_per_group <- proportion_contaminants_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(average_contaminant_proportional_biomass = mean(contaminant_proportional_biomass))
print(average_contaminate_proportional_biomass_per_group)

# t-test for proportion of overlapping taxa 
overlap_in_controlled_studies <- proportion_contaminants_per_bee_id %>%
  group_by(bee_id) %>%
  filter(controlled_for_contamination == "Yes") %>%
  pull(contaminant_proportional_biomass)

overlap_in_uncontrolled_studies <- proportion_contaminants_per_bee_id %>%
  group_by(bee_id) %>%
  filter(controlled_for_contamination == "No") %>%
  pull(contaminant_proportional_biomass)

# Assess the distribution of the data
shapiro.test(proportion_contaminants_per_bee_id$contaminant_proportional_biomass)

#The data is not normally distributed, so going with a wilcox test
mann_whitney_test_result <- wilcox.test(overlap_in_controlled_studies, overlap_in_uncontrolled_studies)
#met with error: Warning message:
In wilcox.test.default(overlap_in_controlled_studies, overlap_in_uncontrolled_studies) :
  cannot compute exact p-value with ties
