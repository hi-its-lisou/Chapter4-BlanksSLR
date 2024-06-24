#load the libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)

#import the data
solitary_bee_bacteria_all<- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()
print(solitary_bee_bacteria_all)

solitary_bee_bacteria_uncontrolled <- solitary_bee_bacteria_all %>%
  filter(controlled_for_contamination == "No")

#Load in function for standard error
standard_error <- function(x) sd(x)/sqrt(length(x)) 

# list the genera from the top 5 bacterial taxa of solitary bees
genus_counts <- table(solitary_bee_bacteria_uncontrolled$genera)
genus_count_data <- data.frame(genera = names(genus_counts), Count = as.numeric(genus_counts))
genus_count_data
# Obtain top 10 reported genera
top_genera <- head(genus_count_data[order(-genus_count_data$Count), ], 10)
top_genera

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
solitary_bee_bacteria_uncontrolled$found_in_vector <- ifelse(solitary_bee_bacteria_uncontrolled$genera %in% vec, "Yes", "No")
# 2) whether it is a contaminant (1 if found in vector, 0 otherwise)
solitary_bee_bacteria_uncontrolled$contaminant <- ifelse(solitary_bee_bacteria_uncontrolled$genera %in% vec, 1, 0)
#write.table(solitary_bee_bacteria_uncontrolled,"solitary_bee_bacteria_uncontrolled.txt",sep= "\t",row.names=FALSE)
getcounts <- table(solitary_bee_bacteria_uncontrolled$found_in_vector)
getcounts #172 do not overlap and 71 overlap with common contaminants

#write_tsv(solitary_bee_bacteria_uncontrolled, "solitary_bee_bacteria_uncontrolled.tsv")

#Create table with proportions and proportional_biomass of overlapping taxa per bee microbiome description (bee ID)
proportion_contaminants_per_bee_id <- solitary_bee_bacteria_uncontrolled %>%
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

#write_tsv(proportion_contaminants_per_bee_id, "proportion_contaminants_per_bee_id")

# calculate the average proportion for overlapping taxa for controlled and uncontrolled studies per bee
average_proportion_per_group <- proportion_contaminants_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_proportion_contaminant = mean(proportion_contaminant, na.rm = TRUE),
    standard_error_proportion_contaminant = standard_error(proportion_contaminant)
  )
print(average_proportion_per_group)


# calculate the average proportional_biomass of overlapping taxa for controlled and uncontrolled studies per bee
average_contaminate_proportional_biomass_per_group <- proportion_contaminants_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_contaminant_proportional_biomass = mean(contaminant_proportional_biomass, na.rm = TRUE),
    standard_error_contaminant_proportional_biomass = standard_error(contaminant_proportional_biomass)
  )
print(average_contaminate_proportional_biomass_per_group)
