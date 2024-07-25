#Load libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(lme4)
library(brms)
library(betareg)




#Import and clean data
df<- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()

df$controlled_for_contamination <- as.factor(df$controlled_for_contamination)
df$paper_id <- as.character(df$paper_id)
df$bee_id <- as.character(df$bee_id)

#df <- filter(df, average_relative_abundance >=0.01)

#Bees in the data 
unique(df$bee_genus)
unique(df$bee_species)

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

### Most reported bacteria across studies ###
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



# Assuming 'df' is your initial dataframe
overlap_per_bee_id <- df %>%
  group_by(bee_id, controlled_for_contamination, paper_id) %>%
  summarise(
    total_genus = n(),
    contaminant_genus = sum(contaminant),
    total_relative_abundance = sum(average_relative_abundance),
    overlap = sum(average_relative_abundance * contaminant)
  ) %>%
  ungroup()

print(overlap_per_bee_id)

# Function to calculate standard error
standard_error <- function(x) sd(x) / sqrt(length(x)) 

# Dataframe with the average relative abundance of PCT per bee
average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_overlap = mean(overlap, na.rm = TRUE),
    standard_error_overlap = standard_error(overlap),
    n = n(),
    lower_ci = average_overlap - qt(0.95, df = n-1) * standard_error_overlap,
    upper_ci = average_overlap + qt(0.95, df = n-1) * standard_error_overlap
  )

print(average_overlap_per_group)

# Plotting the boxplot with error bars for 95% confidence intervals
p3 <- (ggplot(overlap_per_bee_id, 
              aes(x=controlled_for_contamination, 
                  y=overlap, 
                  fill=controlled_for_contamination)) +
         geom_boxplot(position=position_dodge(width=2),
                      width=0.8,
                      outlier.shape = 8, 
                      outlier.size = 1.5) +
         geom_errorbar(data = average_overlap_per_group, 
                       aes(x = controlled_for_contamination, 
                           ymin = lower_ci, 
                           ymax = upper_ci, 
                           y = average_overlap), 
                       width = 0.2,
                       color = "#56B4E9", 
                       size = 1.5) +
         geom_point(size=3) +
                  labs(x="Contamination control",
              y = "Relative Abundance") +
         scale_x_discrete(labels = c("No" = "Uncontrolled (n=25)", "Yes" = "Controlled (n=15)")) +
         scale_fill_manual(values = c("#E69F00", "#009E73"), name = "Contamination control") +  
         theme_light() +
         theme( legend.position = "none",
                axis.title = element_text(size = 12),
                axis.text = element_text(size = 12))
)
p3


ggsave("figures/Figure 4.png", height=6, width=9)





#Test the distribution of data
shapiro.test(overlap_per_bee_id$overlap)

# Not normally distributed, so using a generalised linear mixed model with random varibale set to paper_id
wilcox.test(overlap ~ controlled_for_contamination,
            data = overlap_per_bee_id)


glmm <- glmer(overlap ~ controlled_for_contamination+(1|paper_id), 
              data = overlap_per_bee_id,
              family = binomial)
summary(glmm)


beta <- betareg(overlap ~ controlled_for_contamination,
        data = overlap_per_bee_id)

summary(beta)

### Supplementary material ###
#removing Acinetobacter as a potential contaminant

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

#ggsave("figures/S5.png", height=6, width=9)

# dataframe with the average relative abundance of PCT (minus Acinetobacter) per bee 
average_overlap_per_group <- overlap_per_bee_id %>%
  group_by(controlled_for_contamination) %>%
  summarise(
    average_overlap2 = mean(overlap2, na.rm = TRUE),
    standard_error_overlap2 = standard_error(overlap2)
  )
print(average_overlap_per_group)

#Test the distribution of data
shapiro.test(overlap_per_bee_id$overlap2)

wilcox.test(overlap2 ~ controlled_for_contamination, 
            data = overlap_per_bee_id)

glmm <- glmer(overlap2 ~ controlled_for_contamination+(1|paper_id), 
              data = overlap_per_bee_id,
              family = binomial)
summary(glmm)





######https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8265668/##################################3
library(vegan)

df<- read.delim("data/adhoc_data.txt", header = TRUE, sep = "\t") %>%
  clean_names()
df$controlled_for_contamination <- as.factor(df$controlled_for_contamination)
df$paper_id <- as.factor(df$paper_id)
df$bee_id <- as.character(df$bee_id)

# Pivot the data to create a community matrix
community_matrix <- df %>%
  select(bee_id, genus, average_relative_abundance) %>%
  spread(key = genus, value = average_relative_abundance)
# Replace missing values with zeros
community_matrix[is.na(community_matrix)] <- 0
# Set bee_id as row names
rownames(community_matrix) <- community_matrix$bee_id
community_matrix <- community_matrix[ , -1]
# Add metadata for plotting
metadata <- df %>%
  select(bee_id, controlled_for_contamination, paper_id, bee_genus, bee_species) %>%
  distinct()

# Compute Bray-Curtis dissimilarity matrix
dissimilarity_matrix <- vegdist(community_matrix, method = "bray")


# Perform NMDS
nmds <- metaMDS(dissimilarity_matrix, k = 2)
# Merge metadata with NMDS scores
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$bee_id <- rownames(nmds_scores)
nmds_data <- merge(nmds_scores, metadata, by = "bee_id")

# Plot NMDS
manual_nmds <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = controlled_for_contamination)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = c("#A72525", "355E3B"), name = "Contamination control") +
    labs(title = "Bray-Curtis beta diversity per bee",
       x = "NMDS1",
       y = "NMDS2")


# Perform PERMANOVA
adonis_result <- adonis2(dissimilarity_matrix ~ controlled_for_contamination, data = metadata)
print(adonis_result)

# Perform SIMPER analysis
simper_result <- simper(community_matrix, metadata$controlled_for_contamination)
summary(simper_result)
#Flavobacterium, Stenotrophomonas, Atopobium




# Calculate Shannon diversity index
shannon_diversity <- diversity(community_matrix, index = "shannon")
# Calculate Evenness
evenness <- shannon_diversity / log(specnumber(community_matrix))
# Combine alpha diversity indices with metadata
alpha_diversity_df <- data.frame(
  bee_id = rownames(community_matrix),
  shannon_diversity = shannon_diversity,
  evenness = evenness
)

# Join with metadata
alpha_diversity_df <- alpha_diversity_df %>%
  left_join(metadata, by = "bee_id")


# Plot Shannon Diversity Index
ggplot(alpha_diversity_df, aes(x = controlled_for_contamination, y = shannon_diversity, fill = controlled_for_contamination)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Contamination Control",
       x = "Controlled for Contamination",
       y = "Shannon Diversity Index") +
  scale_fill_manual(values = c("#A72525", "355E3B"), name = "Contamination control") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot Evenness
ggplot(alpha_diversity_df, aes(x = controlled_for_contamination, y = evenness, fill = controlled_for_contamination)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Evenness by Contamination Control",
       x = "Controlled for Contamination",
       y = "Evenness") +
  scale_fill_manual(values = c("#A72525", "355E3B"), name = "Contamination control") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




