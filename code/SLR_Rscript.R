##############
### Set up ###
##############

#load libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)

#Load and clean the metadataset 
SLR_DF <- read.delim("data/slr_metadata.txt", header = TRUE, sep = "\t") %>%
  clean_names()
SLR_DF
colnames(SLR_DF)

########################
### Research subject ### 
########################

#Insect orders represented in the metadata
order_counts <- table(SLR_DF$insect_order)
order_data <- data.frame(insect_order = names(order_counts), Count = as.numeric(order_counts))
order_data <- order_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
order_data

# Target microbial community of insect
SLR_DF <- SLR_DF %>%
  mutate(target_community_of_microbiome = if_else(
    target_community_of_microbiome %in% c("Body", "Gut"), 
    target_community_of_microbiome, 
    "Other")) # Changing the labels of the target community, so that if it's not Body, Gut, it's 'Other'

target_microbiome_data <- table(SLR_DF$target_community_of_microbiome)
target_microbiome_data <- data.frame(target_community_of_microbiome = names(target_microbiome_data), Count = as.numeric(target_microbiome_data))
target_microbiome_data <- target_microbiome_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
target_microbiome_data

#Development stage used across the studies
Development_stage_data <- table(SLR_DF$development_stage)
Development_stage_data <- data.frame(development_stage = names(Development_stage_data), Count = as.numeric(Development_stage_data))
Development_stage_data <- Development_stage_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
Development_stage_data

#Target region of the 16S rRNA gene
region_V_data <- table(SLR_DF$x16s_r_rna_gene_target_region)
region_V_data <- data.frame(x16s_r_rna_gene_target_region = names(region_V_data), Count = as.numeric(region_V_data))
region_V_data <- region_V_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
region_V_data

#################################
### Specimen handling methods ### 
#################################

# The proportion of studies that sampled whole insects or dissected gut
unique(SLR_DF$sampled)
SLR_DF <- SLR_DF %>%
  mutate(sampled = if_else(
    sampled %in% c("Gut", "Whole"), 
    sampled, 
    "Other"))

sample_method_data <- table(SLR_DF$sampled)
sample_method_data <- data.frame(sampled = names(sample_method_data), Count = as.numeric(sample_method_data))
sample_method_data <- sample_method_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
sample_method_data

#How many papers used a surface sterilisation step
Surface_sterilised_data <- table(SLR_DF$surface_sterilized)
Surface_sterilised_data <- data.frame(surface_sterilized = names(Surface_sterilised_data), Count = as.numeric(Surface_sterilised_data))
Surface_sterilised_data <- Surface_sterilised_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
Surface_sterilised_data

################################
### Use of negative controls ### 
################################
negative_controls <- table(SLR_DF$did_they_use_negative_controls)
negative_controls <- data.frame(did_they_use_negative_controls = names(negative_controls), Count = as.numeric(negative_controls))
negative_controls <- negative_controls %>% mutate(Percentage = (Count / sum(Count)) * 100)
negative_controls

controlled_contam <- table(SLR_DF$was_there_a_comparison_with_their_samples_to_control_for_contamination)
controlled_contam <- data.frame(was_there_a_comparison_with_their_samples_to_control_for_contamination = names(controlled_contam), Count = as.numeric(controlled_contam))
controlled_contam <- controlled_contam %>% mutate(Percentage = (Count / sum(Count)) * 100)
controlled_contam

### Figure 3 - Part A
# Modified from Welsh & Eisenhofer (2024) The prevalence of controls in phyllosphere microbiome research: a methodological review. https://github.com/brady-welsh/2023_Controls_MR/blob/main/code/Figure_5.md

#Colour scheme 
colscheme <- c("Controlled for contamination using blanks" = "355E3B", "Used blanks, but did not control contamination" = "#014D72")

# Counts of Publications per year
n_Publications <- SLR_DF %>%
  group_by(year) %>%
  summarise(publications = n()) 
n_Publications

####
# Dataframe for the number of studies that use negative controls each year (but did not control for contam)
neg_prop <- SLR_DF %>%
  mutate(n_yes = ifelse(did_they_use_negative_controls == "Yes", 1, 0)) %>%
  group_by(year) %>%
  summarise(n_studies = length(year),
            n_yes = sum(n_yes),
            .groups = "drop") %>% 
  rowwise() %>%
  mutate(proportion = n_yes / n_studies,
         did_they_use_negative_controls = "Yes")

# Dataframe for the number of studies that use negative controls to control contamination each year
controlcontam_prop <- SLR_DF %>%
  mutate(n_yes = ifelse(was_there_a_comparison_with_their_samples_to_control_for_contamination == "Yes", 1, 0)) %>%
  group_by(year) %>%
  summarise(n_studies = length(year),
            n_yes = sum(n_yes, na.rm = TRUE),
            .groups = "drop") %>% 
  rowwise() %>%
  mutate(proportion = n_yes / n_studies,
         was_there_a_comparison_with_their_samples_to_control_for_contamination = "Yes")

# Add a factor level for "Number of Publications" in the control for contamination column for the legend
neg_prop$dummy <- "Used blanks, but did not control contamination"
controlcontam_prop$dummy <- "Controlled for contamination using blanks"
n_Publications$dummy <- "Number of Publications"

#Plot it
(neg_prop %>%
  ggplot(aes(x = as.factor(year), y = proportion * 100, 
             group = did_they_use_negative_controls,
             colour = "Used blanks, but did not control contamination")) +
   geom_line(linewidth = 2) +
   geom_line(
     data = controlcontam_prop, 
     aes(x = as.factor(year), 
         y = proportion * 100,
         group = was_there_a_comparison_with_their_samples_to_control_for_contamination,
         colour = "Controlled for contamination using blanks"),
     linewidth = 2,
     position = position_jitterdodge(0.25)) +
   geom_line(
     data = n_Publications, 
     aes(x = as.factor(year), 
         y = publications, 
         group = dummy, colour = dummy, linetype = dummy),
     linewidth = 1,
     linetype = "dashed"
     ) +
   scale_y_continuous(
    name = "Percentage of studies",
    sec.axis = sec_axis(~ ., name = "Number of Publications"),
    labels = scales::label_percent(scale = 1)
    ) +
  scale_color_manual(
    values = c(colscheme, "Number of Publications" = "grey")
    ) +
  scale_linetype_manual(
    values = c("Number of Publications" = "dashed", 
               "Used blanks, but did not control contamination" = "solid", 
               "Controlled for contamination using blanks" = "solid")
    ) +
 theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.y.right = element_line(color = "grey", linetype = "dashed", linewidth = 1),
        axis.ticks.y.right = element_line(color = "grey"),
        axis.text.y.right = element_text(color = "grey"),
        axis.title.y.right = element_text(color = "grey")
  ) +
  labs(x = "Year", y = "Percentage of studies"))

ggsave("figures/Figure3a.png", height=6, width=9)

### Figure 3 - Part B
#new dataframe to split negative control use further according to year
summary_df <- SLR_DF %>%
  count(year, did_they_use_negative_controls,
        contamination_Control = was_there_a_comparison_with_their_samples_to_control_for_contamination) %>%
  arrange(year)

# Create an interaction variable by combining control variables into a single string
summary_df$interaction_level <- paste(
  summary_df$did_they_use_negative_controls,
  summary_df$contamination_Control,
  sep = "."
)

# Define the levels for the interaction variable
interaction_levels <- c( "No.NA",
                         "Yes.No",
                         "Yes.Yes")
# Defning the labels for the interactions
interaction_labels <- c("Did not report using a negative control",
                        "Used blanks, but did not control contamination",
                        "Controlled for contamination using blanks")

# Convert the interaction variable to factor with the specific levels and labels
summary_df$interaction_level <- factor(
  summary_df$interaction_level,
  levels = interaction_levels,
  labels = interaction_labels
)

#Plot it
(ggplot(summary_df, aes(x = as.factor(year), y = n, fill = interaction_level)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = " ",
    x = "Year",
    y = "Number of Publications",
    fill = "Comparison with Samples for Contamination Control"
  ) +
  scale_fill_manual(
    values = c(
      "Did not report using a negative control" = "#6B1818",
      "Used blanks, but did not control contamination" = "#014D72",
      "Controlled for contamination using blanks" = "355E3B"
    )) +
 theme_bw() +
  theme (legend.title = element_blank(),
         legend.text = element_text(size = 12),
         legend.position = "bottom", 
         axis.title = element_text(size = 14),
         axis.text = element_text(size = 14),
         axis.text.x = element_text(angle = 45, hjust = 1)
         ))

ggsave("figures/Figure3b.png", height=6, width=9)

#Plot the number of publications per year to obtain p value
SLR_DF %>%
  group_by(year) %>%
  summarise(publications = n()) %>%
  ggplot(aes(x = year, y = publications, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal()+
  stat_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "blue") +
  geom_text(aes(x = max(year), y = max(publications), 
                label = sprintf("R^2 = %.2f, p = %.3f", 
                                cor(year, publications)^2, 
                                summary(lm(publications ~ year))$coefficients[8])),
            hjust = 1, vjust = 1, color = "blue")

#Check the distribution of citation data with a histogram
ggplot(SLR_DF, aes(x = is_referenced_by_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Reference Counts",
    x = "Reference Count",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    title = element_text(size = 16)
  )# Does not follow a normal distribution

# Mann-Whitney U test on whether use use negative controls and the number of references
wilcox.test(is_referenced_by_count ~ did_they_use_negative_controls, data = SLR_DF)
wilcox.test(is_referenced_by_count ~ was_there_a_comparison_with_their_samples_to_control_for_contamination=="Yes", data = SLR_DF)

##########################
### Limit of detection ### 
##########################

SLR_DF <- SLR_DF %>%
  mutate(do_they_determine_the_limit_of_detection = if_else(
    do_they_determine_the_limit_of_detection %in% c("Yes"), 
    do_they_determine_the_limit_of_detection, 
    "No")) #Needed to merge NAs with No

LoD_data <- table(SLR_DF$do_they_determine_the_limit_of_detection)
LoD_data <- data.frame(do_they_determine_the_limit_of_detection = names(LoD_data), Count = as.numeric(LoD_data))
LoD_data <- LoD_data %>% mutate(Percentage = (Count / sum(Count)) * 100)
LoD_data

#########################################################################
### Acknowledging alternative, non-symbiotic sources of microbial DNA ### 
#########################################################################
# Proportion of studies that acknowledge mitochondria or chloroplasts
acks_mito_chloro <- table(SLR_DF$animal_plant_dna)
acks_mito_chloro <- data.frame(animal_plant_dna = names(acks_mito_chloro), Count = as.numeric(acks_mito_chloro))
acks_mito_chloro <- acks_mito_chloro %>% mutate(Percentage = (Count / sum(Count)) * 100)
acks_mito_chloro

# Proportion of studies that acknowledge transient bacteria
ack_transient <- table(SLR_DF$transient_bacteria)
ack_transient <- data.frame(transient_bacteria = names(ack_transient), Count = as.numeric(ack_transient))
ack_transient <- ack_transient %>% mutate(Percentage = (Count / sum(Count)) * 100)
ack_transient

# Proportion of studies that acknowledge relic DNA
ack_relic <- table(SLR_DF$relic_dna)
ack_relic <- data.frame(relic_dna = names(ack_relic), Count = as.numeric(ack_relic))
ack_relic <- ack_relic %>% mutate(Percentage = (Count / sum(Count)) * 100)
ack_relic
