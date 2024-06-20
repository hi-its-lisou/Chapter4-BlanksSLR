#load libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)

#Load the data
SLR_DF <- read.delim("data/slr_metadata.txt", header = TRUE, sep = "\t") %>%
  clean_names()
SLR_DF
colnames(SLR_DF)

################################################################################
#Insect orders represented in the metadata
unique(SLR_DF$insect_order)
order_counts <- table(SLR_DF$insect_order)
order_data <- data.frame(insect_order = names(order_counts), Count = as.numeric(order_counts))

#Development stage used across the studies
unique(SLR_DF$development_stage)
Development_stage_counts <- table(SLR_DF$development_stage)
Development_stage_data <- data.frame(development_stage = names(Development_stage_counts), Count = as.numeric(Development_stage_counts))
Development_stage_data

#How many papers used a surface sterilisation step
Surface_sterilised_counts <- table(SLR_DF$surface_sterilized)
Surface_sterilised_data <- data.frame(surface_sterilized = names(Surface_sterilised_counts), Count = as.numeric(Surface_sterilised_counts))
Surface_sterilised_data
(161/245)*100

# The proportion of studies that sampled whole insects or dissected gut
sampled_counts <- table(SLR_DF$sampled)
Sampled_data <- data.frame(sampled = names(sampled_counts), Count = as.numeric(sampled_counts))
Sampled_data
(131/245)*100
(95/245)*100
((1+3+13+2)/245)*100

#count for the types of methods used to control for contamination
contam_control_how <- SLR_DF %>% summarise(count = n(), .by = c(how))
contam_control_how

################################################################################
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
  )                                     # Does not follow a normal distribution

# Mann-Whitney U test on whether use use negative controls and the number of references
wilcox.test(is_referenced_by_count ~ did_they_use_negative_controls, data = SLR_DF)
wilcox.test(is_referenced_by_count ~ was_there_a_comparison_with_their_samples_to_control_for_contamination=="Yes", data = SLR_DF)

################################################################################
#Get counts of publications per year
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

################################################################################
### Figure 3 - Part A ###
# Modified from Welsh & Eisenhofer (2024) The prevalence of controls in phyllosphere microbiome research: a methodological review.
#https://github.com/brady-welsh/2023_Controls_MR/blob/main/code/Figure_5.md

#Colour scheme 
colscheme <- c("Controlled contamination" = "#34A853", "Negative controls" = "#FF9900")

# Counts of studies per year
n_studies <- SLR_DF %>%
  group_by(year) %>%
  summarise(publications = n()) 
n_studies

#count of studies that used negative controls grouped by year
neg_prop <- SLR_DF %>%
  summarise(count = n(), .by = c(year, did_they_use_negative_controls)) %>%
  group_by(year) %>%
  mutate(proportion = count / sum(count)) %>%
  filter(did_they_use_negative_controls == "Yes")

#counts of studies that used sequenced blanks to control for contamination
controlcontam_prop <- SLR_DF %>%
  summarise(count = n(), .by = c(year, was_there_a_comparison_with_their_samples_to_control_for_contamination)) %>%
  group_by(year) %>%
  mutate(proportion = count / sum(count)) %>%
  filter(was_there_a_comparison_with_their_samples_to_control_for_contamination == "Yes")

# Add a factor level for "Number of Studies" in the control for contamination column for the legend
neg_prop$dummy <- "Negative controls"
controlcontam_prop$dummy <- "Controlled contamination"
n_studies$dummy <- "Number of Studies"

#Plot it
(neg_prop %>%
  ggplot(aes(x = as.factor(year), y = proportion * 100, 
             group = did_they_use_negative_controls,
             colour = "Negative controls")) +
  geom_line(linewidth = 2) +
  geom_line(data = controlcontam_prop, aes(x = as.factor(year), y = proportion * 100,
                                           group = was_there_a_comparison_with_their_samples_to_control_for_contamination,
                                           colour = "Controlled contamination"),
            linewidth = 2) +
  geom_point(data = controlcontam_prop, aes(x = as.factor(year), y = proportion * 100,
                                            group = was_there_a_comparison_with_their_samples_to_control_for_contamination,
                                            colour = "Controlled contamination"),
             size = 4) +
  geom_line(data = n_studies, 
            aes(x = as.factor(year), y = publications, group = dummy, colour = dummy, linetype = dummy), 
            linewidth = 1,
            linetype = "dashed") +
  scale_y_continuous(
    name = "Percentage of studies",
    sec.axis = sec_axis(~ ., name = "Number of Studies"),
    labels = scales::label_percent(scale = 1)
  ) +
  scale_color_manual(values = c(colscheme, "Number of Studies" = "grey")) +
  scale_linetype_manual(values = c("Number of Studies" = "dashed", "Negative controls" = "solid", "Controlled contamination" = "solid")) +
  theme_classic() +
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

ggsave("figures/Figure3a.png", height=9, width=12)

################################################################################
### Figure 3 - Part B ###
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
interaction_labels <- c("No Negative Controls",
                        "Used negative controls, but did not control for contamination",
                        "Used negative controls and controlled for contamination")

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
      "No Negative Controls" = "red",
      "Used negative controls, but did not control for contamination" = "111111",
      "Used negative controls and controlled for contamination" = "355E3B"
    )) +
  theme (legend.title = element_blank(),
         legend.text = element_text(size = 12),
         legend.position = "bottom"))

ggsave("figures/Figure3b.png", height=9, width=12)

################################################################################
# Proportion of studies that acknowledge relic DNA
acks1 <- table(SLR_DF$relic_dna)
acks_data1 <- data.frame(relic_dna = names(acks1), Count = as.numeric(acks1))
acks_data1
(8/245)*100

# Proportion of studies that acknowledge transient bacteria
acks2 <- table(SLR_DF$transient_bacteria)
acks_data2 <- data.frame(transient_bacteria = names(acks2), Count = as.numeric(acks2))
acks_data2
(50/245)*100

# Proportion of studies that acknowledge mitochondria or chloroplasts
acks3 <- table(SLR_DF$animal_plant_dna)
acks_data3 <- data.frame(animal_plant_dna = names(acks3), Count = as.numeric(acks3))
acks_data3
(85/245)*100

