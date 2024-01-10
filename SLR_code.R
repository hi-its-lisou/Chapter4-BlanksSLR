#importing the dataframe
setwd("D:/Research/PhD/Data/SLR")
SLR_DF <- read.delim("slr_metadata.txt", header = TRUE, sep = "\t")
SLR_DF

#load libraries
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)

#cleaning the names
clean_names(data)
SLR_DF

#Insect orders represented in the metadata
unique(SLR_DF$Insect_order)

order_counts <- table(SLR_DF$Insect_order)
order_data <- data.frame(Insect_order = names(order_counts), Count = as.numeric(order_counts))

ggplot(order_data, aes(x = "", y = Count, fill = Insect_order)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = , fill = "Insect Order") +
  scale_fill_brewer(palette="Set3") +
  theme(text = element_text(size = 20))


################################################################################
#Development stage used
unique(SLR_DF$Development_stage)

Development_stage_counts <- table(SLR_DF$Development_stage)
Development_stage_data <- data.frame(Development_stage = names(Development_stage_counts), Count = as.numeric(Development_stage_counts))
Development_stage_data

ggplot(Development_stage_data, aes(x = "", y = Count, fill = Development_stage)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = "", fill = "Development stage sampled") +
  scale_fill_brewer(palette = "Set1") +
  geom_text(aes(label = paste0(round(Count / sum(Count) * 100), "%")), 
            position = position_stack(vjust = 0.5)
            , size= 6)+
  theme(text = element_text(size=20))

################################################################################
#Surface sterilisation
Surface_sterilized_counts <- table(SLR_DF$Surface_sterilized)
Surface_sterilized_data <- data.frame(Surface_sterilized = names(Surface_sterilized_counts), Count = as.numeric(Surface_sterilized_counts))
Surface_sterilized_data
(167/254)*100

################################################################################
#Get counts of publications per year
SLR_DF %>%
  group_by(Year) %>%
  summarise(publications = n()) %>%
ggplot(aes(x = Year, y = publications, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal()+
  stat_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "blue") +
  geom_text(aes(x = max(Year), y = max(publications), 
                label = sprintf("R^2 = %.2f, p = %.3f", 
                                cor(Year, publications)^2, 
                                summary(lm(publications ~ Year))$coefficients[8])),
            hjust = 1, vjust = 1, color = "blue")

#Negative control use through the years
unique(SLR_DF$Year)
unique(SLR_DF$Did_they_use_negative_controls)

summary_df <- SLR_DF %>%
  count(Year, Did_they_use_negative_controls) %>%
  arrange(Year)
summary_df

ggplot(summary_df, aes(x = as.factor(Year), y = n, fill = Did_they_use_negative_controls)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Publication Counts Over Years",
       x = "Year",
       y = "Number of Publications",
       fill = "Used Negative Controls") +
  scale_fill_manual(values = c("No" = "darkred", "Yes" = "blue")) +
  theme_minimal()

#Split the studies that do use negative controls into those that control for contamination and those that don't
#new summary dataframe to split negative control use further according to year
summary_df <- SLR_DF %>%
  count(Year, Did_they_use_negative_controls,
        Contamination_Control = Was_there_a_comparison_with_their_samples_to_control_for_contamination) %>%
  arrange(Year)


# the interaction variable needs to be a character
summary_df$interaction_level <- paste(
  summary_df$Did_they_use_negative_controls,
  summary_df$Contamination_Control,
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

# Convert the interaction variable to factor with thee specific  levels and labels
summary_df$interaction_level <- factor(
  summary_df$interaction_level,
  levels = interaction_levels,
  labels = interaction_labels
)

# Stacked bar chart to show the proportion of papers that used negative controls AND those that control for contamination
ggplot(summary_df, aes(x = as.factor(Year), y = n, fill = interaction_level)) +
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
    ),
    guide = guide_legend(title = "Negative Controls")
  ) +
  theme_minimal()

