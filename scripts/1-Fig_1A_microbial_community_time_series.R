# -----------------------------------------------------------------------
# Title: 1-Fig_1A_microbial_community_time_series.R
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: Create the time series plot on microbial community structure 
#              shown in Fig. 1A. This plot displays relative abundance 
#              at the higher taxonomic group level.
# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

## 1. Load Required Libraries
library(tidyverse) # Used for data wrangling (dplyr, tidyr) and plotting (ggplot2)

## 2. Load Metadata
# Information on ASVs (e.g., average relative abundance or taxonomic classification).
# This file provides the link between sequence IDs and higher taxonomic names or other attributes.
taxa_info <- read.csv("data/taxa_information.csv") 

## 3. Data Preparation: Taxonomic Aggregation
# Relative abundance per higher taxonomic group.
data <- read.csv("data/data_sequences_0.1_rel_ab_0.5_occ_binned_4_days_with_temperature.csv") %>%
  # Reshape from wide to long format (time-series style).
  gather(-date, key = "variable", value = "value") %>%
  # Join with taxa_info to get taxonomic labels for each sequence ID.
  left_join(taxa_info, by = c("variable" = "sequence_ID")) %>%
  # Quality Control: remove variables or ASVs that could not be classified.
  filter(!is.na(domain)) %>% 
  # Aggregate the relative abundance of individual ASVs into their 
  # respective higher taxonomic groups (e.g., 'Diatoms', 'Dinoflagellates').
  group_by(date, group, domain) %>%
  summarise(rel_ab = sum(value, na.rm = T), .groups = "drop") 

## 4. Define Visual Attributes
# Assemble attributes for plotting (e.g., color of higher taxonomic groups).
# We want the colors and sorting to be consistent based on abundance and domain.
plotting_attributes <- taxa_info %>%
  select(group, domain, rel_ab_mean, color) %>%
  filter(!is.na(domain)) %>% 
  group_by(group, domain, color) %>%
  summarise(rel_ab = sum(rel_ab_mean, na.rm = T), .groups = "drop") %>%
  # Order the groups by domain and then abundance to make the stack order 
  # of the bar chart logical.
  arrange(domain, rel_ab)

## 5. Formatting for ggplot2
# Apply the specific factor levels to ensure the bars are stacked 
# according to our 'plotting_attributes' order.
plot_data <- data %>%
  mutate(group = factor(group, levels = pull(plotting_attributes, "group")))

## 6. Calculate Secondary Axis (Total Abundance)
# Total relative abundance (normalized between 1 and 0).
# This represents the "abundance" of the total community over time.
scale_values <- function(x){(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))} 

total_abundance_scaled <- ungroup(data) %>%
  group_by(date) %>%
  summarise(rel_ab = sum(rel_ab), .groups = "drop") %>% 
  # Cleaning: replace zero abundance with NA to avoid plotting gaps.
  mutate(rel_ab = ifelse(rel_ab == 0, NA, rel_ab)) %>% 
  # Normalize values to a 0-1 range for the secondary y-axis.
  mutate(rel_ab_scaled = scale_values(rel_ab)) 



## 7. Generate Final Plot
# This plot creates a stacked bar chart (community composition) with a line overlay (total abundance).
ggplot(plot_data, aes(y = rel_ab, x = date(date))) +
  # Create the stacked bar chart. 'position="fill"' ensures the bars reach 100%.
  geom_bar(stat = "identity", aes(color = group, fill = group), position = "fill") +
  # Map specific colors to taxonomic groups from our metadata.
  scale_color_manual(values = pull(plotting_attributes, "color"), name = "") +
  scale_fill_manual(values = pull(plotting_attributes, "color"), name = "") +
  theme_classic() +
  ylab("Community composition (scaled)") +
  xlab("Time") +
  theme(legend.position = "top", 
        text = element_text(size = 12), 
        legend.text = element_text(size = 8), 
        legend.key.size = unit(0.3, "cm")) + 
  # Format x-axis dates.
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  # Define dual y-axis scaling.
  scale_y_continuous(expand = c(0, 0), 
                     sec.axis = sec_axis(~. * 1, name = "Total abundance (scaled)", breaks = seq(0, 1, 0.2)), 
                     breaks = seq(0, 1, 0.2)) +
  # Add the white line representing normalized total abundance.
  geom_line(data = total_abundance_scaled, aes(y = rel_ab_scaled), color = "white", size = .6) +
  # Clean up legend organization.
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

## 8. Save Figure
ggsave("plots/Fig_1A_microbial_community_time_series.pdf", width = 12, height = 3.5)