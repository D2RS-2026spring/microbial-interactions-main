# -----------------------------------------------------------------------
# Title: Keystones chorddiagram (Fig. 2A-C)
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: This script is used to recreate Fig. 3A showing keystone microbial interactions.
# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

# Required libraries
library(tidyverse) # awesome for data wrangling
library(data.table) # read and write large files
library(circlize) # needed for creating chorddiagrams
library(chorddiag)  # to create advanced chorddiagrams, can be installed using devtools::install_github("mattflor/chorddiag")

# -----------------------------------------------------------------------
# 1. Load required data
# -----------------------------------------------------------------------

## Interaction coefficients
mdr_smap_coefficients <- fread("model_out/MDR_smap_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv") %>%
  mutate(interaction=paste(target_2,"-",target_1))

## Strongest interactions
# Import the previously filtered subset of interactions that passed 
# both strength and temporal connectivity thresholds.
interactions <- read.csv("model_out/strongest_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv") %>%
  mutate(interaction=paste(target_2,"-",target_1))

interaction_vector <- pull(interactions,interaction) # strongest pairwise interactions

## Strongest interactions over time
interactions_time <- filter(mdr_smap_coefficients,interaction%in%interaction_vector) # only keep the strongest interactions

## Information on taxa
taxa_info <- read.csv("data/taxa_information.csv") %>%
  arrange(domain,group,sequence_ID) %>%
  select(sequence_ID,group,color)

## Keystone microbes
keystone_taxa <- read.csv("model_out/keystone_taxa_Scripps_Pier_ASVs_tp_2_28092024.csv")  %>%
  arrange(domain,group,sequence_ID)

keystone_taxa_vector <- pull(keystone_taxa,"sequence_ID")

## Create short labels for keystone microbes
# Shortening taxonomic group names (e.g., "Diatoms" to "Diat") and 
# adding numeric IDs to distinguish individual keystone ASVs within the same group.
group_short_IDs <- select(keystone_taxa,sequence_ID,domain,group) %>%
  arrange(domain,sequence_ID) %>%
  mutate(group_short=substr(group, start = 1, stop = 4)) %>%
  group_by(group_short) %>%
  mutate(group_ID=seq(1,n(),1)) %>%
  mutate(group_ID=as.character((group_ID))) %>%
  mutate(group_short_ID=paste0(group_short," ",group_ID)) %>%
  mutate(n_group=n()) %>%
  ungroup() %>%
  mutate(group_short_ID=ifelse(n_group==1,gsub(" 1","",group_short_ID),group_short_ID)) %>%
  select(sequence_ID,group_short_ID)


# -----------------------------------------------------------------------
# 1. Generate the chorddiagram
# -----------------------------------------------------------------------

# Assemble and format the interaction matrix for plotting
# Aggregation step. For the chord diagram, we want keystone ASVs to be 
# individual sectors, but all other microbes are grouped by their higher 
# taxonomic group to avoid visual clutter.
keystone_interaction_matrix <- interactions %>%
  ungroup() %>% # remove any prior grouping
  select(target_2,target_1,MDR_smap_coefficient) %>% # select relevant variables
  left_join(select(taxa_info,-color),by=c("target_2"="sequence_ID")) %>% # add group for target_2
  rename(group_2=group) %>%
  left_join(select(taxa_info,-color),by=c("target_1"="sequence_ID")) %>%  # add group for target_1
  rename(group_1=group) %>%
  mutate(target_1=ifelse(target_1%in%keystone_taxa_vector,target_1,group_1), # replace the target_1 name for all non-keystone (other) microbes by their higher taxonomic group
         target_2=ifelse(target_2%in%keystone_taxa_vector,target_2,group_1)) %>% # replace the target_2 name for all non-keystone (other) microbes by their higher taxonomic group
  group_by(target_2,target_1,group_1,group_2) %>%
  mutate(count=1) %>%
  summarise(MDR_smap_coefficient=sum(MDR_smap_coefficient),links=sum(count)) %>% # get summary statistics
  filter(target_2%in%keystone_taxa_vector) %>% # only keep keystone microbes as target 2
  left_join(group_short_IDs,by=c("target_2"="sequence_ID")) %>% # add the short label for keystone microbes (target_2)
  rename(label_2=group_short_ID) %>%
  left_join(group_short_IDs,by=c("target_1"="sequence_ID")) %>%   # add the short label for keystone microbes (target_1)
  rename(label_1=group_short_ID) %>%
  mutate(label_1=ifelse(is.na(label_1),group_1,label_1)) %>%  # assign group to the non-keystone microbes as a label
  ungroup()


## Chord colors (show the interaction strength and direction)
color_vector <- c('#003c30','#01665e','#35978f','#80cdc1','#c7eae5','#f5f5f5','#f6e8c3','#dfc27d','#bf812d','#8c510a','#543005')

## Calculate the limits for the colors and legend
limit <- filter(keystone_interaction_matrix) %>%
  filter(MDR_smap_coefficient==max(abs(MDR_smap_coefficient))) %>%
  pull(MDR_smap_coefficient)

## Generate the color for the chords (we will extract them from a ggplot object)
# Formatting: Using a dummy ggplot object to map the MDR S-map coefficients 
# to the diverging color palette for the chords.
p_2 <- ggplot(keystone_interaction_matrix,aes(y=MDR_smap_coefficient)) +
  geom_point(x=1,aes(fill=MDR_smap_coefficient),shape=21) +
  scale_fill_gradientn(colours=color_vector,values = scales::rescale(c(-limit, -0.01, 0, 0.01, limit)),
                       limits=c(-limit,limit),name="Interaction\nstrength\n ",
                       breaks=seq(-limit,limit,.1)) +
  theme(legend.position = "top")

g <- ggplot_build(p_2) # used to extract the color information for the chords

## Generate colors for the sectors
group_colors <- filter(unique(select(taxa_info,group,color)),group%in%pull(unique(select(keystone_interaction_matrix,group_1)))) # Extract the colors for all the groups in the group 1 column of the interaction matrix
color_grid <- c(pull(group_colors,"color"),pull(keystone_taxa,"color")) # create a vector of the group colors and keystone microbes colors

# Gather the data for plotting
plot_data <- bind_cols(keystone_interaction_matrix,g$data[[1]][1]) %>% # add color for chords to the interaction data
  select(label_1,label_2,MDR_smap_coefficient,fill) %>% # remove links, because this could confuse the plotting function (it will plot number of links instead of strength)
  arrange(label_2,-MDR_smap_coefficient) # this defines the chord order, separating between positive and negative interactions

## Final plot

## Note: This plot was slightly edited from it's original form for Fig. 3A to make the labels more readable and we manually added the number of links using Adobe Illustrator


pdf("plots/Fig_3A_keystone_chorddiagram.pdf") # save resulting plot as a pdf, be sure to use dev.off() to clear plot history, sometimes this function can give a bug

circos.clear() # remove any prior gap formatting with circos.par()
# Formatting: Customizing 'gap.after' to create large visual breaks between 
# keystone taxa and the aggregated taxonomic groups.
circos.par(gap.after = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,25,5,5,1,1,1,1,1,5,5,1,1,5,1,1,1,25)) # defines the gaps between segments, having larger gaps (5) between higher taxonomic groups of keystones and very large gaps (25) between keystone and other microbes

chordDiagramFromDataFrame(plot_data,
                          directional=1, # have a directional chord diagram
                          order=c(pull(group_colors,"group"),pull(group_short_IDs,"group_short_ID")), # define the order of sectors
                          # some formatting
                          transparency = 0.25,
                          annotationTrack = c("name","grid", "axis"),
                          direction.type = c("arrows", "diffHeight"), 
                          diffHeight  = -0.06,
                          annotationTrackHeight = c(0.05, 0.1),
                          link.arr.type = "big.arrow", 
                          link.sort = "asis", # keep the predefined order of chords
                          link.largest.ontop = TRUE, # plot larger links on top
                          col=plot_data$fill, # color for chords (interaction strength)
                          grid.col = color_grid) # color for sectors (higher taxonomic groups)

dev.off()

## Number of links per keystone microbe
link_summary <- keystone_interaction_matrix %>%
  group_by(label_2) %>%
  summarise(links=sum(links))
