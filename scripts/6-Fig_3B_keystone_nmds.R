# -----------------------------------------------------------------------
# Title: Keystone NMDS (Fig. 3B)
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: This script is used to recreate Fig. 3B the NMDS plot on keystone interactions.
# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

## Required libraries
library(tidyverse) # awesome for data wrangling
library(vegan) # required for the NMDS analysis
library(data.table)

# -----------------------------------------------------------------------
# 1. Load required data
# -----------------------------------------------------------------------

## Interaction coefficients
mdr_smap_coefficients <- fread("model_out/MDR_smap_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv") %>%
  mutate(interaction=paste(target_2,"-",target_1))

## Strongest interactions
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

## Short labels for keystone microbes
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
# 2. NMDS analysis
# -----------------------------------------------------------------------
# Non-metric Multidimensional Scaling (NMDS) is used here to visualize 
# the similarity of "interaction profiles" between keystone taxa. 
# Taxa that interact similarly with the community will cluster together.

set.seed(1234) # random seed for reproducibility in the NMDS iterative process

## Date vector
date_vector <- unique(interactions_time$time)

## NMDS Loop
# We perform a separate ordination for each time point to see how 
# keystone roles shift through time.
for(i in 1:length(date_vector)) { # do the NMDS for each point in time
  
  data_i <- filter(interactions_time,time==date_vector[i]) # get interactions for time point i
  
  ## Shape the data for NMDS
  # We shift coefficients to be positive to satisfy the Bray-Curtis requirements.
  nmds_interactions <-  data_i %>%
    filter(target_2%in%keystone_taxa_vector) %>%
    select(target_2,target_1,MDR_smap_coefficient) %>%
    mutate(MDR_smap_coefficient=MDR_smap_coefficient+abs(min(MDR_smap_coefficient))) %>%
    spread(target_1,MDR_smap_coefficient,fill=0)
  
  nmds_target_2_vector <- pull(nmds_interactions,target_2)
  
  ## Turn data into matrix for nmds function
  nmds_data <- select(nmds_interactions,-target_2) %>%
    as.matrix()
  
  ## Run nmds function using Bray-Curtis dissimilarity
  # metaMDS automatically handles square-root transformations and Wisconsin 
  # double standardization if needed.
  nmds = metaMDS(nmds_data, distance = "bray")
  
  ## Save the NMDS output for time point i
  # Stress measures how well the 2D plot represents the high-dimensional distances.
  nmds_scores_dots_i <- data.frame(scores(nmds)$sites,target=nmds_target_2_vector) %>%
    mutate(time=date_vector[i],
           stress=nmds[["stress"]])
  
  if(i==1){nmds_scores_dots <- nmds_scores_dots_i}else{nmds_scores_dots<-rbind(nmds_scores_dots,nmds_scores_dots_i)}
  
}

## Estimate hull-values (polygons encompassing data points per taxon)
# chull() identifies the outermost points for each keystone taxon 
# across all time points to visualize their "Interaction Niche" area.


for (i in 1:length(keystone_taxa_vector)){
  
  target_i <- keystone_taxa_vector[i]
  
  hull_i <- nmds_scores_dots[nmds_scores_dots$target == target_i, ][chull(nmds_scores_dots[nmds_scores_dots$target == target_i, c("NMDS1", "NMDS2")]), ] 
  
  if(i==1){hull.data <- hull_i} else {hull.data <- rbind(hull.data,hull_i)}  #combine hull
}

## Get text for NMDS plot
nmds_scores_dots_text <- group_by(nmds_scores_dots,target) %>%
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2)) %>%
  left_join(keystone_taxa,by=c("target"="sequence_ID")) %>%
  left_join(group_short_IDs,,by=c("target"="sequence_ID"))

## Assemble the data for plotting
nmds_plot_data <- mutate(nmds_scores_dots,target=factor(target,levels=pull(keystone_taxa,"sequence_ID"))) %>%
  # Quality Control: stress < 0.2 is the standard threshold for a reliable NMDS ordination.
  filter(stress<0.2)  %>% # exclude values where stress is too high --> nmds becomes unreliable
  mutate(ID=seq(1,nrow(.),1))

## Final NMDS plot (Fig. 2B)
ggplot(nmds_plot_data,aes(x=NMDS1,y=NMDS2)) +
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=target,group=target),alpha=0.30) + # add the convex hulls
  geom_point(aes(fill=target,color=target),shape=21,size=2,alpha=.58) +
  geom_label(data=nmds_scores_dots_text,aes(label=group_short_ID),hjust=.5, vjust=.5) + 
  scale_fill_manual(values=pull(keystone_taxa,"color")) +
  scale_color_manual(values=pull(keystone_taxa,"color")) +
  theme_classic()  +
  theme(legend.text.align = 0,plot.background=element_blank(),
        text = element_text(size=12),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1), # size changed to linewidth for ggplot2 3.4.0+ compatibility
        legend.position = "none",
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines')) +
  labs(x="NMDS1", y = "NMDS2") +
  xlim(-0.7,0.7) +
  ylim(-0.7,0.7) 

ggsave("plots/Fig_3B_keystone_nmds.pdf",width=6,height=6)