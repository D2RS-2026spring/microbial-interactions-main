# -----------------------------------------------------------------------
# Title: Microbial interactions (Fig. 2A-C)
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: This script is used to recreate Fig. 2A-C on microbial interactions.
# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

## Required libraries
library(tidyverse) # for data wrangling
library(lubridate) # awesome to work with dates
library(data.table) # read and write large data files
library(ggpubr) # add statistics to plots (e.g., p-values) and combine multiple ggplot objects (e.g., to create a final figure)
library(ggpattern) # add patterns to graph (e.g., a density distribution)
library(moments) # to calculate skewdness and kurtosis of a density function

## Interaction coefficients
mdr_smap_coefficients <- fread("model_out/MDR_smap_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv") # read in the interaction coefficients

## Taxa information
taxa_info <-  read.csv("data/taxa_information.csv")

## Time series length
ts_length <- length(unique(mdr_smap_coefficients$time)) # 453

## Average interactions over time (with 0s removed)
# Calculate the mean strength of each link over the 6-year period and its 
# temporal 'connectivity' (how often that link is active/non-zero).
average_interactions <- mdr_smap_coefficients %>%
  group_by(target_1,target_2,group_1,group_2) %>%
  filter(!MDR_smap_coefficient==0) %>% # remove all 0 interactions
  mutate(count=1) %>% # dummy variable to count the number of links over time
  summarise(MDR_smap_coefficient=mean(MDR_smap_coefficient,na.rm=T), # average pairwise interactions over time
            connectivity=sum(count)/ts_length*100, .groups = "drop") %>% # count how often (%) links appear
  filter(!target_1==target_2) # remove self-interactions

## Calculate a threshold for interactions (only keep the strongest)
# We define 'strong' and 'frequent' interactions using the 50th percentile 
# as a filter to ensure the resulting matrix highlights the most interacting microbes.
inter_threshold <- quantile(abs(average_interactions$MDR_smap_coefficient),na.rm=T,0.50) # threshold for interaction strength, defined by the upper 50% quantile of the distribution

connec_threshold <- quantile(abs(average_interactions$connectivity),na.rm=T,0.50) # threshold for connectivity over time, defined by the upper 50% quantile of the distribution (only keep interactions that are realized during many time points)


# -----------------------------------------------------------------------
# 1.1 Fig. 1A: Interaction matrix
# -----------------------------------------------------------------------

inter_matrix_plot <- average_interactions %>%
  filter(abs(MDR_smap_coefficient)>=inter_threshold&connectivity>=connec_threshold)  # only keep the most relevant (strongest, frequent) interactions in which we have the most confidence

## Save strongest interactions in model_out folder
# write.csv(inter_matrix_plot,"model_out/strongest_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv",row.names=F)

## Effect on other taxa (rows), averaged over time
# Summarize how much each taxon influences the rest of the community.
effect_on_others_rows <- group_by(inter_matrix_plot,target_2,group_2) %>%
  mutate(count=1) %>% # dummy variable to count the number of links per target
  summarise(total_strength=sum(MDR_smap_coefficient,na.rm=T), # calculate sums and average of different interaction metrics per target
            total_strength_abs=sum(abs(MDR_smap_coefficient),na.rm=T),
            MDR_smap_coefficient=mean(MDR_smap_coefficient,na.rm=T),
            connectivity=mean(connectivity),links=sum(count), .groups = "drop") %>%
  # "Interactiveness" combines total net strength and total absolute strength to capture overall interactiveness.
  mutate(interactiveness=sqrt(total_strength_abs^2+total_strength^2)) %>% # calculate interactiveness, which contains information on the total strength and total absolute strength
  left_join(select(taxa_info,sequence_ID,method,domain),by=c("target_2"="sequence_ID")) %>% # add information on the target
  arrange(desc(method),desc(domain),desc(group_2),total_strength) %>% # arrange the data for plotting
  ungroup() %>% # ungroup data to prevent potential issues
  mutate(target_2_ID=seq(1,nrow(.),1)) # for plotting

factor_target_2_ID <- select(effect_on_others_rows,target_2,target_2_ID) # controls the order of target_2 labels

## Horizontal lines to separate higher taxonomic groups
# Logic: These lines visually group the matrix by Domain or Class.
h_lines <- effect_on_others_rows %>%
  ungroup() %>%
  mutate(row=seq(1,nrow(.),1)) %>%
  group_by(group_2) %>%
  summarise(line_index=max(row)+0.5)

## Y-axis labels for higher taxonomic groups
y_axis_breaks_labels <- effect_on_others_rows  %>%
  ungroup() %>%
  mutate(row=seq(1,nrow(.),1)) %>%
  group_by(group_2) %>%
  summarise(upper_limit=max(row),
            lower_limit=min(row)) %>%
  mutate(pos=(upper_limit-lower_limit)/2+lower_limit)


## Effect of other taxa (columns)
effect_of_others_columns <- group_by(inter_matrix_plot,target_1,group_1) %>%
  mutate(count=1) %>%
  summarise(total_strength=sum(MDR_smap_coefficient,na.rm=T),total_strength_abs=sum(abs(MDR_smap_coefficient),na.rm=T),MDR_smap_coefficient=mean(MDR_smap_coefficient,na.rm=T),connectivity=mean(connectivity),links=sum(count), .groups = "drop") %>%
  left_join(select(taxa_info,sequence_ID,method,domain),by=c("target_1"="sequence_ID")) %>% # add information on the target
  arrange(desc(method),desc(domain),desc(group_1),total_strength) %>% # arrange for plotting
  ungroup() %>%
  mutate(target_1_ID=seq(1,nrow(.),1))

factor_target_1_ID <- select(effect_of_others_columns,target_1,target_1_ID)

## Vertical lines to separate higher taxonomic groups
v_lines <- effect_of_others_columns %>%
  ungroup() %>%
  mutate(row=seq(1,nrow(.),1)) %>%
  group_by(group_1) %>%
  summarise(line_index=max(row)+0.5,
            min_row=min(row),
            max_row=max(row)) %>%
  mutate(row_diff=max_row-min_row+1,
         x_label_placement=max_row-row_diff/2) %>%
  arrange(x_label_placement)

## X-axis labels for higher taxonomic groups
x_axis_breaks_labels <- effect_of_others_columns  %>%
  ungroup() %>%
  mutate(row=seq(1,nrow(.),1)) %>%
  group_by(group_1) %>%
  summarise(upper_limit=max(row),
            lower_limit=min(row)) %>%
  mutate(pos=(upper_limit-lower_limit)/2+lower_limit)


# Color for interaction strength and direction
color_vector <- c('#003c30','#01665e','#35978f','#80cdc1','#c7eae5','#f5f5f5','#f6e8c3','#dfc27d','#bf812d','#8c510a','#543005')


# Define number of taxa to select as keystone interactions (10% most interactive)
# Logic: Keystone taxa are defined as the top ten percent (10%) of the interactiveness distribution.
threshold <- round((nrow(effect_on_others_rows)*0.1),0) # threshold for core microbial interactions (top 10%)

keystone_taxa <- taxa_info %>%
  arrange(-interactiveness) %>%
  slice(1:threshold) %>%
  left_join(select(effect_on_others_rows,target_2,target_2_ID),by=c("sequence_ID"="target_2"))

## Save keystone taxa in the folder model_out
# write.csv(select(keystone_taxa,-target_2_ID,"model_out/keystone_taxa_Scripps_Pier_ASVs_tp_2_28092024.csv",row.names=F)

keystone_taxa_vector_names <- pull(keystone_taxa,"sequence_ID")
keystone_taxa_vector_IDs <- pull(keystone_taxa,"target_2_ID")

## Create shaded rectangles for keystone taxa
# Visually highlight the "rows" belonging to keystone taxa in the matrix.
keystone_rectangles <- effect_on_others_rows %>%
  ungroup() %>%
  mutate(row=seq(1,nrow(.),1)) %>%
  right_join(keystone_taxa) %>%
  mutate(fill_color=ifelse(total_strength>0,color_vector[7],color_vector[5]))

## Assemble base data for plotting
inter_matrix_plot_data <- left_join(inter_matrix_plot,factor_target_1_ID) %>%
  left_join(factor_target_2_ID)

# Basis plot (empty)
plot <-  ggplot(inter_matrix_plot_data,aes(y=target_2_ID,x=target_1_ID)) +
  geom_tile(color="white",fill="white")

# Add background rectangles for keystone microbial interactions
# This loop draws transparent highlight bars across the matrix for keystone taxa.
for (i in 1:nrow(keystone_taxa)) {
  rect_index_i <-  keystone_rectangles$row[i]
  fill_color_i <- keystone_rectangles$fill_color[i]
  plot <- plot +
    annotate(geom="rect",xmin=0,xmax=max(inter_matrix_plot_data$target_1_ID)+1,ymin=(rect_index_i-0.5),ymax=(rect_index_i+0.5),fill=fill_color_i,color=fill_color_i,alpha=.4)
}

# Add pairwise interactions
# Heatmap where tile color represents the sign and strength of the interaction.
plot <- plot + geom_tile(color="white",fill="white") +
  geom_tile(data= inter_matrix_plot_data,color="black",aes(fill=MDR_smap_coefficient)) +
  scale_fill_gradientn(colours=color_vector,values = scales::rescale(c(-0.45, -0.01, 0, 0.01, 0.45)),limits=c(-0.45,0.45),name="Interaction\nstrength\n ",breaks=seq(-0.4,0.4,0.2)) + # change colors
  scale_y_continuous(limits=c(0,max(inter_matrix_plot_data$target_2_ID)+1),breaks=y_axis_breaks_labels$pos,labels=y_axis_breaks_labels$group_2, expand = c(0, 0)) + # add higher taxonomic groups as y-labels
  scale_x_continuous(limits=c(0,max(inter_matrix_plot_data$target_1_ID)+1),breaks=x_axis_breaks_labels$pos,labels=x_axis_breaks_labels$group_1, expand = c(0, 0)) # add higher taxonomic groups as x-labels

# Add horizontal lines to separate higher taxonomic groups
for (i in 1:nrow(h_lines)) {
  line_index_i <-  as.numeric(h_lines$line_index[i])
  plot <- plot +
    geom_hline(yintercept=line_index_i,size=0.4,color="black") 
}

# Add vertical lines to separate higher taxonomic groups
for (i in 1:nrow(v_lines)) {
  line_index_i <-  as.numeric(v_lines$line_index[i])
  plot <- plot +
    geom_vline(xintercept=line_index_i,size=0.5,color="black") 
}

# Assemble and format the final plot
plot_inter_matrix <- plot +
  geom_hline(yintercept=56.5,size=1,color="white") + # add thicker lines for domains
  geom_vline(xintercept=51.5,size=1,color="white") + # add thicker lines for domains
  geom_hline(yintercept=56.5,size=0.8,color="black") + # add thicker lines for domains
  geom_vline(xintercept=51.5,size=0.8,color="black") + # add thicker lines for domains
  theme_classic() +
  ylab("Effect on other taxa (rows)") +
  xlab("Effect of other taxa (columns)")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5),
        axis.text.y = element_text(size=5),
        text = element_text(size=12),
        legend.position="none",
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines'))

inter_matrix_legend <- get_legend(plot)    # extract the legend (this will help us to better position it in the final plots)


# -----------------------------------------------------------------------
# 1.2 Fig. 2A: Net interaction strength
# -----------------------------------------------------------------------

## Assemble the data for plotting
# Summarize the total net interaction strength
# for each taxon, highlighting keystone taxa with a unique color.
sum_strength_plot_data <- ungroup(inter_matrix_plot_data) %>% 
  select(target_2_ID,MDR_smap_coefficient) %>%
  group_by(target_2_ID) %>%
  summarise(MDR_smap_coefficient=sum(MDR_smap_coefficient), .groups = "drop") %>%
  mutate(direction=ifelse(MDR_smap_coefficient>0,"positive","negative")) %>%
  mutate(direction=ifelse(direction=="positive"&target_2_ID%in%keystone_taxa_vector_IDs,"positive_core",
                          ifelse(direction=="negative"&target_2_ID%in%keystone_taxa_vector_IDs,"negative_core",direction)))

## Plot the sum of interaction strength (for the rows)
# Bar chart aligned with the rows of the matrix showing the interaction strength.
plot_sum_strength <- ggplot(sum_strength_plot_data,aes(y=as.factor(target_2_ID),x=MDR_smap_coefficient)) +
  geom_bar(stat="identity",alpha=.8,aes(fill=direction,color=direction)) +
  scale_fill_manual(values=c(color_vector[5],color_vector[1],color_vector[7],color_vector[11])) +
  scale_color_manual(values=c(color_vector[1],color_vector[1],color_vector[11],color_vector[11])) +
  theme_classic() +
  scale_y_discrete(breaks = NULL) +
  xlab("") +
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines')) 

## Add the sum of interaction strength and the legend to the interaction matrix plot
plots_inter_matrix <- ggarrange(plot_inter_matrix,plot_sum_strength,labels=c("A"),ncol=2,nrow=1,widths=c(2.1,0.4),align="h")

## First part of the final figure (Fig. 1A)
column_1 <- ggarrange(plots_inter_matrix,inter_matrix_legend,ncol=2,nrow=1,widths=c(2.5,0.2))


# -----------------------------------------------------------------------
# 2. Fig. 2B: Net interaction strength and pairwise interaction strength
# -----------------------------------------------------------------------

## We first estimate the density distributions for the net interaction strength and all pairwise interaction strengths and then add them to the same plot

## Calculate the median of the net interaction strengths (from the rows)
median_distribution <- round(median(effect_on_others_rows$total_strength,na.rm=T),3)

## We want to color negative and positive links differently, thus we have to first extract and save the density curves
p <- ggplot(data=effect_on_others_rows,aes(x=total_strength))  +
  geom_density(aes(y = ..scaled..)) # keep density between 1 and 0
g <- ggplot_build(p)

## Extract density information
plot_data <- data.frame(y=g$data[[1]]$scaled,x=g$data[[1]]$x)

## Calculate the median of all pairwise interaction strengths
median_distribution2 <- round(median(inter_matrix_plot$MDR_smap_coefficient,na.rm=T),3)

## Again, we want to color negative and positive links differently, thus we have to first extract and save the density curves
p2 <- ggplot(inter_matrix_plot,aes(x=MDR_smap_coefficient)) +
  geom_density(aes(y = ..scaled..)) # keep density between 1 and 0
g2 <- ggplot_build(p2)

## Extract density information for pairwise interactions
plot_data_2 <- data.frame(y=g2$data[[1]]$scaled,x=g2$data[[1]]$x)

## Create fake data for the pattern legend
legend_data <- data.frame(pattern=c("pattern_1","pattern_2"),color=c("color_1","color_2"),value=c(1,2),variable=c("A","B"))

## Create fake plot for the pattern legend
legend_plot <- ggplot(legend_data,aes(y=value,x=variable)) +
  geom_bar_pattern(stat="identity",aes(pattern=pattern),fill = "white", color = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.01) +
  scale_pattern_manual(values=c("none","stripe"),labels=c(bquote(paste(sum(),"interaction strength (M=",.(median_distribution),")")),bquote(paste("Interaction strength (M=",.(median_distribution2),")"))),name="") +
  theme(legend.text=element_text(size=8),legend.key.size = unit(.6, 'cm'))

## Save the legend
legend <- get_legend(legend_plot) # extract the legend

## Assemble the final plot
# Combine the solid density (net strength) and the striped density (pairwise strength) 
# into one visualization.
plot_density <- ggplot() + 
  geom_ribbon(data=filter(plot_data,x>0),aes(x=x,y=y,ymin=0,ymax=y),color=color_vector[11],fill=color_vector[7],alpha=.8) +
  geom_ribbon(data=filter(plot_data,x<0),aes(x=x,y=y,ymin=0,ymax=y),color=color_vector[1],fill=color_vector[5],alpha=.8) +
  geom_ribbon_pattern(data=filter(plot_data_2,x>0),aes(x=x,y=y,ymin=0,ymax=y),alpha=0,color=color_vector[11],size=.8,
                      pattern_fill = color_vector[10],
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.03) +
  geom_ribbon_pattern(data=filter(plot_data_2,x<0),aes(x=x,y=y,ymin=0,ymax=y),alpha=0,color=color_vector[1],fill=color_vector[1],size=.8,
                      pattern_fill = color_vector[2],
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.03) +
  geom_vline(xintercept=0,color="white") +
  geom_vline(xintercept=0,linetype="dashed") +
  theme_classic() +
  scale_y_continuous(name="Density (scaled)",limits = c(0,1),expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0)) +
  xlab("Interaction strength") +
  annotation_custom(legend,
                    xmin = 0, ymin = 0.75) + # ad custom legend
  ## Format final plot
  theme(plot.margin=unit(c(.5,1,.5,.5),"cm"),text = element_text(size=12),plot.title = element_text(hjust = 0.5))


# -----------------------------------------------------------------------
# 3. Fig. 2C: Relative abundance and interactiveness
# -----------------------------------------------------------------------


## Get information on average relative abundance for all taxa
abundance_data <- taxa_info %>%
  select(sequence_ID,rel_ab_mean) %>%
  rename(target_2=sequence_ID) 

## Get information on net interaction strength and interactiveness for all taxa
total_strength_data <- effect_on_others_rows %>%
  select(target_2,interactiveness,total_strength) 

## Assemble the data for plotting
# Explore the relationship between a taxon's abundance
# and its interactiveness. Keystone taxa are highlighted in colored bars.
rank_plot <- left_join(total_strength_data,abundance_data, by = "target_2") %>%
  mutate(direction=ifelse(total_strength>0,"positive","negative")) %>%
  mutate(keystone_taxa=ifelse(target_2%in%unique(keystone_taxa$sequence_ID),"yes","no")) %>%
  mutate(fill_keystone_taxa=ifelse(keystone_taxa=="yes"&total_strength>0,color_vector[8],
                                   ifelse(keystone_taxa=="yes"&total_strength<0,color_vector[4],NA))) %>%
  mutate(color_keystone_taxa=ifelse(keystone_taxa=="yes"&total_strength>0,color_vector[11],
                                    ifelse(keystone_taxa=="yes"&total_strength<0,color_vector[1],NA))) %>%
  arrange(-interactiveness) %>%
  mutate(rank=seq(1,nrow(.),1))

## Produce and save the final plot
plot_rank_strength_abundance <- ggplot(rank_plot,aes(y=interactiveness,x=rank)) +
  geom_vline(xintercept = round(nrow(rank_plot)*0.1,0)) + # Vertical line at the keystone decile
  geom_bar(aes(y=rel_ab_mean/5,fill=direction,color=direction),stat="identity",size=0.15) +
  geom_bar(data=filter(rank_plot,rank%in%seq(1,threshold,1)),aes(y=rel_ab_mean/5),stat="identity",size=0.15,fill="white",color="white") +
  geom_bar(data=filter(rank_plot,rank%in%seq(1,threshold,1)),aes(y=rel_ab_mean/5),stat="identity",size=0.15,color=rank_plot$color_keystone_taxa[1:threshold],fill=rank_plot$fill_keystone_taxa[1:threshold]) +
  geom_point(aes(fill=direction,color=direction),shape=21) +
  geom_point(data=filter(rank_plot,rank%in%seq(1,threshold,1)),shape=21,fill="white",color="white") +
  geom_point(data=filter(rank_plot,rank%in%seq(1,threshold,1)),shape=21,color=rank_plot$color_keystone_taxa[1:threshold],fill=rank_plot$fill_keystone_taxa[1:threshold]) +
  ylab("Interactiveness") +
  scale_y_continuous(expand = c(0, 0),sec.axis = sec_axis(~.*5, name="Relative abundance (%)")) +
  theme_classic() +
  xlab("Rank") +
  scale_fill_manual(values=c(color_vector[5],color_vector[7])) +
  scale_color_manual(values=c(color_vector[1],color_vector[11])) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none") +
  # format final plot
  theme(plot.margin=unit(c(.5,1,.5,.5),"cm"),text = element_text(size=12),legend.position = "none")


## Combine Figure 2 B and C into one column
column_2 <- ggarrange(plot_density,plot_rank_strength_abundance,labels=c("B","C"),ncol=1,nrow=2,align="hv")


# Create final figure and save it in plots
ggarrange(column_1,column_2,widths=c(1.8,1.2),ncol=2,nrow=1,align="hv")

ggsave("plots/Fig_2_microbial_interactions.pdf",height=7.3,width=13)


# -----------------------------------------------------------------------
# 4. Summary and test statistics
# -----------------------------------------------------------------------

## Summary
summary_stats <- effect_on_others_rows %>%
  mutate(direction=ifelse(total_strength>0,"positive","negative")) %>%
  mutate(count=1) %>%
  group_by(direction) %>%
  summarise(count=sum(count),
            strength=mean(total_strength))


## Calculate skewdness and kurtosis for Figure 2B (tutorial: https://www.statology.org/skewness-kurtosis-in-r/)
# Skewness and Kurtosis help characterize the 
# asymmetry and the "heavy-tailedness" of the interaction strength distribution.

## Total interaction strength
skewness(effect_on_others_rows$total_strength) # calculate skewness, here 0.402
kurtosis(effect_on_others_rows$total_strength) # calculate kurtosis, here 4.32

# Jarque-Bera test for normality
jarque.test(effect_on_others_rows$total_strength)

# Pairwise interactions
skewness(inter_matrix_plot$MDR_smap_coefficient) # calculate skewness, here 0.82
kurtosis(inter_matrix_plot$MDR_smap_coefficient) # calculate kurtosis, here 9.63

jarque.test(inter_matrix_plot$MDR_smap_coefficient)
