# -----------------------------------------------------------------------
# Title: Microbial interactions and water temperature (Fig. 4)
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: This script is used to recreate Fig. 4 showing how microbial interactions changes as a function of water temperature.
# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

## Required Libraries

library(tidyverse) # data formatting
library(data.table) # read and write large data files
library(ggpubr) # arrange multiple plots in the same panel
library(ggridges) # Plot multiple density distributions on the same plot
library(ggtext) # Add text to plots

###################################

# -----------------------------------------------------------------------
# 1. Load required data
# -----------------------------------------------------------------------

## Data (for water temperature)
data <- read.csv("data/data_sequences_0.1_rel_ab_0.5_occ_binned_4_days_with_temperature.csv") %>%
  mutate(time=date(date)) %>%
  select(time,temperature)

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

# 'label_edited' uses Markdown syntax (*Italics* and **Bold**) which 
# will be rendered by the ggtext package in the plot themes.
keystone_taxa <- keystone_taxa %>%
  left_join(group_short_IDs) %>%
  mutate(label_edited=paste0("*",name_edited,"* (**",group_short_ID,"**)"))

# -----------------------------------------------------------------------
# 2. Fig. 4A-B: Temperature dependency of interactiveness and percent 
# facilitation of selected keystone taxa
# -----------------------------------------------------------------------

time_var_inter_matrix <- interactions_time %>%
  filter(!MDR_smap_coefficient==0) %>% # remove all 0 interactions
  mutate(count=1) %>%
  mutate(facilitation=ifelse(MDR_smap_coefficient>0,MDR_smap_coefficient,NA), # record if a link is facilitative (coefficient > 0) or inhibiting (coefficient < 0)
         inhibition=ifelse(MDR_smap_coefficient<0,MDR_smap_coefficient,NA)) %>%
  # Aggregate data per target and time point
  group_by(time,target_2,group_2) %>%
  summarise(total_strength=sum(MDR_smap_coefficient,na.rm=T), 
            total_strength_abs=sum(abs(MDR_smap_coefficient),na.rm=T),
            MVD_smap_coefficient=mean(MDR_smap_coefficient,na.rm=T),
            links=sum(count),
            facilitation=sum(facilitation,na.rm = T),
            inhibition=sum(inhibition,na.rm=T)) %>%
  # Add information on water temperature
  mutate(time=date(time)) %>%
  left_join(data) %>%
  mutate(temperature=round(temperature,0)) %>% # round water temperature values (to create discrete bins)
  # Interactiveness is a summary metric of net and absolute interaction strengths.
  mutate(interactiveness=sqrt(total_strength^2+total_strength_abs^2)) %>% # calculate interactiveness
  # Percent facilitation measures the relative proportion of positive vs negative influence.
  mutate(perc_facilitation=facilitation/(abs(inhibition)+facilitation)*100) %>%
  ungroup() %>%
  select(temperature,target_2,group_2,interactiveness,perc_facilitation) %>%
  # Aggregate data per target and temperature bin
  group_by(temperature,target_2,group_2) %>%
  summarise(mean_interactiveness=mean(interactiveness),
            sd_interactiveness=sd(interactiveness),
            n_obs=length(interactiveness),
            mean_perc_facilitation=mean(perc_facilitation),
            sd_perc_facilitation=sd(perc_facilitation))


## Gather plot data
keystone_taxa_var_inter_matrix <- time_var_inter_matrix %>%
  filter(target_2%in%pull(keystone_taxa,"sequence_ID")) %>%
  mutate(se_interactiveness=sd_interactiveness/sqrt(n_obs), # calculate standard error of the mean (SE)
         se_perc_facilitation=sd_perc_facilitation/sqrt(n_obs)) %>%
  left_join(select(keystone_taxa,sequence_ID,label_edited),by=c("target_2"="sequence_ID")) %>%
  mutate(label_edited=factor(label_edited,levels=pull(keystone_taxa,"label_edited")))

## Interactiveness

## Plot for all keystones
ggplot(keystone_taxa_var_inter_matrix,aes(y=mean_interactiveness,x=temperature)) +
  geom_errorbar(aes(ymin=mean_interactiveness-se_interactiveness,ymax=mean_interactiveness+se_interactiveness,color=label_edited),width=0.2) +
  geom_line(aes(color=label_edited)) +
  geom_point(fill="white",color="white",shape=21,size=3) +
  geom_point(aes(fill=label_edited),color="black",shape=21,size=3,alpha=.8) +
  theme_classic() +
  facet_wrap(~label_edited,scale="free") +
  xlab("Temperature (°C)") +
  ylab("Interactiveness") +
  scale_color_manual(values=pull(keystone_taxa,"color")) +
  scale_fill_manual(values=pull(keystone_taxa,"color")) +
  scale_x_continuous(breaks=seq(13,26,2)) +
  theme(legend.text = element_markdown(),
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines'),
        text = element_text(size=12),
        strip.text.x = element_markdown(size=8),
        strip.background.x=element_rect(color = NA),
        legend.position = "none")

## Selection of keystones for Fig. 4
plot_keystone_interactiveness_subset <- ggplot(filter(keystone_taxa_var_inter_matrix,target_2%in%c("Proteobacteria_22_1","Proteobacteria_7166_1","Dinoflagellata_1198_1","Dinoflagellata_1793_1") ),aes(y=mean_interactiveness,x=temperature)) +
  geom_errorbar(aes(ymin=mean_interactiveness-se_interactiveness,ymax=mean_interactiveness+se_interactiveness,color=label_edited),width=0.2) +
  geom_line(aes(color=label_edited)) +
  geom_point(fill="white",color="white",shape=21,size=3) +
  geom_point(aes(fill=label_edited),color="black",shape=21,size=3,alpha=.8) +
  theme_classic() +
  facet_wrap(~label_edited,scale="free",nrow=1) +
  xlab("Temperature (°C)") +
  ylab("Interactiveness") +
  scale_color_manual(values=pull(filter(keystone_taxa,sequence_ID%in%c("Proteobacteria_22_1","Proteobacteria_7166_1","Dinoflagellata_1198_1","Dinoflagellata_1793_1")),"color")) +
  scale_fill_manual(values=pull(filter(keystone_taxa,sequence_ID%in%c("Proteobacteria_22_1","Proteobacteria_7166_1","Dinoflagellata_1198_1","Dinoflagellata_1793_1")),"color")) +
  scale_x_continuous(breaks=seq(13,26,2)) +
  theme(legend.text = element_markdown(),
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines'),
        text = element_text(size=12),
        strip.text.x = element_markdown(size=8),
        strip.background.x=element_rect(color = NA),
        legend.position = "none")

## Facilitation (%)

## Color vector
color_vector <- c('#003c30','#01665e','#35978f','#80cdc1','#c7eae5','#f5f5f5','#f6e8c3','#dfc27d','#bf812d','#8c510a','#543005')

## Plot for all keystones
ggplot(keystone_taxa_var_inter_matrix,aes(y=mean_perc_facilitation,x=temperature)) +
  annotate(geom="rect",xmin=-Inf,xmax=Inf,ymin=50,ymax=Inf,fill=color_vector[7],color=color_vector[7],alpha=.4) +
  annotate(geom="rect",xmin=-Inf,xmax=Inf,ymin=50,ymax=-Inf,fill=color_vector[5],color=color_vector[5],alpha=.4) +
  geom_hline(yintercept=50,linetype="dashed") +
  geom_errorbar(aes(ymin=mean_perc_facilitation-se_perc_facilitation,ymax=mean_perc_facilitation+se_perc_facilitation,color=label_edited),width=0.2) +
  geom_line(aes(color=label_edited)) +
  geom_point(fill="white",color="white",shape=21,size=3) +
  geom_point(aes(fill=label_edited),color="black",shape=21,size=3,alpha=.8) +
  theme_classic() +
  facet_wrap(~label_edited,scale="free") +
  xlab("Temperature (°C)") +
  ylab("Facilitation and Inhibition") +
  scale_color_manual(values=pull(keystone_taxa,"color")) +
  scale_fill_manual(values=pull(keystone_taxa,"color")) +
  scale_x_continuous(breaks=seq(13,26,2)) +
  theme(legend.text = element_markdown(),
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines'),
        text = element_text(size=12),
        strip.text.x = element_markdown(size=8),
        strip.background.x=element_rect(color = NA),
        legend.position = "none")


## Selection of keystones for Fig. 4
plot_keystone_facilitation_subset <- ggplot(filter(keystone_taxa_var_inter_matrix,target_2%in%c("Archaea_164_1","Proteobacteria_11431_1","Bacteroidetes_919_2","Proteobacteria_11429_2")),aes(y=mean_perc_facilitation,x=temperature)) +
  annotate(geom="rect",xmin=-Inf,xmax=Inf,ymin=50,ymax=Inf,fill=color_vector[7],color=color_vector[7],alpha=.4) +
  annotate(geom="rect",xmin=-Inf,xmax=Inf,ymin=50,ymax=-Inf,fill=color_vector[5],color=color_vector[5],alpha=.4) +
  geom_hline(yintercept=50,linetype="dashed") +
  geom_errorbar(aes(ymin=mean_perc_facilitation-se_perc_facilitation,ymax=mean_perc_facilitation+se_perc_facilitation,color=label_edited),width=0.2) +
  geom_line(aes(color=label_edited)) +
  geom_point(fill="white",color="white",shape=21,size=3) +
  geom_point(aes(fill=label_edited),color="black",shape=21,size=3,alpha=.8) +
  theme_classic() +
  facet_wrap(~label_edited,scale="free",nrow=1) +
  xlab("Temperature (°C)") +
  ylab("Facilitation (%)") +
  scale_color_manual(values=pull(filter(keystone_taxa,sequence_ID%in%c("Archaea_164_1","Proteobacteria_11431_1","Bacteroidetes_919_2","Proteobacteria_11429_2")),"color")) +
  scale_fill_manual(values=pull(filter(keystone_taxa,sequence_ID%in%c("Archaea_164_1","Proteobacteria_11431_1","Bacteroidetes_919_2","Proteobacteria_11429_2")),"color")) +
  scale_x_continuous(breaks=seq(13,26,2)) +
  theme(legend.text = element_markdown(),
        plot.margin = unit(c(0.2,0.1,0,0.1), 'lines'),
        text = element_text(size=12),
        strip.text.x = element_markdown(size=8),
        strip.background.x=element_rect(color = NA),
        legend.position = "none")

# -----------------------------------------------------------------------
# 3. Fig. 4C: Turnover in keystones along a gradient of water temperature
# -----------------------------------------------------------------------

## Get temperature levels
temperature_levels <- data.frame(temperature=unique(pull(time_var_inter_matrix,"temperature")))

## Select keystones per temperature level
# At each 1°C bin, we rank the taxa by interactiveness and select 
# the top 16 to visualize the succession of keystone taxa.
for (i in 1:nrow(temperature_levels)){
  
  time_var_inter_matrix_i <- filter(time_var_inter_matrix,temperature==temperature_levels$temperature[i]) %>%
    ungroup() %>%
    arrange(-mean_interactiveness) %>%
    slice(1:16) # 16 is our threshold for the number of keystones (=10% of all taxa)
  
  temp_keystone_taxa_i <- select(time_var_inter_matrix_i,target_2,temperature)
  
  if(i==1){temp_keystone_taxa<-temp_keystone_taxa_i}else{temp_keystone_taxa<-rbind(temp_keystone_taxa,temp_keystone_taxa_i)}
  
}

## Define order and assemble plot data
factor_levels <- taxa_info %>%
  filter(sequence_ID%in%pull(temp_keystone_taxa,"target_2")) %>%
  pull(sequence_ID) %>%
  rev()

groups_colors <- drop_na(unique(select(taxa_info,group,color)))

plot_data <- temp_keystone_taxa %>%
  left_join(taxa_info,by=c("target_2"="sequence_ID")) %>%
  left_join(groups_colors) %>%
  mutate(target_2=factor(target_2,levels=factor_levels)) 

axis_labels <- data.frame(sequence_ID=factor_levels) %>%
  left_join(group_short_IDs) %>%
  mutate(group_short_ID=ifelse(is.na(group_short_ID)," ",group_short_ID)) %>%
  pull("group_short_ID")


## Final plot (Fig. 4C)
# Heatmap-style plot showing when each keystone ASV "appears" as a 
# keystone taxa across the 13-26°C range.
plot_keystone_temperature <- ggplot(plot_data,aes(y=target_2,x=temperature)) +
  geom_hline(aes(yintercept=target_2),color="white") +
  geom_hline(data=filter(plot_data,target_2%in%pull(keystone_taxa,"sequence_ID")),aes(yintercept=target_2),size=3,color="gray") +
  geom_tile(fill="white") +
  geom_point(color=pull(plot_data,"color"),size=2) +
  theme_classic() +
  scale_x_continuous(breaks=seq(13,26,2)) +
  scale_y_discrete(labels=axis_labels) +
  theme(panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_line(size=1),
        legend.position="none",
        text = element_text(size=12)) +
  ylab("") +
  xlab("Temperature (°C)")


# -----------------------------------------------------------------------
# 4. Fig. 4D-E:  Density distributions as a function of temperature
# -----------------------------------------------------------------------

## Statistical test: Kolmogorov-Smirnov
# We perform pairwise KS tests to see if the shape of the 
# interactiveness distribution changes significantly between temperatures.
temp_combinations <- expand.grid(temp_1=unique(time_var_inter_matrix$temperature),temp_2=unique(time_var_inter_matrix$temperature)) %>%
  mutate(ks_test=NA, variable="Interactiveness")

for(i in 1:nrow(temp_combinations)){ 
  
  temp_1_i <- temp_combinations$temp_1[i]
  temp_2_i <- temp_combinations$temp_2[i]
  
  distribution_1 <- filter(time_var_inter_matrix,temperature==temp_1_i) %>% pull("mean_interactiveness")
  distribution_2 <- filter(time_var_inter_matrix,temperature==temp_2_i) %>% pull("mean_interactiveness")
  
  ks_test <- ks.test(distribution_1,distribution_2)
  temp_combinations$ks_test[i] <- ks_test[["p.value"]]
}

# Assign letters to show significant differences
temp_letters <- data.frame(temperature=unique(time_var_inter_matrix$temperature),temp_letter=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n"),test_letters=NA)

for (i in 1:nrow(temp_letters)){
  temp_i <- temp_letters$temperature[i]
  ks_test_temp_i <- filter(temp_combinations,temp_1==temp_i) %>%
    mutate(sig=ifelse(ks_test<0.05,1,0)) %>%
    mutate(temp_letter_2=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n")) %>%
    filter(sig==0) %>%
    pull("temp_letter_2") 
  temp_letters$test_letters[i] <- str_c(ks_test_temp_i , collapse = "")
}

## Color palette setup
color_vector <- c('#053061','#2166ac','#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f')
colors <- colorRampPalette(color_vector)

## Final plot (Fig. 4D: Interactiveness Ridges)
# geom_density_ridges_gradient visualizes the community-wide shift in
# interaction strength across the thermal gradient.
plot_interactiveness <- ggplot(time_var_inter_matrix, aes(x = mean_interactiveness, y = as.character(temperature))) +
  geom_density_ridges_gradient(aes(,fill=as.character(temperature)),scale = 3, rel_min_height = 0.01,quantile_lines=TRUE, quantile_fun=function(mean_interactiveness,...)median(mean_interactiveness),alpha=.8) +
  geom_label(data=temp_letters,aes(label=test_letters,y=as.character(temperature),x=1.25), angle = 90,label.size=NA) +
  xlim(0,1.5) +
  scale_fill_manual(values=colors(14)) +
  scale_y_discrete(breaks=seq(13,26,2)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=12),
        panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_line(size=1)) +
  ylab("Temperature (°C)") +
  xlab("Interactiveness") +
  coord_flip()


## Final plot (Fig. 4E: Facilitation Ridges)
plot_perc_facilitation <- ggplot(time_var_inter_matrix, aes(x = mean_perc_facilitation, y = as.character(temperature))) +
  geom_density_ridges_gradient(aes(,fill=as.character(temperature)),scale = 3, rel_min_height = 0.01,quantile_lines=TRUE, quantile_fun=function(perc_facilitation,...)median(perc_facilitation),alpha=.8) +
  geom_label(data=temp_letters,aes(label=test_letters,y=as.character(temperature),x=-7), angle = 90,label.size = NA) +
  geom_vline(xintercept=50,linetype="dashed") +
  scale_y_discrete(breaks=seq(13,26,2)) +
  scale_fill_manual(values=colors(14)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=12),
        panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_line(size=1)) +
  scale_x_continuous(breaks=seq(0,100,25)) +
  ylab("Temperature (°C)") +
  xlab("Facilitation (%)") +
  coord_flip()


## Assemble the final Figure 4
# Logic: Arranging all sub-plots into a single comprehensive panel for publication.
row_1 <- ggarrange(plot_keystone_interactiveness_subset,plot_keystone_facilitation_subset,nrow=2,ncol=1,labels=c("A","B"))
row_2 <- ggarrange(plot_keystone_temperature,plot_interactiveness,plot_perc_facilitation,labels=c("C","D","E"),ncol=3,nrow=1,widths=c(1.6,1.2,1.2))

ggarrange(row_1,row_2,nrow=2,ncol=1)

## Save the final figure
ggsave("plots/Fig_4_temperature_dependency.pdf",height=12,width=13)