###########################################################
# Analyzing ASVs
# 
# Author: JS Huisman
###########################################################

library(tidyverse)
library(ggplot2)
library(phyloseq)
library(stringdist) # for the one base away analysis
library(patchwork)


proc_data_dir = '../data/1_community_experiment'
fig_data = '../figures/1_community_experiment'


theme_set(theme_minimal() + theme(text = element_text(size = 20)))

color_list <- c('Estuary' = '#A72626', 
                'Marine 2' = '#F5BD23', 
                'Marine 1' = '#356FEA', 
                'Brackish' = '#004D07')

color_salt2 <- c('16%' = '#7EFAD5',
                 '31%' = '#1E88E5',
                 '46%' = '#6B02B5')

## LOAD DATA: ps_noRare was preprocessed to ASVs that appear only once across all samples #############################

rel_abund = read_csv(paste0(proc_data_dir, '/rel_abund_df_norare.csv'))

rel_abund <- rel_abund %>%
  mutate(inoc = factor(inoc, levels = c('charles', 'ica', 'nahant', 'fucus')),
         inoc_map = case_when(
           inoc == 'charles' ~ 'Brackish',
           inoc == 'ica' ~ 'Estuary',
           inoc == 'nahant' ~ 'Marine 1',
           inoc == 'fucus' ~ 'Marine 2'
         ),
         family = factor(Family), 
         salt_num = as.numeric(gsub('%', '', salt)) )

rel_abund_1perc = rel_abund %>%
  filter(rel_abund >= 0.01) %>%
  mutate(family = factor(Family))

inoc_label = c('charles' = 'Brackish', 'ica' = 'Estuary',
           'nahant' = 'Marine 1', 'fucus' = 'Marine 2')
salt_label = c('16%' = '16 (g/L)', '31%' = '31 (g/L)',
               '46%' = '46 (g/L)')

mcn_df = read_csv(paste0(proc_data_dir, '/mcn_df_norare.csv')) %>%
  mutate(inoc = factor(inoc, levels = c('charles', 'ica', 'nahant', 'fucus')),
         inoc_map = case_when(
           inoc == 'charles' ~ 'Brackish',
           inoc == 'ica' ~ 'Estuary',
           inoc == 'nahant' ~ 'Marine 1',
           inoc == 'fucus' ~ 'Marine 2'
         ) )

ps_norare = readRDS(paste0(proc_data_dir, '/ps_norare.rds'))

diversity_df <- read_csv(paste0(proc_data_dir, '/diversity_df_norare.csv')) %>%
  mutate(salt = gsub('%', '', salt),
         inoc_map = case_when(
           inoc == 'charles' ~ 'Brackish',
           inoc == 'ica' ~ 'Estuary',
           inoc == 'nahant' ~ 'Marine 1',
           inoc == 'fucus' ~ 'Marine 2'
         ) )

## ANALYZE DATA ##########################################################################################################

############################################################################################
## ASV abundance on day 6 - colored by family - Fig. S1 ####

abund_15C = ggplot(rel_abund %>% filter(rel_abund >= 0.01, day == 'D6', temp == '15C')) +
  geom_bar(aes(x = rep, y = rel_abund, fill = family), color = 'black', stat = 'identity', show.legend = F) +
  facet_grid(rows = vars(inoc), cols = vars(salt), labeller = labeller(inoc = inoc_label, salt = salt_label)) +
  labs(x = 'Replicate', y = 'Relative abundance', fill = 'Family', title = '15C') +
  scale_fill_viridis_d(option = 'viridis')

abund_20C = ggplot(rel_abund %>% filter(rel_abund >= 0.01, day == 'D6', temp == 'roomT')) +
  geom_bar(aes(x = rep, y = rel_abund, fill = family), color = 'black', stat = 'identity', show.legend = T) +
  facet_grid(rows = vars(inoc), cols = vars(salt), labeller = labeller(inoc = inoc_label, salt = salt_label)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = 'Replicate', y = 'Relative abundance', fill = 'Family', title = '20C') +
  scale_fill_viridis_d(option = 'viridis')

abund_15C + abund_20C +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom',
        panel.spacing.y = unit(1, 'lines'))

ggsave(paste0(fig_data, '/FigS1_ASV_abund_1perc.pdf'), width = 14, height = 8)


## Panel 1A, example ASV abundance from roomT, rep R1
ggplot(rel_abund %>% filter(day == 'D6', temp == 'roomT', rep == 'R1')) +
  geom_bar(aes(x = factor(salt_num), y = rel_abund, fill = OTU), stat = 'identity', show.legend = F) +
  facet_wrap(nrow = 1, vars(inoc), labeller = labeller(inoc = inoc_label)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = 'Salinity (g/L)', y = 'Relative abundance', fill = 'ASV') +
  scale_fill_viridis_d(option = 'viridis')

ggsave(paste0(fig_data, '/Fig1A_abundance.pdf'), width = 8, height = 4)

##############################################
## Dynamics of ASVs - Fig. S2 ####

## Per inoculum and temperature, across the 7 dilution cycles
plot_list = list()
i = 1
for (inoc_cond in c('charles', 'ica', 'nahant', 'fucus')){
  for(temp_cond in c('15C', 'roomT')){
    temp_name = ifelse(temp_cond == 'roomT', '20C', '15C')
    
    new_plot <- ggplot(rel_abund_1perc %>% filter(temp == temp_cond, inoc == inoc_cond) %>% 
             mutate(day_cont = as.numeric(gsub('D', '', day)))) +
      geom_bar(aes(x = day_cont, y = rel_abund, fill = family), stat = 'identity', show.legend = T) +
      facet_grid(vars(rep), vars(salt), labeller = labeller(salt = salt_label)) +
      labs(x = 'Dilution cycle', y = 'Rel. abundance', fill = 'Family', title = paste0(inoc_label[inoc_cond], ', ', temp_name)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) + 
      scale_x_continuous(breaks = c(1, 3, 5, 6, 7)) + 
      scale_fill_viridis_d(option = 'viridis', drop = F) +
      theme(legend.position = 'none',
            plot.title = element_text(size = 20),
            panel.spacing.x = unit(2, 'lines'))
    
    plot_list[[i]] = new_plot
    i = i+1
    #ggsave(paste0(fig_data, '/abund_asvs_time_', inoc_cond, '_', temp_name, '.pdf'), plot = new_plot, width = 10, height = 5)
  }
}

# Combine all plots using patchwork - create Fig. S2
(plot_list[[1]] + theme(axis.title.x = element_blank()) ) + (plot_list[[2]] + theme(axis.title.x = element_blank())) +
(plot_list[[3]] + theme(strip.text.x = element_blank(), axis.title.x = element_blank())) + (plot_list[[4]] + theme(strip.text.x = element_blank(), axis.title.x = element_blank())) + 
(plot_list[[5]] + theme(strip.text.x = element_blank(), axis.title.x = element_blank())) + (plot_list[[6]] + theme(strip.text.x = element_blank(), axis.title.x = element_blank())) +
(plot_list[[7]] + theme(strip.text.x = element_blank())) + (plot_list[[8]] + theme(strip.text.x = element_blank())) +
  plot_layout(ncol = 2, guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom',
       plot.title = element_blank(),
       panel.spacing.x = unit(2, 'lines'))

ggsave(paste0(fig_data, '/FigS2_abund_asvs_time_all.pdf'), width = 14, height = 18)


############################################################################################
## Diversity

# final diversity across conditions - at both temperatures
ob_div_plot <- ggplot(diversity_df %>% filter(salt != 'inoculum', day == 'D6', diversity_metric == 'Observed')) +
  geom_smooth(aes(x = as.numeric(salt), y = div_value, color = inoc_map, group = interaction(inoc)), method = 'lm', se = F, linewidth = 2) +
  geom_point(aes(x = as.numeric(salt), y = div_value, color = inoc_map, shape = temp), alpha = 1, size = 4) +
  #facet_wrap(vars(temp), labeller = labeller(temp = c('15C' = '15C', 'roomT' = '20C'))) +
  scale_color_manual(values = color_list) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  coord_cartesian(ylim = c(0, 210)) +
  labs(x = 'Salinity (g/L)', y = 'Observed diversity', color = 'Inoculum', shape = 'Temperature') +
  theme(text = element_text(size = 20),
        legend.position = 'none')

ob_div_plot
ggsave(paste0(fig_data, '/Fig1D_richness_d6.pdf'), height = 4, width = 4)
