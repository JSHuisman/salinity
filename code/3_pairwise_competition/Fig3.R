######################################################################################################################
## Analyzing competition data from Oct 2024
## Author: JS Huisman
######################################################################################################################

library(tidyverse)
#library(readxl)
library(lubridate)

proc_data_dir = '../data/3_pairwise_competition'
fig_dir = '../figures/3_pairwise_competition'

theme_set(theme_minimal() + theme(text = element_text(size = 20)))

######################################################################################################################

unique_pairs = c('40-2', '40-8', '40-39', '40-II.F2', 
                 '6-25', '2-25', '2-32', '2-II.F2')
faster_pair = c("6-25" = 6, 
                '2-25' = 2, '2-32' = 32, "2-II.F2" = 2, 
                "40-2" = 2, "40-8" = 8, '40-39' = 40, "40-II.F2" = 40)

competition_data <- read_csv(file = paste0(proc_data_dir, '/competition_plating_data.csv')) 
summarized_comp_data <- read_csv(file = paste0(proc_data_dir, '/summarized_data.csv'))
rep_summarized_comp_data <- read_csv(file = paste0(proc_data_dir, '/rep_summarized_data.csv'))

average_growth_df <- read_csv('../data/1.2_isolates/average_growthrates.csv')

######################################################################################################################

color_list_2 = c('#6E969B', '#004D07', '#82E2F3', '#F5BD23', '#741B7B',  '#C0E49B',  '#A72626',
                 '#B77406',  '#356FEA', '#000000')


salinity_colors <- c('16' = '#7EFAD5', '31' = '#1E88E5',
                     '46' = '#6B02B5', '61' = 'darkblue')

species_colors = setNames(color_list_2[1:8] , c('2', '6', '8', '25', '32', '39', '40', 'II.F2'))
species_colors_byname = setNames(color_list_2[1:8] , c('Pa', 'Sh', 'Po', 'Ar', 'Pn', 'Pf', 'Sx', 'Pp'))

pair_colors = setNames(c(color_list_2[2], color_list_2[1],
                         color_list_2[5], color_list_2[1], color_list_2[1],
                         color_list_2[3], color_list_2[7], color_list_2[7]), names(faster_pair))

species_abbreviations = c('2'= 'Pa', '40'='Sx', '25'='Ar', '8'='Po', '6'='Sh',
                          '32'='Pn', '39'='Pf', 'II.F2'='Pp')

pair_abbreviations = c('2-25' = 'Pa-Ar', '2-32'='Pa-Pn', '2-II.F2'='Pa-Pp', '40-2'='Sx-Pa', 
                       '40-39'='Sx-Pf', '40-8'='Sx-Po', '40-II.F2'='Sx-Pp', '6-25'='Sh-Ar')

##### Ratios after cycle 7  ###############################################

# # Ratios across all pairs after cycle 7 - not included in manuscript
# ggplot(competition_data  %>% filter(day_ID == 'D7', competition != 100) ) +
#     geom_point(aes(x = salt , y = ratio_fast, shape = competition), color = 'lightgrey' , size = 3 )+ 
#     geom_line(aes(x = salt , y = ratio_fast, linetype = competition, 
#                   group = interaction(bio_rep, competition)), color = 'lightgrey', show.legend = F)+ 
#     geom_line(data= summarized_comp_data %>% filter(day_ID == 'D7', pair_type != 'mono'), 
#               aes(y = mean_ratio_fast, x = salt, color = pair), linewidth = 1, show.legend = F ) +
#     facet_wrap(vars(pair), nrow = 2, labeller = labeller(pair = pair_abbreviations)) +
#     coord_cartesian(ylim = c(0, 1)) +
#     labs(x = 'Salinity (g/L)', y = 'Ratio faster grower', shape = 'Inoculation ratio') +
#   scale_color_manual(values = color_list_2) +
#   scale_x_continuous(breaks=c(16, 31, 46, 61)) +
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = 'bottom',
#         panel.spacing = unit(2, 'lines')) 
#   
# ggsave(paste0(fig_dir, '/Ratios_D7.png'), width = 10, height = 8)


## For each pair individually, plot the ratios of both strains - Fig 3B/D + S14
for(pair_id in unique_pairs){

    # reshape the dataframe for plotting; select only D7 and the specific pair
    # raw data
    long_pair_data <- competition_data %>% 
      filter(day_ID == 'D7', pair == pair_id) %>%
      pivot_longer(c(CFU_1, CFU_2), values_to = 'CFU', names_prefix = 'CFU_') %>%
      mutate(strain = ifelse(name == 1, strain1, strain2),
             ratio = CFU /(CFU_total))
    
    # Average calculated across all inoculation ratios, all bio replicates
    mean_long_pair_data <- long_pair_data %>% 
      group_by(salinity, strain) %>%
      summarize(mean_ratio = mean(ratio, na.rm = T),
                sd_ratio = sd(ratio, na.rm = T)) %>%
      mutate(strain_label = species_abbreviations[strain])
    
    ggplot(mean_long_pair_data) +
      geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio, 
                      x = salinity, color = strain_label), width = 2, linewidth = 1, show.legend = F) +
      #geom_point(data = long_pair_data, aes(x = salinity, y = ratio, color = strain), size = 3)+
      geom_point(aes(x = salinity, y = mean_ratio, color = strain_label), size = 3)+
      geom_line(aes(y = mean_ratio, x = salinity, color = strain_label, group = strain_label), linewidth = 2) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_continuous(breaks=c(16, 31, 46, 61)) +
      labs(x = 'Salinity (g/L)', y = 'Rel. abundance', color = "Isolate") +
      theme(axis.title = element_text(size = 25),
        legend.position = c(0.8, 0.5)) + 
      scale_color_manual(values = species_colors_byname) 
    
    ggsave(paste0(fig_dir, '/Ratios_', pair_id, '.pdf'), width = 8, height = 5) 

}

## Summary figure with change in faster grower relative to 30 g/L - Fig 3F
# reshape the dataframe for plotting; select only D7; calculate average per pair across all inoculation ratios, all bio replicates
ratio_fast_df <- competition_data %>% 
  filter(day_ID == 'D7', competition != '100') %>% 
  select(salinity, rep, competition, pair, ratio_fast) %>%
  pivot_wider(id_cols = c('pair', 'rep', 'competition'), names_from = salinity, names_prefix = 'salt_', values_from = ratio_fast) %>%
  mutate(#salt_46 = ifelse(pair == '40-II.F2', 1-salt_46, salt_46),
         #salt_61 = ifelse(pair == '40-II.F2', 1-salt_61, salt_61),
         diff_31_61 = salt_61 - salt_31,
         diff_16_61 = salt_61 - salt_16) 

mean_ratio_fast <- ratio_fast_df %>%
  group_by(pair) %>%
  summarize(mean_diff_6131 = mean(diff_31_61, na.rm = T),
            sd_diff_6131 = sd(diff_31_61, na.rm = T),
            mean_diff_6116 = mean(diff_16_61, na.rm = T),
            sd_diff_6116 = sd(diff_16_61, na.rm = T)) %>%
  mutate(pair_name = pair_abbreviations[pair])
#Note: In the 40-II.F2 competition the identity of the faster grower switches at 46 and 61 relative to 16 and 31

ggplot(mean_ratio_fast) +
  geom_point(aes(y = pair_name, x = mean_diff_6131, color = pair_name), size = 4, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = pair_name, xmin = mean_diff_6131 - sd_diff_6131, xmax = mean_diff_6131 + sd_diff_6131, color = pair_name), 
                linewidth = 1, width = 0.4, position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(y = '', x = 'Change in relative abundance \n from 31 to 61 g/L', color = 'Pair') +
  coord_cartesian(xlim = c(-0.8, 1)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = color_list_2) +
  theme(text = element_text(size = 36),
        axis.text = element_text(size = 20),
        legend.position = 'none')

ggsave(paste0(fig_dir, '/Fig3F.pdf'), width = 6, height = 5)

##### Development of growth rate as a function of salinity - Fig 3C/E + S14 ####################

# Dataframe with the measured growth rate of all competing species
# at the experimentally used salinities
growth_map_df <- average_growth_df %>% 
  filter(isolate_id %in% c('2', '25', '6', '32', 'II.F2', '40', '8', '39'),
         well_salt %in% c(15, 30, 45, 65)) %>%
  mutate(salinity = case_when(
    well_salt == 15 ~ "16",
    well_salt == 30 ~ "31",
    well_salt == 45 ~ "46",
    well_salt == 65 ~ "61"
  )) %>%
  select(isolate_id, salinity, growthrate, carr_cap) 

for(pair_id in unique_pairs){

  # Pair - community growth rate ; using experimentally measured strain abundances
  # Compute abundance weighted mean growth rate at a given salinity
    cgr_df <- competition_data  %>% 
      mutate(salinity = as.character(salinity)) %>%
      filter(day_ID == 'D7', pair == pair_id) %>%
      left_join(growth_map_df, by = c("salinity", "strain1" = "isolate_id")) %>%
      left_join(growth_map_df, by = c("salinity", "strain2" = "isolate_id"), suffix = c("1", "2")) %>%
      rowwise() %>%
      mutate(pair_gr = (ratio_strain1*growthrate1 + ratio_strain2*growthrate2),
             pair_carrcap = (ratio_strain1*carr_cap1 + ratio_strain2*carr_cap2),
             fast_strain = factor(faster_pair[pair]) ) %>%
      mutate(salinity = as.numeric(salinity)) %>%
      group_by(salinity) %>%
      summarize(mean_pair_gr = mean(pair_gr, na.rm = T),
                sd_pair_gr = sd(pair_gr, na.rm = T),
                mean_pair_carrcap = mean(pair_carrcap, na.rm = T),
                sd_pair_carrcap = sd(pair_carrcap, na.rm = T))
    
    ggplot(average_growth_df %>% 
             filter(isolate_id %in% str_split_1(pair_id, '-') )  %>%
             mutate(strain_label = species_abbreviations[isolate_id]) ) +
      geom_ribbon(aes(x = well_salt, ymin = growthrate - gr_sd, 
                      ymax = growthrate + gr_sd, fill = strain_label), alpha = 0.3, linewidth = 2) +
      geom_line(aes(x = well_salt, y = growthrate, color = isolate_id), linewidth = 0.5, alpha = 0.3, show.legend = F )+ 
      geom_point(data = cgr_df, aes(x = salinity, y = mean_pair_gr), size = 4, show.legend = F )+ 
      geom_errorbar(data = cgr_df, aes(x = salinity, ymin = mean_pair_gr - sd_pair_gr, 
                                        ymax = mean_pair_gr + sd_pair_gr), width = 2, show.legend = F )+ 
      geom_line(data = cgr_df, aes(x = salinity, y = mean_pair_gr), linewidth = 2, show.legend = F )+ 
      labs(x = 'Salinity (g/L)', y = 'Growthrate (1/h)', shape = 'Inoculation \nratio', 
           fill = "Isolate", color = "Isolate", linetype = 'Inoculation ratio') +
      scale_fill_manual(values = species_colors_byname) +
      scale_color_manual(values = species_colors) +
      coord_cartesian(ylim = c(0, 1.4)) +
      theme(axis.title = element_text(size = 25), 
            legend.position = c(0.8, 0.8))
    
    ggsave(paste0(fig_dir, '/CGR_', pair_id, '.pdf'), width = 8, height = 5) 
    
}



## Summary figure with change in CGR relative to 31 g/L - Fig 3G
# reshape the dataframe for plotting; select only D7; calculate average per pair across all inoculation ratios, all bio replicates
cgr_df <- competition_data  %>% 
  mutate(salinity = as.character(salinity)) %>%
  filter(day_ID == 'D7', competition != '100') %>%
  left_join(growth_map_df, by = c("salinity", "strain1" = "isolate_id")) %>%
  left_join(growth_map_df, by = c("salinity", "strain2" = "isolate_id"), suffix = c("1", "2")) %>%
  rowwise() %>%
  mutate(pair_gr = (ratio_strain1*growthrate1 + ratio_strain2*growthrate2),
         mean_gr = (growthrate1 + growthrate2)/2)  %>% 
  select(salinity, rep, competition, pair, pair_gr, mean_gr) %>%
  pivot_wider(id_cols = c('pair', 'rep', 'competition'), names_from = salinity, names_prefix = 'salt_', values_from = c('pair_gr', 'mean_gr')) %>%
  mutate(diff_31_61 = pair_gr_salt_61 - pair_gr_salt_31,
         iso_diff_31_61 = mean_gr_salt_61 - mean_gr_salt_31,
         rel_diff_31_61 = (diff_31_61 - iso_diff_31_61)/mean_gr_salt_31)

mean_cgr_df <- cgr_df %>%
  group_by(pair) %>%
  summarize(mean_diff_6131 = mean(diff_31_61, na.rm = T),
            sd_diff_6131 = sd(diff_31_61, na.rm = T),
            mean_rel_diff_6131 = mean(rel_diff_31_61, na.rm = T),
            sd_rel_diff_6131 = sd(rel_diff_31_61, na.rm = T),
            iso_mean_diff_6131 = mean(iso_diff_31_61, na.rm = T)) %>%
  mutate(pair_name = pair_abbreviations[pair])
#Note: In the 40-II.F2 competition the identity of the faster grower switches at 46 and 61 relative to 16 and 31

ggplot(mean_cgr_df) +
  #geom_point(aes(y = pair_name, x = mean_rel_diff_6131, color = pair_name), size = 4, position = position_dodge(width = 1)) +
  #geom_errorbar(aes(y = pair_name, xmin = mean_rel_diff_6131 - sd_rel_diff_6131, xmax = mean_rel_diff_6131 + sd_rel_diff_6131, color = pair_name), 
  #              linewidth = 1, width = 0.4, position = position_dodge(width = 1)) +
  geom_point(aes(y = pair_name, x = mean_diff_6131, color = pair_name), size = 4, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = pair_name, xmin = mean_diff_6131 - sd_diff_6131, xmax = mean_diff_6131 + sd_diff_6131, color = pair_name), 
                linewidth = 1, width = 0.4, position = position_dodge(width = 1)) +
  geom_point(aes(y = pair_name, x = iso_mean_diff_6131, color = pair_name), size = 4, position = position_dodge(width = 1), alpha = 0.5, shape = 'triangle') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(y = '', x = 'Change in growth rate \n from 31 to 61 g/L', color = 'Pair') +
  coord_cartesian(xlim = c(-0.8, 1)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = color_list_2) +
  theme(text = element_text(size = 36),
        axis.text = element_text(size = 20),
        legend.position = 'none')

ggsave(paste0(fig_dir, '/Fig3G_CGR.pdf'), width = 6, height = 5)



