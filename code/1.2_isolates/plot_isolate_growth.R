################################################################################
## Plot growth related figures
## Author: JS Huisman
################################################################################

library(tidyverse)
library(readxl)
library(gcplyr)
library(patchwork)

proc_data_dir = '../data/1.2_isolates/'
fig_dir = '../figures/1.2_isolates/'

theme_set(theme_minimal() + theme(text = element_text(size = 20)))

color_list <- c('Estuary' = '#A72626', 
                'Marine 2' = '#F5BD23', 
                'Marine 1' = '#356FEA', 
                'Brackish' = '#004D07')


color_list_2 = c('#C0E49B', '#F5BD23', '#82E2F3', '#B77406', '#A72626', '#6E969B', 
                 '#741B7B', '#004D07', '#356FEA', '#000000')

################################################################################

all_growth_df <- read_csv(paste0(proc_data_dir, 'all_growthdata.csv'))
average_growth_df <- read_csv(paste0(proc_data_dir, 'average_growthrates.csv'))

# Number of growth measurements per isolate
n_growth_meas <- all_growth_df %>% 
  group_by(isolate_id) %>%
  dplyr::count() %>%
  mutate(n = floor(n/12))

n3_isolate_ids <- n_growth_meas[n_growth_meas$n>=3, 'isolate_id']$isolate_id

# select only isolates with at least 3 independent measurements
n3_all_growth_df <- all_growth_df %>%
  filter(isolate_id %in% n3_isolate_ids)

n3_mean_growth_df <- average_growth_df %>%
  filter(isolate_id %in% n3_isolate_ids)

################################################################################
## Plot some example growth curves - Fig. 1C ####

isolate_selection = c('29', '16', '6')

ggplot(average_growth_df %>% filter(isolate_id %in% isolate_selection ) ) +
  geom_ribbon(aes(x = well_salt, ymin = growthrate - gr_sd,
                    ymax = growthrate + gr_sd, fill = genus), alpha = 0.5) +
  geom_line(aes(x = well_salt, y = growthrate, color = genus), linewidth = 2) + 
  labs(x = 'Salinity (g/L)', y = 'Growthrate (1/h)', color = '', fill = '') +
  scale_color_manual(values = c('#C0E49B', '#6E969B', '#741B7B')) +
  scale_fill_manual(values = c('#C0E49B', '#6E969B', '#741B7B')) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 15))

ggsave(paste0(fig_dir, '/Fig1B.pdf'), height = 4, width = 6)

################################################################################
## Plot all growth curves - Fig. S9 ####

ggplot(n3_mean_growth_df %>% filter(!is.na(genus))) +
  geom_ribbon(aes(x = well_salt, ymin = growthrate - gr_sd, ymax = growthrate + gr_sd,
                fill = isolate_id), alpha = 0.5, show.legend = F) + 
  #geom_line(aes(x = well_salt, y = max_percap, color = isolate_id, group = interaction(isolate_id, replicate)), 
  #          alpha = 0.5, show.legend = F) + 
  geom_line(aes(x = well_salt, y = growthrate, color = isolate_id), linewidth = 1, show.legend = F) + 
  facet_wrap(vars(family), nrow = 5) +
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(x = 'Salinity', y = 'Growth rate', color = 'Replicate') +
  #scale_color_manual(color_list_2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = 'bottom')


ggsave(paste0(fig_dir, 'FigS8_growthcurves_by_family.pdf'), height = 10, width = 14)


###############################################################################
#### Plot correlation between measured growth rate and estimated rrn - Fig. S16 ####

average_growth_df <- read_csv('../data/1.2_isolates/average_growthrates.csv')
isolate_rrndb <- read_csv(file = '../data/1.2_isolates/isolate_rrndb.csv')
isolate_info <- read_csv('../data/1.2_isolates/paper_isolate_ids.csv')

# find rmax for each isolate
max_df <- average_growth_df %>%
  filter(isolate_id %in% n3_isolate_ids) %>%
  select("id" = isolate_id, "gr_salt" = well_salt, 'gr' = growthrate, gr_sd) %>%
  group_by(id) %>%
  summarise(gr_max = max(gr, na.rm = T),
            gr_sd = gr_sd[gr == gr_max],
            gr_salt = factor('max',  levels = c(5, 35, 'max') )) %>%
  rename(gr_max = 'gr')

# Combine growth data at 5 and 35 g/L (only isolates with >= 3 measurements) and rrndb
isolate_growth_rrndb = average_growth_df %>%
  filter(isolate_id %in% n3_isolate_ids) %>%
  select("id" = isolate_id, "gr_salt" = well_salt, 'gr' = growthrate, gr_sd) %>%
  filter(gr_salt %in% c(5, 35)) %>%
  mutate(gr_salt = factor(gr_salt, levels = c(5, 35,  'max'))) %>%
  bind_rows(max_df) %>%
  inner_join(isolate_rrndb, by = 'id') %>%
  left_join(isolate_info %>% select('id' = isolate_id, inoc), by = 'id')  %>%
  mutate(oceanic = ifelse(inoc == 'charles', 'Brackish', 'Marine'))

# color by isolate origin
ggplot(isolate_growth_rrndb ) +
  geom_point(aes(x = copy_number_NCBI, y = gr, color = oceanic, shape = oceanic), size = 5) +
  geom_errorbar(aes(x = copy_number_NCBI, ymin = pmax(gr - gr_sd, 0), ymax = gr + gr_sd, color = oceanic)) +
  stat_smooth(aes(x = copy_number_NCBI, y = gr, group = oceanic, linetype = oceanic, color = oceanic, fill = oceanic), method = lm) +
  facet_wrap(vars(gr_salt), labeller = labeller(gr_salt = c('5' = '5 (g/L)', '35' = '35 (g/L)', 'max' = 'Sopt')), nrow = 2) +
  scale_color_manual(values = c('#A72626', '#6E969B')) +
  scale_fill_manual(values = c('#A72626', '#6E969B')) +
  labs(x = 'rRNA copy number', y = 'Growth rate (1/h)', color = 'Isolate origin', fill = 'Isolate origin',
       linetype = 'Isolate origin', shape = 'Isolate origin') +
  theme(legend.position = 'bottom',
        panel.spacing.x = unit(2, 'lines'))

ggsave(paste0(fig_dir, '/FigS15_growth_vs_mcn_bysalinity.pdf'), width = 12, height = 10)


