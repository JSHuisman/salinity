######################################################################################################################
## Plot OD data from pairwise competitions
## Author: JS Huisman
######################################################################################################################

library(tidyverse)
library(readxl)
library(lubridate)

proc_data_dir = '../data/3_pairwise_competition'
fig_dir = '../figures/3_pairwise_competition'

theme_set(theme_minimal() + theme(text = element_text(size = 20)))

od_data <- read_csv(file = paste0(proc_data_dir, '/od_data.csv'))

salinity_colors <- c('16' = '#7EFAD5',
                     '31' = '#1E88E5',
                     '46' = '#6B02B5',
                     '61' = 'darkblue')

##### OD data over the experiment timecourse ###############################################

# plot OD over time for each monoculture
ggplot(od_data %>% filter(competition ==100)) +
  geom_point(aes(y = od, x = day_ID, color = as.factor(salinity))) +
  geom_line(aes(y = od, x = day_ID, color = as.factor(salinity), linetype = as.factor(bio_rep),
                group = interaction(salinity, competition, bio_rep))) +
  facet_wrap(vars(pair), nrow = 2) +
  scale_color_manual(values = salinity_colors) +
  labs(x = "Cycle", y = "OD", color = "Salinity (g/L)", linetype = "Biological \nreplicate")

ggsave(paste0(fig_dir, '/OD_over_time_monoculture.png'), width = 10, height = 6)

# plot OD at D7 for each monoculture
ggplot( )+
  geom_point(data = od_data %>% filter(competition ==100, day_ID %in% c( 'D7')),
             aes(y = od, x = as.factor(salinity), color = as.factor(salinity)), size = 3) +
  facet_wrap(vars(pair), nrow = 2) +
  scale_color_manual(values = salinity_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Salinity (g/L)", y = "OD", color = "Salinity (g/L)") +
  theme(legend.position = "bottom",
        text = element_text(size = 20))

ggsave(paste0(fig_dir, '/OD_D7_monoculture.png'), width = 10, height = 5)


# plot OD over time for each pair
ggplot(od_data %>% filter(competition !=100)) +
  geom_point(aes(y = od, x = day_ID, color = as.factor(salinity), shape = as.factor(bio_rep))) +
  geom_line(aes(y = od, x = day_ID, color = as.factor(salinity), linetype = competition,
                group = interaction(salinity, competition, bio_rep))) +
  facet_wrap(vars(pair)) +
  scale_color_manual(values = salinity_colors) +
  labs(x = "Cycle", y = "OD", color = "Salinity (g/L)", linetype = "Inoculation \nratio", 
       shape = "Biological \nreplicate")

ggsave(paste0(fig_dir, '/OD_over_time_pair.png'), width = 10, height = 8)


# plot OD at D7 for each pair
ggplot( )+
  geom_point(data = od_data %>% filter(competition !=100, day_ID %in% c( 'D7')),
             aes(y = od, x = as.factor(salinity), color = as.factor(salinity), shape = competition), size = 2) +
  facet_wrap(vars(pair), nrow = 2) +
  scale_color_manual(values = salinity_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Salinity (g/L)", y = "OD", color = "Salinity (g/L)", shape = "Inoculation ratio") +
  theme(legend.position = "bottom",
        text = element_text(size = 20))

ggsave(paste0(fig_dir, '/OD_D7_pair.png'), width = 10, height = 5)

