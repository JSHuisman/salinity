###########################################################
## Plot daily dilutions results
## Author: JS Huisman
###########################################################

library(tidyverse)
library(ggplot2)
library(readxl)
library(patchwork)

proc_data_dir = '../data/1_community_experiment/'

# global rules for plotting
theme_set(theme_minimal() + theme(text = element_text(size = 20)))

###########################################################

color_list <- c('Estuary' = '#A72626', 
                'Marine 2' = '#F5BD23', 
                'Marine 1' = '#356FEA', 
                'Brackish' = '#004D07',
                'blank' = '#82E2F3') #6E969B

color_salt2 <- c('16' = '#7EFAD5',
                '31' = '#1E88E5',
                '46' = '#6B02B5')

##### OD #################################################
## read OD data
all_od_data <- read_csv(paste0(proc_data_dir, 'all_od_data.csv'),
                        col_types = 'ffdfd')

blank_df <- all_od_data %>%
  filter(cond == 'blank') %>%
  group_by(day, temp, salt) %>%
  summarise(mean_blank = mean(od),
            .groups = 'keep')

od_data <- all_od_data %>%
  filter(cond != 'blank') %>%
  left_join(blank_df, by = c('day', 'temp', 'salt')) %>%
  mutate(clean_od = od - mean_blank,
         inoc_map = case_when(
           cond == 'charles' ~ 'Brackish',
           cond == 'ica' ~ 'Estuary',
           cond == 'nahant' ~ 'Marine 1',
           cond == 'fucus' ~ 'Marine 2',
           cond == 'blank' ~ 'Blank'
         ),
         temp = ifelse(temp == '15C', '15C', '20C'),
         cycle = gsub('d', 'C', day) )

## plot day 6; split by salinity, colour by temp and condition - lines
od_salt_D6_plot <- ggplot(od_data %>% filter(day == 'd6') %>%
         mutate(cond = factor(cond, levels = c('charles', 'ica', 'nahant', 'fucus')),
                salt_num = as.numeric(salt)) )+
  geom_smooth(aes(x = salt_num, y = clean_od), method='lm', color = 'black', alpha = 0.2) +
  geom_line(aes(x = salt_num, y = clean_od, group = interaction(temp, cond), linetype = temp, color = inoc_map), linewidth = 0.5, stat = "summary", fun = "mean") +
  geom_point(aes(x = salt_num, y = clean_od, colour = inoc_map, shape = temp), size = 3) +
  facet_wrap(vars(inoc_map), nrow = 1) +
  scale_color_manual(values = color_list) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  coord_cartesian(ylim = c(0, 0.9)) +
  labs(x = 'Salinity (g/L)', y = 'OD', colour = 'Inoculum', linetype = 'Temperature', shape  = 'Temperature') +
  theme_bw() + 
  theme(legend.position = 'right',
        text = element_text(size = 20),
        strip.text.x.top = element_blank(),
        panel.spacing = unit(1.5, 'lines'))

od_salt_D6_plot
ggsave('../figures/1_community_experiment/odvssalt_d6.png', width  = 12, height = 4)

## Fig S4A: just day 6;  colour by condition
ggplot(od_data %>% filter(day == 'd6') %>%
                            mutate(cond = factor(cond, levels = c('charles', 'ica', 'nahant', 'fucus')),
                                   salt_num = as.numeric(salt)) )+
  geom_smooth(aes(x = salt_num, y = clean_od, color = inoc_map), method='lm', se = F, alpha = 0.2) +
  geom_point(aes(x = salt_num, y = clean_od, colour = inoc_map, shape = temp), size = 3) +
  scale_color_manual(values = color_list) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  coord_cartesian(ylim = c(0, 0.9)) +
  labs(x = 'Salinity (g/L)', y = 'OD', colour = 'Inoculum', linetype = 'Temperature', shape  = 'Temperature') +
  theme(legend.position = 'right',
        text = element_text(size = 20),
        strip.text.x.top = element_blank(),
        panel.spacing = unit(1.5, 'lines'))

ggsave('../figures/1_community_experiment/FigS4A.pdf', width  = 6, height = 4)

# Load isolate growth data
average_growth_df <- read_csv('../data/1.2_isolates/average_growthrates.csv')
isolate_info <- read_csv('../data/1.2_isolates/paper_isolate_ids.csv')

isolate_growth_df <- average_growth_df %>% 
  left_join(isolate_info, by = c('isolate_id', 'class', 'order', 'family', 'genus', 'species')) %>% 
  filter(!is.na(inoc),
         n_growth_meas >= 3)

mean_isolate_carrcap_df <- isolate_growth_df %>%
  pivot_wider(id_cols = c('isolate_id', 'inoc'), names_from = well_salt, 
              names_prefix = 'salt_', values_from = carr_cap) %>%
  mutate(diff_30_45 = salt_45 - salt_30,
         diff_15_30 = salt_30 - salt_15,
         diff_15_45 = salt_45 - salt_15) %>%
  group_by(inoc) %>%
  summarize(mean_diff_4530 = mean(diff_30_45, na.rm = T),
            sd_diff_4530 = sd(diff_30_45, na.rm = T),
            mean_diff_4515 = mean(diff_15_45, na.rm = T),
            sd_diff_4515 = sd(diff_15_45, na.rm = T),
            mean_diff_3015 = mean(diff_15_30, na.rm = T),
            sd_diff_3015 = sd(diff_15_30, na.rm = T)) %>%
  mutate(inoc_map = case_when(
    inoc == 'charles' ~ 'Brackish',
    inoc == 'ica' ~ 'Estuary',
    inoc == 'nahant' ~ 'Marine 1',
    inoc == 'fucus' ~ 'Marine 2'
  )) %>%
  pivot_longer(cols = -c('inoc', 'inoc_map'), names_pattern = '(.*)_diff_(.*)',
               names_to = c('metric', 'window')) %>%
  pivot_wider(id_cols = c('inoc', 'inoc_map', 'window'), names_from = metric)

mean_community_cc <- od_data %>% filter(day == 'd6') %>%
  group_by(day, temp, salt, cond) %>%
  mutate(rep = row_number(),
         inoc = factor(cond, levels = c('charles', 'ica', 'nahant', 'fucus')),
         salt_num = as.numeric(salt)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('rep', 'inoc', 'temp', 'inoc_map'), names_from = salt_num,
              names_prefix = 'salt_', values_from = clean_od) %>%
  mutate(diff_30_45 = salt_46 - salt_31,
         diff_15_30 = salt_31 - salt_16,
         diff_15_45 = salt_46 - salt_16) %>%
  group_by(inoc, inoc_map) %>%
  summarize(mean_diff_4530 = mean(diff_30_45, na.rm = T),
            sd_diff_4530 = sd(diff_30_45, na.rm = T),
            mean_diff_4515 = mean(diff_15_45, na.rm = T),
            sd_diff_4515 = sd(diff_15_45, na.rm = T),
            mean_diff_3015 = mean(diff_15_30, na.rm = T),
            sd_diff_3015 = sd(diff_15_30, na.rm = T)) %>%
  pivot_longer(cols = -c('inoc', 'inoc_map'), names_pattern = '(.*)_diff_(.*)',
               names_to = c('metric', 'window')) %>%
  pivot_wider(id_cols = c('inoc', 'inoc_map', 'window'), names_from = metric)

# Put the IGR and CGR together for plotting
plot_df = bind_rows(mean_isolate_carrcap_df, mean_community_cc, .id = 'metric')  %>% 
  mutate(metric = factor(ifelse(metric == 1, 'Isolates', 'Community'), levels = c('Isolates', 'Community'))) %>%
  filter(window == '4515')

ggplot(plot_df) +
  geom_point(aes(x = inoc_map, y = mean, color = inoc_map, shape = metric, alpha = metric), size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = inoc_map, ymin = mean - sd, ymax = mean + sd, color = inoc_map, alpha = metric), position = position_dodge(width = 0.5), 
                linewidth = 2, width = 0) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = '', y = latex2exp::TeX('\\Delta carrying capacity'), color = 'Inoculum', alpha = '', shape = '') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  scale_alpha_discrete(range = c(0.4, 1)) +
  scale_color_manual(values = color_list) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.position = 'inside',
        legend.position.inside = c(0.8, 0.2)) +
  guides(color = "none")

ggsave('../figures/1_community_experiment/FigS4C.pdf', width = 6, height = 4)


# family level difference in carrying capacity?

color_list_families = c("Pseudoalteromonadaceae" = '#6E969B', 
                        "Moraxellaceae"  = '#C0E49B', "Pseudomonadaceae"  = '#82E2F3', 
                        "Flavobacteriaceae" = '#F5BD23',"Shewanellaceae"= '#741B7B',
                        "Oceanospirillaceae" = '#004D07',
                        "Rhodobacteraceae"  = '#A72626',
                        "Bacillaceae" = '#B77406',  "Colwelliaceae" = '#356FEA')

carrcap_family_df <- average_growth_df %>%
  filter(n_growth_meas >= 3,
         !is.na(family)) %>%
  group_by(class, order, family, genus, species, isolate_id) %>%
  summarize(max_cc = max(carr_cap, na.rm = T),
            cc_diff = carr_cap[well_salt == 45] -  carr_cap[well_salt == 15]) %>%
  group_by(family) %>%
  summarize(n_isolate = n(),
            mean_max_cc = mean(max_cc, na.rm = T),
            sd_max_cc = sd(max_cc, na.rm = T),
            mean_cc_diff = mean(cc_diff, na.rm = T),
            sd_cc_diff = sd(cc_diff, na.rm = T)) %>%
  filter(n_isolate >= 3) %>%
  mutate(ordered_family_id = factor(family, levels = c("Pseudomonadaceae" , "Moraxellaceae", "Oceanospirillaceae" ,
                                                       "Shewanellaceae","Pseudoalteromonadaceae" ,"Colwelliaceae", 
                                                       "Flavobacteriaceae", "Bacillaceae" ,  "Rhodobacteraceae"  )))

max_cc_plot <- ggplot(carrcap_family_df) +
  geom_point(aes(y = ordered_family_id, x = mean_max_cc, color = ordered_family_id), size = 6, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = ordered_family_id, xmin = pmax(mean_max_cc - sd_max_cc, 0), 
                    xmax = mean_max_cc + sd_max_cc, color = ordered_family_id), 
                linewidth = 2, width = 0, position = position_dodge(width = 1)) +
  labs(y = '', x = latex2exp::TeX('Max OD600'), color = 'Family') +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = color_list_families) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'none')

max_cc_plot

cc_diff_plot <- ggplot(carrcap_family_df) +
  geom_point(aes(y = ordered_family_id, x = mean_cc_diff, color = ordered_family_id), size = 6, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = ordered_family_id, xmin = mean_cc_diff - sd_cc_diff, 
                    xmax = mean_cc_diff + sd_cc_diff, color = ordered_family_id), 
                linewidth = 2, width = 0, position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(y = '', x = latex2exp::TeX('\\Delta OD600 (45-15 g/L)'), color = 'Family') +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_color_manual(values = color_list_families) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'none')

cc_diff_plot


(max_cc_plot + theme(axis.text.y = element_text(size = 20))) +
  (cc_diff_plot + theme(axis.text.y = element_blank()))
ggsave('../figures/1_community_experiment/FigS4B.pdf', width = 9, height = 4)

