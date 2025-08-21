######################################################################################################################
## Estimate and plot community growth rates
## Author: JS Huisman
######################################################################################################################
library(tidyverse)
library(ggplot2)
library(patchwork)

proc_data_dir = '../data/1.2_isolates/'
fig_dir = '../figures/1_community_experiment/'

seq_data_dir = '../data/1_community_experiment'

theme_set(theme_minimal() + 
            theme(text = element_text(size = 20)))

color_list <- c('Estuary' = '#A72626', 
                'Marine 2' = '#F5BD23', 
                'Marine 1' = '#356FEA', 
                'Brackish' = '#004D07')

color_list_2 = c('#6E969B', '#004D07', '#82E2F3', '#F5BD23', '#741B7B',  '#C0E49B',  '#A72626',
                 '#B77406',  '#356FEA', '#000000')

######################################################################################################################

# Load growth data
all_results_df <- read_csv(paste0(proc_data_dir, 'all_growthdata.csv'))
average_growth_df <- read_csv(paste0(proc_data_dir, 'average_growthrates.csv'))

# Load sequencing data, blast, mcn
# ps_noRare was preprocessed to ASVs that appear only once across all samples
rel_abund = read_csv(paste0(seq_data_dir, '/rel_abund_df_norare.csv'))

rel_abund <- rel_abund %>%
  mutate(inoc = factor(inoc, levels = c('charles', 'ica', 'nahant', 'fucus')))

mcn_df = read_csv(paste0(seq_data_dir, '/mcn_df_norare.csv')) %>%
  mutate(inoc = factor(inoc, levels = c('charles', 'ica', 'nahant', 'fucus')))

######################################################################################################################
# Map ASVs to isolates and their growth rates ####

asv_isolate_map <- read_csv("../data/1.2_isolates/asv_isolate_map.csv")

# split the list of isolates matching an ASV into rows
long_asv_isolate_map <- asv_isolate_map %>%
  rowwise() %>%
  mutate(matching_isolates = list(str_split_1(matching_isolates, '/')) ) %>%
  unnest_longer(matching_isolates) %>%
  select("OTU" = QueryID, "isolate_id" = matching_isolates)

# for each isolate, determine max growth rate and max density across all salinities
max_gr_cc_df <- average_growth_df %>%
  ungroup() %>%
  select(isolate_id, "gr_salt" = well_salt, 'gr' = growthrate, 'cc' = carr_cap, n_growth_meas) %>%
  group_by(isolate_id) %>%
  summarise(max_gr = max(gr, na.rm = T),
            s_max_gr = gr_salt[which.max(gr)],
            max_cc = max(cc, na.rm = T),
            s_max_cc = gr_salt[which.max(cc)],
            n_growth = min(n_growth_meas))

# for each ASV determine the isolate with most growth measurements (from all equivalently best matching)
# this will be the stand-in isolate id of the ASV (for plotting)
asv_best_isolate <- long_asv_isolate_map %>%
  left_join(max_gr_cc_df, by = 'isolate_id') %>%
  group_by(OTU) %>%
  arrange(-n_growth) %>%
  filter(row_number()==1) %>%
  arrange(OTU) %>%
  ungroup() %>%
  select(OTU, isolate_id, n_growth) %>%
  complete(rel_abund['OTU']) %>%
  mutate(meas_isolate_id = ifelse(is.na(n_growth), NA, isolate_id))
  

# per asv, average across the growth rates of all best-matching isolates (include only isolates with at least 3 measurements)
asv_growth_df = average_growth_df %>%
  ungroup() %>%
  select(isolate_id, "gr_salt" = well_salt, 'gr' = growthrate, 'cc' = carr_cap, n_growth_meas) %>%
  filter(gr_salt %in% c(15, 30, 35, 45)) %>%
  pivot_wider(id_cols = c('isolate_id'), names_from = 'gr_salt', values_from = c('gr', 'cc')) %>% 
  right_join(long_asv_isolate_map, by = 'isolate_id') %>%
  left_join(max_gr_cc_df, by = 'isolate_id') %>%
  filter(n_growth >= 3) %>%
  group_by(OTU) %>%
  summarize(across(-c(isolate_id, n_growth), ~ mean(.x, na.rm = TRUE)),
            n_growth = sum(n_growth, na.rm = T)) 

# combine growth info with rel. abundance
asv_isolate_abund_gr = asv_growth_df %>%
  right_join(rel_abund, by = 'OTU') %>%
  left_join(asv_best_isolate, by = "OTU") %>%
  mutate(salt = gsub('%', '', salt) ) 
  

###########################################################
# plot fraction of the communities at C6 covered by isolates - Supp Fig. S6

isolate_abund_15C = ggplot() +
  geom_bar(data = asv_isolate_abund_gr %>% filter(day == 'D6', temp == '15C'),
           aes(x = rep, y = rel_abund, fill = meas_isolate_id), color = 'black', linewidth = 0.1, stat = 'identity') +
  facet_grid(rows = vars(inoc), cols = vars(salt)) +
  labs(x = 'Replicate', y = 'Relative abundance', fill = 'Representative isolate', title = '15C') +
  scale_fill_viridis_d(option = 'viridis', na.value="grey") +
  theme(legend.position = 'bottom')

isolate_abund_roomT = ggplot(asv_isolate_abund_gr %>% filter(day == 'D6', temp == 'roomT') ) +
  geom_bar(aes(x = rep, y = rel_abund, fill = meas_isolate_id), color = 'black', linewidth = 0.1, stat = 'identity', show.legend = T) +
  labs(x = 'Replicate', y = 'Relative abundance', fill = 'Representative isolate', title = '20C') +
  facet_grid(rows = vars(inoc), cols = vars(salt)) +
  scale_fill_viridis_d(option = 'viridis', na.value="grey")

isolate_abund_15C + isolate_abund_roomT +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'right')

ggsave(paste0(fig_dir, 'FigS6.pdf'), width = 16)

## Which percentage of communities is covered
abund_covered_df <- asv_isolate_abund_gr %>% 
  filter(day == 'D6',
         !is.na(meas_isolate_id),
         !is.na(rel_abund)) %>%
  select(OTU, rep, day, inoc, salt, temp, GenusSpecies, copy_number_NCBI, rel_abund, isolate_id, meas_isolate_id) %>%
  group_by(inoc, salt, temp, rep) %>%
  #arrange(-rel_abund) %>% View()
  summarize(abund_covered = sum(rel_abund)) %>% 
  arrange(abund_covered)

abund_covered_df

abund_covered_df %>%
  ungroup() %>%
  summarize(mean_covered = mean(abund_covered)*100,
            median_covered = median(abund_covered)*100,
            sd_covered = sd(abund_covered)*100)

## Which isolates are used to cover the community abundance?
asv_isolate_abund_gr %>% 
  filter(day == 'D6',
         !is.na(meas_isolate_id),
         !is.na(rel_abund)) %>%
  pull(meas_isolate_id) %>% unique() %>% sort()

######################################################################################################################
##### community growth rates ######

# compute weighted mean growth rates of the experimental communities
# ascribing each ASV to different growth rates (i.e. measured at 15, 30, 45 g/L or their best growth rate across all salinities)
community_gr = asv_isolate_abund_gr %>% 
  filter(day == 'D6') %>%
  group_by(rep, day, inoc, salt, temp) %>%
  summarise(community_gr_15 = weighted.mean(gr_15, rel_abund, na.rm = T),
            community_gr_30 = weighted.mean(gr_30, rel_abund, na.rm = T) ,
            community_gr_35 = weighted.mean(gr_35, rel_abund, na.rm = T) ,
            community_gr_45 = weighted.mean(gr_45, rel_abund, na.rm = T) ,
            community_cc_15 = weighted.mean(cc_15, rel_abund, na.rm = T),
            community_cc_30 = weighted.mean(cc_30, rel_abund, na.rm = T) ,
            community_cc_35 = weighted.mean(cc_35, rel_abund, na.rm = T) ,
            community_cc_45 = weighted.mean(cc_45, rel_abund, na.rm = T),
            comm_max_gr = weighted.mean(max_gr, rel_abund, na.rm = T),
            comm_max_cc = weighted.mean(max_cc, rel_abund, na.rm = T)) %>%
  pivot_longer(cols = starts_with('community_'), names_pattern = '(.*)(_15|_30|_35|_45)', 
               names_to = c('type', 'gr_salt')) %>%
  mutate(gr_salt = as.numeric(gsub('_', '', gr_salt)),
         salt_num = as.numeric(gsub('%', '', salt)),
         inoc_map = case_when(
           inoc == 'charles' ~ 'Brackish',
           inoc == 'ica' ~ 'Estuary',
           inoc == 'nahant' ~ 'Marine 1',
           inoc == 'fucus' ~ 'Marine 2'
         )) %>%
  pivot_wider(names_from = type) %>%
  dplyr::rename(gr= community_gr, carr_cap = community_cc)

# isolate growth at the community salinity
community_gr_by_salt = community_gr %>%
  filter(case_when(
    salt == '16' ~ gr_salt == 15,
    salt == '31' ~ gr_salt == 30,
    salt == '46' ~ gr_salt == 45
  ))

## Map every isolate onto growth at 30 g/L
g30_salt_D6_plot <- ggplot(community_gr %>% filter(gr_salt == 30) ) +
  geom_smooth(aes(x = salt_num, y = gr, color = inoc_map, fill = inoc_map ), method='lm', se = F, linewidth = 2, alpha = 0.2,  show.legend = F) +
  geom_point(aes(x = salt_num, y = gr, color = inoc_map, group = interaction(salt_num, inoc), shape = temp), size = 4,  show.legend = F) +
  #facet_wrap(vars(inoc_map), ncol = 4) +
  coord_cartesian(ylim = c(0,1.4)) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  labs(x = 'Salinity g/L', y = 'CCI', colour = 'Inoculum', fill = 'Inoculum', linetype = 'Temperature', shape  = 'Temperature') +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  theme(legend.position = 'right',
        panel.spacing = unit(1.5, 'lines'))

g30_salt_D6_plot
ggsave(paste0(fig_dir, 'CCI_30gL.pdf'), width = 5, height = 4)

## Growth at their actual community salinity #####
cgr_salt_D6_plot <- ggplot(community_gr_by_salt  ) +
  geom_smooth(aes(x = salt_num, y = gr, group = inoc, color = inoc_map, fill = inoc_map ), method='lm', alpha = 0.2, se = F, linewidth = 2) +
  geom_point(aes(x = salt_num, y = gr, color = inoc_map, shape = temp, group = interaction(salt_num, inoc)), 
             size = 4, show.legend = F) +
  #facet_wrap(vars(inoc_map), ncol = 4) +
  coord_cartesian(ylim = c(0,1.4)) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  labs(x = 'Salinity g/L', y = 'CGR (1/h)', color = 'Inoculum', fill = 'Inoculum', shape = 'Temperature') +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  theme(panel.spacing = unit(1.5, 'lines'),
        legend.position = 'none')

cgr_salt_D6_plot
ggsave(paste0(fig_dir, 'Fig1E_CGR.pdf'), width = 4, height =4)


# cgr_salt_D6_plot + g30_salt_D6_plot + 
#   plot_layout(guides = 'collect', nrow = 1) &
#   theme(text = element_text(size = 25),
#         axis.text = element_text(size = 20),
#         legend.position = 'bottom',
#         plot.margin = unit(c(1,2,1,1), 'lines'))
# ggsave(paste0(fig_dir, 'CCI30_CGR_serial_dilution.pdf'), width = 10, height = 7)

############################################
mean_CCI_diff <- community_gr %>%
  filter(gr_salt == 30) %>%
  pivot_wider(id_cols = c('rep', 'inoc', 'temp', 'inoc_map'), names_from = salt_num,
              names_prefix = 'salt_', values_from = gr) %>%
  group_by(inoc, inoc_map) %>%
  summarize(mean_16 = mean(salt_16, na.rm = T),
            sd_16 = sd(salt_16, na.rm = T),
            mean_31 = mean(salt_31, na.rm = T),
            sd_31 = sd(salt_31, na.rm = T),
            mean_46 = mean(salt_46, na.rm = T),
            sd_46 = sd(salt_46, na.rm = T)) %>%
  mutate(mean_diff_4530 = mean_46 - mean_31,
         sd_diff_4530 = sqrt(sd_46^2+sd_31^2),
         se_diff_4530 = sd_diff_4530/sqrt(6),
         mean_diff_4515 = mean_46 - mean_16,
         sd_diff_4515 = sqrt(sd_46^2+sd_16^2),
         se_diff_4515 = sd_diff_4515/sqrt(6)
         ) 

## CCI
ggplot(mean_CCI_diff) +
  geom_point(aes(y = inoc_map, x = mean_diff_4515, color = inoc_map), size = 6, position = position_dodge(width = 1)) +
    geom_errorbar(aes(y = inoc_map, xmin = mean_diff_4515 - se_diff_4515, 
                    xmax = mean_diff_4515 + se_diff_4515, color = inoc_map), 
                linewidth = 2, width = 0, position = position_dodge(width = 1)) +
  
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(y = '', x = latex2exp::TeX('$Delta CCI_{45-15}$') , color = 'Pair') +
  coord_cartesian(xlim = c(-0.2, 0.5)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = color_list) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.position = 'none')
ggsave('../figures/1_community_experiment/Fig1F_CCI_summary.pdf', width = 4, height = 4)


###########################################################
## MCN plot in same style
# ggplot(mcn_df %>% filter(day == 'D6', temp == 'roomT') %>%
#          mutate(inoc = factor(inoc, levels = c('charles', 'ica', 'nahant', 'fucus')),
#                 salt_num = as.numeric(gsub('%', '', salt)),
#                 inoc_map = case_when(
#                   inoc == 'charles' ~ 'Brackish',
#                   inoc == 'ica' ~ 'Estuary',
#                   inoc == 'nahant' ~ 'Marine 1',
#                   inoc == 'fucus' ~ 'Marine 2'
#                 )) ) +
#   geom_point(aes(x = salt_num, y = mcn, color = inoc_map), size = 3) +
#   facet_wrap(vars(inoc_map), ncol = 4) +
#   coord_cartesian(ylim = c(0,11)) +
#   scale_x_continuous(breaks = c(16, 31, 46)) +
#   labs(x = 'Salinity (g/L)', y = 'Mean copy number', color = 'Inoculum') +
#   scale_color_manual(values = color_list) +
#   #theme_bw() +
#   theme(legend.position = 'bottom',
#         text = element_text(size = 20),
#         panel.spacing = unit(1.5, 'lines')) 
# 
# ggsave(paste0(fig_dir, 'mcn_afo_salinity.pdf'), width = 8)

######################################################################################################################
# Extract additional features from the growth curves - Fig. S7

gr_traits = average_growth_df %>%
  filter(n_growth_meas >= 3,
         isolate_id != 'I.D4') %>%
  ungroup() %>%
  select(isolate_id, "gr_salt" = well_salt, 'gr' = growthrate) %>%
  group_by(isolate_id) %>%
  summarise(smin = gr_salt[min(which(gr > 0))],
            smax = gr_salt[max(which(gr > 0))],
            sopt = gr_salt[which.max(gr)],
            gr_max = max(gr),
            smin_halfmax = gr_salt[min(which(gr > gr_max/2))],
            smax_halfmax = gr_salt[max(which(gr > gr_max/2))],
            range_halfmax = smax_halfmax - smin_halfmax) 

community_gr_traits <- gr_traits %>%
  left_join(asv_isolate_abund_gr, by = 'isolate_id') %>% 
  filter(!is.na(rel_abund)) %>%
  filter(day == 'D6') %>% 
  group_by(rep, day, inoc, salt, temp) %>%
  summarise(community_smin = weighted.mean(smin, rel_abund, na.rm = T),
            community_smax = weighted.mean(smax, rel_abund, na.rm = T) ,
            community_sopt = weighted.mean(sopt, rel_abund, na.rm = T) ,
            community_maxgr = weighted.mean(gr_max, rel_abund, na.rm = T) ,
            community_range_halfmax = weighted.mean(range_halfmax, rel_abund, na.rm = T) ) %>%
  mutate(salt_num = as.numeric(gsub('%', '', salt)) ,
         inoc = factor(inoc, levels = c('charles', 'ica', 'nahant', 'fucus')),
        inoc_map = case_when(
                         inoc == 'charles' ~ 'Brackish',
                         inoc == 'ica' ~ 'Estuary',
                         inoc == 'nahant' ~ 'Marine 1',
                         inoc == 'fucus' ~ 'Marine 2'
                       ),
        temp = ifelse(temp == '15C', '15C', '20C'))

# Minimum salinity
community_smin <- ggplot(community_gr_traits ) +
  geom_point(aes(x = salt_num, y = community_smin, color = inoc_map, group = interaction(salt_num, inoc_map), shape = temp), size = 3, show.legend = T) +
  facet_wrap(vars(inoc_map), ncol = 4) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  labs(x = 'Salinity (g/l)', y = 'Community Smin', color = 'Inoculum', shape = 'Temperature') +
  scale_color_manual(values = color_list) +
  theme(legend.position = 'bottom',
        text = element_text(size = 20),
        panel.spacing = unit(1.5, 'lines'))

community_smin
#ggsave(paste0(fig_dir, 'community_smin.png'), width = 8, height = 4)

# Maximum salinity
community_smax <- ggplot(community_gr_traits ) +
  geom_point(aes(x = salt_num, y = community_smax, color = inoc_map, group = interaction(salt_num, inoc_map), shape = temp), size = 3, show.legend = F) +
  facet_wrap(vars(inoc_map), ncol = 4) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = 'Salinity (g/l)', y = 'Community Smax', color = 'Inoculum', shape = 'Temperature') +
  scale_color_manual(values = color_list) +
  theme(legend.position = 'bottom',
        text = element_text(size = 20),
        panel.spacing = unit(1.5, 'lines'))

community_smax
#ggsave(paste0(fig_dir, 'community_smax.png'), width = 8, height = 4)

# Salinity optimum
community_sopt <- ggplot(community_gr_traits ) +
  geom_point(aes(x = salt_num, y = community_sopt, color = inoc_map, group = interaction(salt_num, inoc_map), shape = temp), size = 3, show.legend = F) +
  facet_wrap(vars(inoc_map), ncol = 4) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  labs(x = 'Salinity (g/l)', y = 'Community Sopt', color = 'Inoculum', shape = 'Temperature') +
  scale_color_manual(values = color_list) +
  theme(legend.position = 'bottom',
        text = element_text(size = 20),
        panel.spacing = unit(1.5, 'lines'))

community_sopt
#ggsave(paste0(fig_dir, 'community_sopt.png'), width = 8, height = 4)

# Max growth rate
community_maxgr <- ggplot(community_gr_traits ) +
  geom_point(aes(x = salt_num, y = community_maxgr, color = inoc_map, group = interaction(salt_num, inoc_map), shape = temp), size = 3, show.legend = F) +
  facet_wrap(vars(inoc_map), ncol = 4) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  labs(x = 'Salinity (g/l)', y = 'Community max growthrate', color = 'Inoculum', shape = 'Temperature') +
  scale_color_manual(values = color_list) +
  theme(legend.position = 'bottom',
        text = element_text(size = 20),
        panel.spacing = unit(1.5, 'lines'))

community_maxgr
#ggsave(paste0(fig_dir, 'community_maxgr.png'), width = 8, height = 4)

# Salinity range
community_range_halfmax <- ggplot(community_gr_traits ) +
  geom_point(aes(x = salt_num, y = community_range_halfmax, color = inoc_map, group = interaction(salt_num, inoc_map), shape = temp), size = 3, show.legend = F) +
  facet_wrap(vars(inoc_map), ncol = 4) +
  scale_x_continuous(breaks = c(16, 31, 46)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = 'Salinity (g/l)', y = 'Community \nrange halfmax', color = 'Inoculum', shape = 'Temperature') +
  scale_color_manual(values = color_list) +
  theme(legend.position = 'bottom',
        text = element_text(size = 20),
        panel.spacing = unit(1.5, 'lines'))

community_range_halfmax
#ggsave(paste0(fig_dir, 'community_range_halfmax.png'), width = 8, height = 4)

## Put all these panels together 

community_smin + community_smax +
community_sopt + community_range_halfmax +
plot_layout(nrow = 4, guides = 'collect') &
  theme(legend.position = 'bottom')

ggsave(paste0(fig_dir, 'FigS7_community_traits.png'), width = 11, height = 12)

