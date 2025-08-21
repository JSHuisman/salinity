######################################################################################################################
## Pre-process competition data from Oct 2024
## Author: JS Huisman
######################################################################################################################

library(tidyverse)
library(readxl)
library(lubridate)

data_dir = '../raw_data/competition/2024_Oct'
proc_data_dir = '../data/3_pairwise_competition'

######################################################################################################################

unique_pairs = c('40-2', '40-8', '40-39', '40-II.F2',
                 '6-8', '6-25', '2-25', '2-32', '2-II.F2')
faster_pair = c("6-8" = 8, "6-25" = 6,
                '2-25' = 2, '2-32' = 32, "2-II.F2" = 2,
                "40-2" = 2, "40-8" = 8, '40-39' = 40, "40-II.F2" = 40)

###### Combine files from the different plating dates ################################################################

# Read in all plate counts
D0_counts = read_xlsx(paste0(data_dir, '/plating_2024-10-25.xlsx'), na = c('', 'TMTC'))%>%
  mutate(plate_rep = as.character(plate_rep),
         volume_plated_ul = '10')
D1_counts = read_xlsx(paste0(data_dir, '/plating_2024-10-27.xlsx'), na = c('', 'TMTC')) %>%
  mutate(plate_rep = as.character(plate_rep),
         volume_plated_ul = '10')
D3_counts = read_xlsx(paste0(data_dir, '/plating_2024-10-31.xlsx'), na = c('', 'TMTC')) %>%
  mutate(plate_rep = as.character(plate_rep),
         volume_plated_ul = '10')
D5_counts = read_xlsx(paste0(data_dir, '/plating_2024-11-04.xlsx'), na = c('', 'TMTC')) %>%
  mutate(plate_rep = as.character(plate_rep),
         volume_plated_ul = '10')
D7_counts = read_xlsx(paste0(data_dir, '/plating_2024-11-08.xlsx'), na = c('', 'TMTC', '-')) %>%
  mutate(plate_rep = as.character(plate_rep))
#D7_counts_extra = read_xlsx(paste0(data_dir, '/plating_extra_2024-11-08.xlsx'), na = c('', 'TMTC', '-'))

file_mapping = c('0' = 'C0', '1' = 'C1', '3' = 'C3', '5' = 'C5', '7' = 'C7')

# Duplicate D0 so all salinities are coupled to a copy of the D0 counts
D0_counts_extra <- crossing(D0_counts, salt = c('16', '31', '46', '61')) %>%
  select(-salinity) %>%
  rename(salinity = salt)

#D7_counts_extra
count_data = bind_rows(lst(D0_counts_extra, D1_counts, D3_counts, D5_counts, D7_counts), .id = 'count_file') %>%
  mutate(day_ID = gsub('_counts|_counts_extra', '', count_file))

###### Pick dilution and add error bars ################################################################

# Consolidate the plate counts:
# remove wells that were visibly contaminated
# sum across multiple plates counted at the same dilution (on the same day / of the same well)
    # D0 had 3 plates; in most other cases it will be just 1 [these are tech. replicates]
# On day 7; 2-25 and 6-25 should use dilution d7, all other pairs d6 
consolidated_count_data <- count_data %>%
  select(-count_file, -row, -col) %>%
  filter(is.na(contam),
         pair != '6-8') %>%
  group_by(date_of_plating, day_ID, dilution, salinity, 
           well, competition, pair, bio_rep) %>%
  summarize(total_count_1 = sum(competitor_1, na.rm = T),
            total_count_2 = sum(competitor_2, na.rm = T),
            total_volume = sum(as.numeric(volume_plated_ul), na.rm = T),
            .groups = 'drop') %>%
  filter(day_ID == 'D0' | 
           (pair %in% c('2-25', '6-25') & dilution == 'd7') | 
           ((!pair %in% c('2-25', '6-25')) & dilution == 'd6') )


# Transform counts into CFU/mL and compute ratios
combined_data = consolidated_count_data %>%
  mutate(date_of_plating = as_date(date_of_plating),
         dil_factor = 10^as.numeric(gsub('d', '', dilution))) %>%
  group_by(day_ID, salinity, competition, pair) %>%
  mutate(rep = as.factor(row_number()),
         across(starts_with('n_'), function(x){ifelse(is.na(x), 0, x)} )) %>% #this adds a replicate number
  rowwise() %>%
  mutate(salt = as.factor(salinity),
         total_count = sum(total_count_1, total_count_2, na.rm = T),
         CFU_1 = (total_count_1 * dil_factor*1000)/as.numeric(total_volume),
         sd_CFU_1 = (sqrt(total_count_1) * dil_factor*1000)/as.numeric(total_volume), # factor 1000 because the volume is in uL; CFU is CFU/mL
         CFU_2 = (total_count_2 * dil_factor*1000)/as.numeric(total_volume),
         sd_CFU_2 = (sqrt(total_count_2) * dil_factor*1000)/as.numeric(total_volume),
         CFU_total = (total_count * dil_factor*1000)/as.numeric(total_volume),
         sd_CFU_total = (sqrt(total_count) * dil_factor*1000)/as.numeric(total_volume),
         ratio_strain1 = total_count_1/total_count,
         ratio_strain2 = total_count_2/total_count,
         strain1 = ifelse(competition!=100, str_split_1(pair, '-')[1], pair),
         strain2 = ifelse(competition!=100, str_split_1(pair, '-')[2], pair),
         ratio_fast = case_when(
                      competition == 100 ~ ratio_strain1,
                      strain1 == faster_pair[pair] ~ ratio_strain1,
                      strain1 != faster_pair[pair] ~ ratio_strain2)) %>%
  ungroup()

write_csv(combined_data, file = paste0(proc_data_dir, '/competition_plating_data.csv'))


##### Different ways of grouping/summarizing the data ###############################################

# summarize across biological replicates but not inoculation ratios
rep_summarized_data <- combined_data %>%
  group_by(day_ID, salinity, salt, competition, pair, strain1, strain2) %>%
  summarize(mean_CFU_1 = mean(CFU_1, na.rm = T),
            mean_CFU_2 = mean(CFU_2, na.rm = T),
            mean_ratio_fast = mean(ratio_fast, na.rm = T))

write_csv(rep_summarized_data, file = paste0(proc_data_dir, '/rep_summarized_data.csv'))


# summarize across all biological replicates and inoculation ratios
summarized_data <- combined_data %>%
  mutate(pair_type = ifelse(competition == '100', 'mono', 'pair')) %>%
  group_by(day_ID, salinity, salt, pair, pair_type, strain1, strain2) %>%
  summarize(mean_CFU_1 = mean(CFU_1, na.rm = T),
            mean_CFU_2 = mean(CFU_2, na.rm = T),
            mean_ratio_fast = mean(ratio_fast, na.rm = T))

write_csv(summarized_data, file = paste0(proc_data_dir, '/summarized_data.csv'))


##### Read OD data  ###############################################

init_od_df = tibble(filename = '2024-10-25_initOD', date_of_plating = as_date('2024-10-25'),
                    day_ID = 'D0', salinity = c('16', '31', '46', '61'))

days = c('2024-10-27_D1_','2024-10-29_D2_','2024-10-31_D3_',
         '2024-11-02_D4_','2024-11-04_D5_','2024-11-06_D6_','2024-11-08_D7_')
cond = c('16gl', '31gl', '46gl', '61gl')

od_files = crossing(days, cond) %>%
  rowwise() %>%
  mutate(filename = paste0(days, cond),
         date_of_plating = as_date(str_split_1(days, '_')[1]),
         day_ID = str_split_1(days, '_')[2],
         salinity = gsub('gl', '', cond)) %>%
  select(-days, -cond) %>%
  bind_rows(init_od_df) %>%
  arrange(day_ID, salinity)

# Loop through all measurement dates & individual plates
OD_results_df = tibble()
for (exp_id in 1:dim(od_files)[1]){

  raw_od_data = read_xlsx(paste0(data_dir,'/', od_files[exp_id, 'filename'], '.xlsx'), skip = 25, n_max = 8)

  long_od_data <- raw_od_data %>%
    pivot_longer(cols = -1, names_to = 'col') %>%
    rename('row' = `<>`, 'od' = value) %>%
    mutate(well = paste0(row, col),
           col = as.numeric(col)) %>%
    bind_cols(od_files[exp_id,])

  OD_results_df = bind_rows(OD_results_df, long_od_data)
}

plate_map_long <- read_xlsx(paste0(data_dir,'/plate_map_long.xlsx'))

contam_info <- count_data %>%
  #filter(!is.na(contam)) %>%
  select(date_of_plating, well, contam, salinity, day_ID) %>%
  mutate(date_of_plating = as_date(date_of_plating))


od_meta_df <- OD_results_df %>%
  left_join(plate_map_long, by = c('row', 'col', 'well')) %>%
  left_join(contam_info, by = c('well', 'date_of_plating', 'salinity', 'day_ID')) %>%
  filter(pair != '6-8',
         is.na(contam)) %>%
  select(-contam)

write.csv(od_meta_df, file = paste0(proc_data_dir, '/od_data.csv'))

