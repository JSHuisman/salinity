######################################################################################################################
## Estimate growth rates from all OD data
## Author: JS Huisman
######################################################################################################################
library(tidyverse)
library(ggplot2)
library(readxl)
library(gcplyr)
library(patchwork)

raw_data_dir = '../raw_data/isolates/'
raw_growthdata_dir = '../raw_data/isolates/growth_measurements/od_48hr_growth/'
proc_data_dir = '../processed_data/isolates/growth_data/'
fig_dir = '../figures/growth_rates/growth_per_plate/'

theme_set(theme_bw() + 
            theme(text = element_text(size = 20)))

######################################################################################################################
## Loading metadata ####

# Isolate information
isolate_info <- read_csv(paste0('../processed_data/isolates/paper_isolate_ids.csv')) %>%
  mutate(name = paste0(genus, ' ', species))


# Growth rate measurement information
meas_info <- read_xlsx(paste0(raw_data_dir, 'growth_measurements/metadata_growth_measurement.xlsx'), na = c("", "-"))
meas_info <- meas_info %>%
  mutate(isolate_id = ifelse(isolate_id == 'i59', 'isolate_59', isolate_id)) %>%
  left_join(isolate_info, by = 'isolate_id') %>%
  filter(as.logical(include_plate),
         !is.na(isolate_id)) %>%
  mutate(measurement_date = as_date(measurement_date))

meas_dates = unique(meas_info$measurement_date)

unique_experiments = meas_info %>%
  select(measurement_date, plate_id) %>%
  distinct() 

# salinities measured (kept constant across all experiments)
salt_values = c(0, 5, 15, 20, 30, 35, 40, 45, 50, 65, 80, 100)
col_names <- setNames(salt_values, factor(1:12))

######################################################################################################################
## Read all OD files, plot raw OD plate-wise, compute growth rates from OD data ####
# Computes gcplyr growth rates - with background/blank subtraction

# Flag to run all analyses if the analysis pipeline changed
redo_all = T

# Loop through all measurement dates & individual plates
all_results_df = tibble()
for (exp_id in 1:dim(unique_experiments)[1]){
  
  exp_date = unique_experiments$measurement_date[exp_id]
  
  ## Exclude wells based on experiment info 
  # parse which wells are listed to exclude; remove NA entries
  exclusion_col = meas_info %>% filter(measurement_date == exp_date, 
                                       plate_id == unique_experiments$plate_id[exp_id]) %>% pull(exclude_wells)
  exclude_wells = unlist(str_split(exclusion_col[!is.na(exclusion_col)], ', '))
  
  # Load OD data
  if (is.na(unique_experiments$plate_id[exp_id])){
    plate = 1
    filename = paste0('48h_growthrate_', exp_date, '.xlsx')    
  } else {
    plate = unique_experiments$plate_id[exp_id]
    filename = paste0('48h_growthrate_', exp_date, '_plate_', plate, '.xlsx')
  }
  
  # run the analysis only if there is new data and/or we want to rerun the analysis
  if(file.exists(paste0(proc_data_dir, 'growthdata_', exp_date, '_plate_', plate, '.xlsx')) & !redo_all){
    results_df <- read_csv(paste0(proc_data_dir, 'growthdata_', exp_date, '_plate_', plate, '.xlsx'))
    results_df <- results_df %>% mutate(isolate_id = factor(isolate_id))
    all_results_df = bind_rows(all_results_df, results_df)
    next
  }
  
  # read in the raw OD data - there are two types of files, before/after the Tecan switched to xml output
  if (is.na(unique_experiments$plate_id[exp_id])){
    raw_od_data <- try(read_xlsx(paste0(raw_growthdata_dir, filename), skip = 33, n_max = 100))
    if ('try-error' %in% class(raw_od_data)){
      next
    }
    
    first_od_data <- raw_od_data %>%
      pivot_longer(cols = -`Cycle Nr.`, names_to = 'cycle') %>%
      pivot_wider(names_from = `Cycle Nr.`) %>%
      rename('time' = `Time [s]`, 'temp' = `Temp. [°C]`) %>%
      pivot_longer(cols = -c('cycle', 'time', 'temp'), names_to = 'well', values_to = 'od') %>%
      mutate(well_row = gsub('[0-9]*$', '', well),
             well_col = factor(gsub('[A-H]', '', well), levels = 1:12),
             well_salt = salt_values[as.numeric(gsub('[A-H]', '', well))],
             time_hr = time/3600,
             cycle = as.numeric(cycle))
  } else {
    raw_od_data <- try(read_xlsx(paste0(raw_growthdata_dir, filename)))
    if ('try-error' %in% class(raw_od_data)){
      next
    }
    first_od_data <- raw_od_data
  }

  
  ## Combine with measurement info
  first_od_data <- first_od_data  %>%
    mutate(measurement_date = as_date(exp_date),
           plate_id = unique_experiments$plate_id[exp_id]) %>%
  left_join(meas_info %>% dplyr::rename(well_row = plate_row), 
            by = c('measurement_date', 'plate_id', 'well_row')) %>% 
    filter(isolate_id != '-') %>%
    mutate(plate_id = ifelse(is.na(plate_id), 1, plate_id),
           isolate_id = factor(isolate_id)) %>%
    group_by(well) %>%
    arrange(time_hr) %>%
    mutate(smooth_od = zoo::rollapply(od, 7, mean, align = "center", na.rm = T, fill = NA) )%>%
    ungroup()
  
  # gcplyr derivatives (as background subtraction we take the minimal OD value measured)
  od_analysis_df <- first_od_data %>%
    filter(time_hr >= 3,
           !well %in% exclude_wells) %>% # with our protocol, we do not expect maximal growth within 3 hours
    group_by(well) %>%
    mutate(deriv = gcplyr::calc_deriv(x = time_hr, y = smooth_od, trans_y = 'log',
                                      percapita = TRUE, blank = min(first_od_data$smooth_od, na.rm = T),
                                      window_width_n = 21),
           doub_time = doubling_time(y = deriv))
  
  # put all analyses together
  results_df <- od_analysis_df %>%
    group_by(well, well_row, well_salt, isolate_id) %>%
    summarize(min_dens = first_minima(y = smooth_od, return = 'y'),
              max_dens = max_gc(smooth_od, na.rm = TRUE), #Maximum OD attained
              dens_diff = max_dens - min_dens,
              #first_maxima_y = first_maxima(x = time_hr, y = smooth_od, 
              #                              return = "y"),
              #first_max_time = extr_val(time_hr, which(smooth_od == first_maxima_y)[1]),
              lag_time = lag_time(x = time_hr, y = smooth_od, deriv = deriv, y0 = min_dens,
                                  blank = min(first_od_data$smooth_od, na.rm = T), na.rm = T),
              max_percap = max_gc(deriv, na.rm = TRUE),
              max_percap_time = extr_val(time_hr, which_max_gc(deriv)) 
              ) %>%
    mutate(measurement_date = as_date(exp_date),
           plate_id = plate) %>%
    mutate(max_percap = ifelse(max_dens <= 0.15 | dens_diff < 0.05, 0, max_percap), # remove cases that are not significantly different from 0
           max_percap_time = ifelse(max_dens <= 0.15 | dens_diff < 0.05, NA, max_percap_time),
           lag_time = ifelse(max_dens <= 0.15 | dens_diff < 0.05, NA, lag_time)) 
  
  # write out the processed data
  write_csv(results_df, paste0(proc_data_dir, 'growthdata_', exp_date, '_plate_', plate, '.xlsx'))
  all_results_df = bind_rows(all_results_df, results_df)
  
  ## PLOTTING ####
  # plotting the od_data - log plot
  ggplot(first_od_data) +
    geom_line(aes(x = time_hr, y = smooth_od)) +
    geom_vline(data = results_df, aes(xintercept= max_percap_time), linetype = 'dashed', color = "red") +
    coord_trans(y = 'log') +
    labs(x = 'Time (hr)', y = 'Log OD600') +
    scale_x_continuous(breaks = c(0, 24, 48)) +
    facet_grid(rows = vars(isolate_id), cols = vars(well_salt)) +
    theme(panel.spacing = unit(1, 'lines')) 
  ggsave(paste0(fig_dir, exp_date, '_plate_', plate, '_log_growth.png'), width = 16, height = 9)
  
  ## Plot gcplyr derivatives and max growth rate
  ggplot(od_analysis_df) +
    geom_line(aes(x = time_hr, y = deriv)) +
    geom_point(data = results_df, 
               aes(x = max_percap_time, y = max_percap),
               size = 2, color = "red") +
    labs(x = 'Time (hr)', y = 'Derivative of Log OD600') +
    scale_x_continuous(breaks = c(0, 24, 48)) +
    facet_grid(rows = vars(isolate_id), cols = vars(well_salt), scales = 'free_y') +
    theme(panel.spacing = unit(1, 'lines')) 
  ggsave(paste0(fig_dir, exp_date, '_plate_', plate, '_gcplyr_derivs.png'), width = 16, height = 9)
  
}

## add measurement info
all_results_df <- all_results_df %>%
  left_join(meas_info %>% dplyr::rename(well_row = plate_row) %>% 
              mutate(plate_id = ifelse(is.na(plate_id), 1, plate_id),
                     isolate_id = factor(isolate_id)), 
            by = c('measurement_date', 'plate_id', 'well_row', 'isolate_id')) %>% 
  mutate(exp_id = paste0(measurement_date, '_', plate_id)) %>%
  mutate(replicate = factor(exp_id)) 

write_csv(all_results_df, paste0(proc_data_dir, 'all_growthdata.csv'))

###########################################################
# write out summary values for each isolate

average_growth_df = all_results_df %>% 
  select(measurement_date, plate_id, isolate_id, well_salt, max_dens, lag_time, max_percap, class, order, family, genus, species) %>%
  group_by(isolate_id, well_salt, class, order, family, genus, species) %>%
  summarize(lag_time = mean(lag_time), 
            growthrate = mean(max_percap),
            gr_sd = sd(max_percap),
            carr_cap = mean(max_dens),
            carr_cap_sd = sd(max_dens),
            n_growth_meas = n()
  )

write_csv(average_growth_df, paste0(proc_data_dir, 'average_growthrates.csv'))

###########################################################
# update isolate file with information on # growth measurements

isolate_info <- read_csv(paste0('../processed_data/isolates/', 'paper_isolate_ids.csv'))

isolate_info_w_ngrowth <- all_results_df %>% 
  group_by(isolate_id) %>%
  summarize(n_growth_meas = floor(n() / 12) ) %>%
  right_join(isolate_info, by = "isolate_id") %>%
  mutate(n_growth_meas = ifelse(is.na(n_growth_meas), 0, n_growth_meas)) %>%
  arrange(isolate_id)

write_csv(isolate_info_w_ngrowth, paste0(proc_data_dir, 'paper_isolate_ids_ngrowth.csv'))

