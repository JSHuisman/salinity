###########################################################
## Loading OD and pH data for the serial dilution experiment
## Author: JS Huisman
###########################################################

library(tidyverse)
library(ggplot2)
library(readxl)

data_dir = '../raw_data/daily_dilutions/'

###########################################################

## specifying which wells correspond to which conditions
raw_plate_map <- read_xlsx(paste0(data_dir, 'plate_map.xlsx'),
                           skip = 1, col_names = c('row_id', paste0('col_', 1:12)))

plate_map <- raw_plate_map %>%
  pivot_longer(cols = 2:13, names_to = 'col_id',values_to = "cond") %>%
  mutate(salt = factor(case_when(
    col_id %in% paste0('col_', 1:4) ~ 16,
    col_id %in% paste0('col_', 5:8) ~ 31,
    col_id %in% paste0('col_', 9:12) ~ 46
  )))

## make list of OD files
od_files <- list.files(path = paste0(data_dir, 'od_data/'), pattern = '.xlsx')

## load OD data
all_od_data <- data.frame()

for (filename in od_files){
  raw_od_data <- read_xlsx(paste0(data_dir, 'od_data/', filename),
                           range = 'A25:M32', col_names = c('row_id', paste0('col_', 1:12)))
  
  # deduce day and temp from filename
  file_temp = str_split_1(filename, '_')[3]
  file_day = sub('.xlsx', '', str_split_1(filename, '_')[4])
  
  # match conditions to od values using the plate_map
  od_data <- raw_od_data %>% 
    pivot_longer(cols = 2:13, names_to = 'col_id', values_to = "od") %>%
    right_join(plate_map, by = c('row_id', 'col_id')) %>%
    mutate(temp = file_temp,
           day = file_day) %>%
    select(day, temp, salt, cond, od)
  
  all_od_data <- bind_rows(all_od_data, od_data)
}

write_csv(all_od_data, '../processed_data/daily_dilutions/all_od_data.csv')



######################################################################################################################
## Loading pH data

## specifying which wells correspond to which conditions
raw_plate_map <- read_xlsx(paste0(data_dir, 'plate_map.xlsx'),
                           skip = 1, col_names = c('row_id', paste0('col_', 1:12)))

plate_map <- raw_plate_map %>%
  pivot_longer(cols = 2:13, names_to = 'col_id',values_to = "cond") %>%
  mutate(salt = factor(case_when(
    col_id %in% paste0('col_', 1:4) ~ 16,
    col_id %in% paste0('col_', 5:8) ~ 31,
    col_id %in% paste0('col_', 9:12) ~ 46
  )))

## make list of pH files
ph_files <- list.files(path = paste0(data_dir, 'ph_data/'), pattern = '.xlsx')

## load pH data
all_ph_data <- data.frame()

for (filename in ph_files){
  raw_ph_data <- read_xlsx(paste0(data_dir, 'ph_data/', filename), na = "blank",
                           skip = 1, col_names = c('row_id', paste0('col_', 1:12)))
  
  # deduce day and temp from filename
  file_temp = str_split_1(filename, '_')[3]
  file_day = str_split_1(filename, '_')[4]
  
  # match conditions to ph values using the plate_map
  ph_data <- raw_ph_data %>% 
    pivot_longer(cols = 2:13, names_to = 'col_id', values_to = "ph") %>%
    right_join(plate_map, by = c('row_id', 'col_id')) %>%
    mutate(temp = file_temp,
           day = file_day) %>%
    select(day, temp, salt, cond, ph)
  
  all_ph_data <- bind_rows(all_ph_data, ph_data)
}



write_csv(all_ph_data, '../processed_data/daily_dilutions/all_ph_data.csv')













