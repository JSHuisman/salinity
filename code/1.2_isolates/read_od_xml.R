###########################################################
## Reading growth data from xml files
## Author: JS Huisman
##
## The script assumes 48hr growth data from a Tecan 
## Infinite M Nano
###########################################################

library(tidyverse)
#library(ggplot2)
library(readxl)
library(xml2)
library(openxlsx)
#library(xmltools)

raw_data_dir = '../raw_data/isolates/growth_measurements/od_files_from_platereader/'
write_data_dir = '../raw_data/isolates/growth_measurements/od_48hr_growth/'
#fig_dir = '../figures/growth_rates/'

###########################################################

# List all files in the raw data directory - currently selecting data from 2023-2025!
od_files = grep('^2025.*.xml', list.files(raw_data_dir), value = T)
dates = str_extract(od_files, '^202[3-5]-[0-9][0-9]-[0-9][0-9]')
plates = gsub('.xml', '',str_extract(od_files, '[0-9].xml$'))

for (i in 1:length(od_files)){
# read in OD data 
file = paste0(raw_data_dir, od_files[i])
xml_data <- suppressWarnings(read_xml(file))

## Structure of the xml: 
# At the highest level: [1] Plate name, [2] Software info, [3] Script info, [4] Measurement
    #xml_children(xml_data)
# Below [4]: Time, Parameters, as many data cycles as timepoints (960)
    #xml_children(xml_children(xml_data)[4])
# Below [3 - 963]: measurements for 96 wells at that timepoint

# read out cycles
if (length(xml_children(xml_data)) >3){
  cycles = as.numeric(xml_attr(xml_children(xml_children(xml_data)[4]), 'Cycle')[3:962])  
} else{
  next
}

# Some xml files are empty due to a plate reader crash - skip these
if (any(is.na(cycles))){
  next
}

# read out temperature during the experiment
temperatures = as.numeric(xml_attr(xml_children(xml_children(xml_data)[4]), 'Temperature')[3:962])

# Iterate through xml file and save od values
exp_info_df = tibble('cycle' = cycles, 'time' = (cycles-1)*180, 'temp' = temperatures, 'time_hr' = ((cycles-1)*180)/3600)
od_df = data.frame()
for (cycle in cycles){
  wells = xml_attr(xml_children(xml_children(xml_children(xml_data)[4])[cycle+2]), 'Pos')
  od_values = xml_double(xml_children(xml_children(xml_children(xml_data)[4])[cycle+2]))
  new_df = tibble('cycle' = cycle, 'well' = wells, 'od' = od_values)
  
  od_df = bind_rows(od_df, new_df)
}

## Add some generally useful additional columns about the experiment
salt_values = c(0, 5, 15, 20, 30, 35, 40, 45, 50, 65, 80, 100)

final_od_df <- as_tibble(od_df) %>%
  left_join(exp_info_df, by = 'cycle') %>%
  mutate(well_row = gsub('[0-9]*$', '', well),
         well_col = factor(gsub('[A-H]', '', well), levels = 1:12),
         well_salt = salt_values[as.numeric(gsub('[A-H]', '', well))])
#final_od_df

# Write OD values to file
write_csv(final_od_df, paste0(write_data_dir, 'od_', dates[i], '_plate_', plates[i], '.csv'))
write.xlsx(final_od_df, paste0(write_data_dir, '48h_growthrate_', dates[i], '_plate_', plates[i], '.xlsx'))

# plotting the od_data
# ggplot(final_od_df) +
#   geom_line(aes(x = time_hr, y = od)) +
#   facet_grid(rows = vars(well_row), cols = vars(well_salt))
# 
# ggsave(paste0(fig_dir, 'all_growth_curves_', dates[i], '_plate', plates[i], '.png'), width = 16, height = 9)

}



