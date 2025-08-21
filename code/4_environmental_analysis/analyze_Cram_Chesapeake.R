###########################################################
## Analyzing the data from the Chesapeake bay
## Author: JS Huisman
###########################################################

library(tidyverse)
library(readxl)
library(dada2)

main_data_dir <- "../sequencing/dilution_experiment_Illumina/Jana_Martina_data"

data_dir = '../sequencing/environmental_data/Cram_Chesapeake'
proc_data_dir = paste0(data_dir, '/processed_data')

fig_dir = '../figures/env_data'


microbe_abund <- read.csv('../sequencing/environmental_data/Cram_Chesapeake/emi16557-sup-0002-tables1.csv')
taxa_df <- read.csv('../sequencing/environmental_data/Cram_Chesapeake/emi16557-sup-0003-tables2.csv')
sample_df <- read.csv('../sequencing/environmental_data/Cram_Chesapeake/emi16557-sup-0004-tables3.csv')

theme_set(theme_minimal() + theme(text = element_text(size = 20)))
###########################################################

#Steps:
# use taxa_df with the rRNA db to assign a copy number to each ASV
# use microbe abund with the copy number-taxa, to compute MCN per sample
# use sample df to link to environmental parameters


## Remove non-bacterial + mitochondrial ASVs ######################

asv_to_keep <- taxa_df %>%
  filter(Kingdom == 'Bacteria',
      Family  != "Mitochondria",
      !is.na(Phylum)
  ) %>% pull(ASV)

subset_microbe_abund <- microbe_abund %>%
  filter(ASV %in% asv_to_keep)

## Load rRNA database #############################
rrndb_NCBI <- read_delim('../processed_data/sequencing/rrnDB-5.8_pantaxa_stats_NCBI.tsv', delim = '\t') %>%
  select (rank, name, copy_number_NCBI = mean)
rrndb_names <- c("superkingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species')

## Calculate rel abundance & mean copynumber ###############

df_names <- c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

# turn ps object into a dataframe
df_all <- subset_microbe_abund %>%
  left_join(taxa_df, by = 'ASV') %>%
  left_join(sample_df, by = 'ID') %>%
  mutate(Species = NA)

## we loop through the taxon levels from most to least specific
df_rrndb <- df_all
df_rrndb <- df_rrndb %>%
  dplyr::rename(name = df_names[7] ) %>%
  left_join(rrndb_NCBI %>% filter(rank == rrndb_names[7]), by = 'name') %>%
  dplyr::rename(!!df_names[7] := name) # this just renames the column back to the taxon level

for(level in 6:1){
  small_rrndb <- rrndb_NCBI %>% 
    filter(rank == rrndb_names[level]) 
  
  df_rrndb <- df_rrndb %>%
    dplyr::rename(name = df_names[level] ) %>%
    left_join(small_rrndb, by = 'name') %>%
    mutate(copy_number_NCBI.x = coalesce(copy_number_NCBI.x,copy_number_NCBI.y)) %>%
    mutate(rank.x = coalesce(rank.x,rank.y)) %>%
    select (-copy_number_NCBI.y, -rank.y) %>%
    dplyr::rename(copy_number_NCBI=copy_number_NCBI.x, rank =rank.x, !!df_names[level] := name)
}


# check that we didn't lose anything
dim(df_rrndb)[1] == dim(df_all)[1]
dim(df_rrndb %>% filter(is.na(copy_number_NCBI)))[1]

rel_abund_df <- df_rrndb %>%
  mutate(weigh_abund = reads/copy_number_NCBI) %>%
  group_by(ID) %>%
  mutate(rel_abund = weigh_abund / sum(weigh_abund, na.rm = T))

write_csv(rel_abund_df, paste0(proc_data_dir, '/rel_abund_df.csv'))

mcn_df <- rel_abund_df %>%
  group_by(ID,Station,Depth, CarbonPerLiter_mg, NitrogenPerLiter_mg, Temperature, salt = Salinity) %>%
  summarise(mcn = weighted.mean(copy_number_NCBI, rel_abund, na.rm = T) )

write_csv(mcn_df, paste0(proc_data_dir, '/mcn_df.csv'))

##########
#Initial plot
#mcn_df <- read_csv(paste0(proc_data_dir, '/mcn_df.csv'))

ggplot(mcn_df) +
  geom_point(aes(x = salt, y = mcn, color = ID), show.legend = F) +
  geom_smooth(aes(x = salt, y = mcn), method='lm', color = 'black', alpha = 0.2) +
  labs(color = 'Station', y = 'Mean copynumber (MCN)', x = 'Salinity (g/L)')

ggsave(paste0(fig_dir, '/Cram_mcn_prelim.png'), width = 8, height = 5)

