###########################################################
## Analyzing the data from the Baltic sea
## Author: JS Huisman
###########################################################

library(tidyverse)
library(readxl)
#library(dada2)
library(phyloseq)

main_data_dir <- "../sequencing/dilution_experiment_Illumina/Jana_Martina_data"

data_dir = '../sequencing/environmental_data/Latz_Baltic'
proc_data_dir = paste0(data_dir, '/processed_data')

fig_dir = '../figures/env_data'


seqtab <- read.csv(paste0(data_dir, '/processed_data/seqtab_16S.txt'), header = TRUE, row.names = 1, sep = '\t')
microbe_abund <- read.csv(paste0(data_dir, '/processed_data/seqtab_16S.txt'), header = TRUE, row.names = NULL, sep = '\t') %>%
  rename(ASV = row.names) %>%
  pivot_longer(starts_with('P'), names_to = 'sample', values_to = 'abundance')

taxa <- read.csv(paste0(data_dir, '/processed_data/taxa_16S_SILVA.txt'), row.names = NULL, sep = ' ')%>%
  rename(ASV = row.names)
metadata <- read.csv(paste0(data_dir, '/processed_data/metadata.txt'), sep = '\t') %>%
  rename(sample = sample_alias.library_name)

theme_set(theme_minimal() + theme(text = element_text(size = 20)))
###########################################################

#Steps:
# use taxa_df with the rRNA db to assign a copy number to each ASV
# use microbe abund with the copy number-taxa, to compute MCN per sample
# use sample df to link to environmental parameters

##### Arrange data into phyloseq ##########################################

## Load data - overall phyloseq object called ps ####
samples.out <- rownames(t(seqtab))
samdf <- data.frame(sample = samples.out)
samdf <- left_join(samdf, metadata, by = 'sample')
rownames(samdf) <- samples.out

## THIS DOESN"T WORK AT THE MOMENT
# # construct phyloseq object
# ps <- phyloseq(otu_table(seqtab, taxa_are_rows=TRUE), 
#                sample_data(samdf), 
#                tax_table(taxa))
# 
# # make short names
# dna <- Biostrings::DNAStringSet(taxa_names(ps))
# names(dna) <- taxa_names(ps)
# ps <- merge_phyloseq(ps, dna)
# taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

## Load rRNA database #############################
rrndb_NCBI <- read_delim('../processed_data/sequencing/rrnDB-5.8_pantaxa_stats_NCBI.tsv', delim = '\t') %>%
  select (rank, name, copy_number_NCBI = mean)
rrndb_names <- c("superkingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species')

## Calculate rel abundance & mean copynumber ###############

df_names <- c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

# turn ps object into a dataframe
df_all <- microbe_abund %>%
  left_join(taxa, by = 'ASV') %>%
  left_join(metadata, by = 'sample') 

## we loop through the taxon levels from most to least specific
df_rrndb <- df_all
df_rrndb <- df_rrndb %>%
  dplyr::rename(name = df_names[7] ) %>%
  left_join(rrndb_NCBI %>% filter(rank == rrndb_names[7]), by = 'name') %>%
  dplyr::rename(!!df_names[7] := name) # this just renames the column back to the taxon level

for(level in 6:4){
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
  mutate(weigh_abund = abundance/copy_number_NCBI) %>%
  group_by(sample) %>%
  mutate(rel_abund = weigh_abund / sum(weigh_abund, na.rm = T))

write_csv(rel_abund_df, paste0(proc_data_dir, '/rel_abund_df.csv'))

# make mcn
mcn_df <- rel_abund_df %>%
  group_by(sample) %>%
  summarise(mcn = weighted.mean(copy_number_NCBI, rel_abund, na.rm = T) ) %>%
  left_join(metadata, by = 'sample')

write_csv(mcn_df, paste0(proc_data_dir, '/mcn_df.csv'))

##########
#Initial plot
#mcn_df <- read_csv(paste0(proc_data_dir, '/mcn_df.csv'))

ggplot(mcn_df) +
  geom_point(aes(x = as.numeric(salinity_average), y = mcn, color = as.numeric(temperature_water_CTD_0.1m))) +
  geom_smooth(aes(x = as.numeric(salinity_average), y = mcn), method='lm', color = 'black', alpha = 0.2) +
  labs(color = 'Temperature', y = 'Mean copynumber (MCN)', x = 'Salinity (g/L)')

ggsave(paste0(fig_dir, '/Latz_mcn_prelim.png'), width = 8, height = 5)




