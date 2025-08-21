###########################################################
## Analyzing the data from Beaufort, Alaska (Kellogg et al)
## Author: JS Huisman
###########################################################

library(tidyverse)
library(readxl)
#library(dada2)
library(phyloseq)


main_data_dir <- "../sequencing/dilution_experiment_Illumina/Jana_Martina_data"

data_dir = '../sequencing/environmental_data/Kellogg_Beaufort'
proc_data_dir = paste0(data_dir, '/processed_data')

fig_data = '../figures/env_data'

theme_set(theme_minimal() + theme(text = element_text(size = 20)))

###########################################################
# Extra environmental metadata

#To join the two datasets together, please use the provided site codes (column "site_name" here) 
#and collection dates (column "collection_date" here) in each dataset. Note that the site codes in the sequence data 
#are without hyphens (e.g. JAA) while site codes in the environmental data have hyphens (e.g. JA-A).

# please cite the original GenBank data, journal article, or related Arctic Data Center dataset as appropriate. 

# In each case the first sheet has some metadata information
temp_salinity_data_hobo = read_xlsx(paste0(data_dir, '/metadata/Arctic-Lagoons_Hobo_2011-2014.xlsx'), sheet = 'NSF E Beaufort hobo')
temp_salinity_data_ysi = read_xlsx(paste0(data_dir, '/metadata/Arctic-Lagoons_YSI_2011-2013.xlsx'), sheet = 2)
water_sample_data = read_xlsx(paste0(data_dir, '/metadata/Arctic-Lagoons_discrete_water_samples_2011-2013.xlsx'), sheet = 2)
stable_isotope_data = read_xlsx(paste0(data_dir, '/metadata/Arctic-Lagoons_biota_stable_isotopes_2011-2013.xlsx'), sheet = 2)

##### Load all data ##########################################

## Load data - overall phyloseq object called ps ####
seqtab.nochim = readRDS(paste0(proc_data_dir, '/seqtab.nochim.rds'))
taxa = readRDS(paste0(proc_data_dir, '/taxa.rds'))
metadata = read_csv(paste0(proc_data_dir, '/metadata.csv')) 

samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(sample = samples.out)
samdf <- left_join(samdf, metadata, by = 'sample')
rownames(samdf) <- samples.out

# construct phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# make short names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

## Remove non-bacterial + mitochondrial ASVs ######################
ps <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      !is.na(Phylum)
  )

###### Check read distribution of samples #################

# sample_sum sums reads across all asvs observed in a given sample # thus read depth
sample_sum_df <- data.frame(sum = sample_sums(ps))
# add meta information on the samples
seq_depth_df <- as_tibble(sample_sum_df, rownames = 'sample') %>% 
  left_join(metadata, by = 'sample') 

seq_depth_df %>%
  filter(sum == 0) %>%
  View()

# Plot a histogram of read depths
ggplot(seq_depth_df, aes(x = sum)) +
  #facet_wrap(vars(Organism)) +
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  labs(x = "Read depth", y = 'Number of samples') 
ggsave(paste0(fig_data, '/read_depth_Kellogg.pdf'), width = 10)

## Calculate diversity  ###########################
#D0 = richness

D0 <- estimate_richness(ps,split=TRUE, measures="Observed")

diversity_df = data.frame(sample = row.names(D0), D0) %>%
  left_join(metadata, by = 'sample') %>%
  rename(diversity = Observed)

write_csv(diversity_df, paste0(proc_data_dir, '/diversity_df.csv'))

## Load rRNA database #############################
rrndb_NCBI <- read_delim('../processed_data/sequencing/rrnDB-5.8_pantaxa_stats_NCBI.tsv', delim = '\t') %>%
  select (rank, name, copy_number_NCBI = mean)
rrndb_names <- c("superkingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species')

## Calculate rel abundance & mean copynumber ###############

df_names <- c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

# turn ps object into a dataframe
df_all <- psmelt(ps) %>%
  mutate(Species = NA)

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
  mutate(weigh_abund = Abundance/copy_number_NCBI) %>%
  group_by(Sample) %>%
  mutate(rel_abund = weigh_abund / sum(weigh_abund, na.rm = T))

write_csv(rel_abund_df, paste0(proc_data_dir, '/rel_abund_df.csv'))

mcn_df <- rel_abund_df %>%
  group_by(sample_name, date, geo_loc_name, site_name, env_medium, Organism) %>%
  summarise(mcn = weighted.mean(copy_number_NCBI, rel_abund, na.rm = T) )
write_csv(mcn_df, paste0(proc_data_dir, '/mcn_df.csv'))

##########
#Initial plot

ggplot(mcn_df) +
  geom_point(aes(x = date, y = mcn, color = env_medium))


