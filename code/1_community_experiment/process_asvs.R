###########################################################
# Analyzing ASVs
# 
# Author: JS Huisman
###########################################################

library(tidyverse)
#library(ggplot2)
library(phyloseq)
library(stringdist) # for the one base away analysis

proc_data_dir = paste0('../processed_data/sequencing')

fig_data = '../figures/sequencing_16S'


## LOAD DATA ##########################################################################################################
## Load data - overall phyloseq object called ps ####
seqtab.nochim = readRDS(paste0(proc_data_dir, '/seqtab.nochim.rds'))
taxa = readRDS(paste0(proc_data_dir, '/taxa.rds'))
metadata = read_csv(paste0(proc_data_dir, '/metadata.csv'))

samples.out <- rownames(seqtab.nochim)
plate <- sapply(strsplit(samples.out, "_"), `[`, 1)
well <- sapply(strsplit(samples.out, "_"), `[`, 2)
samdf <- data.frame(plate=plate, well=well)
samdf <- left_join(samdf, metadata, by = c('plate', 'well'))
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

# test for missing values - if this returns character(0) all is well
test <- metadata %>%
  mutate(samplename = paste0(plate, '_', well))
setdiff(test$samplename, rownames(seqtab.nochim))

## Load rRNA database #############################
rrndb_NCBI <- read_delim(paste0(proc_data_dir, "/rrnDB-5.8_pantaxa_stats_NCBI.tsv"), delim = '\t') %>%
  select (rank, name, copy_number_NCBI = mean)
rrndb_names <- c("superkingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species')


## General Preprocessing ##############################################################################################
## Remove weird ASVs ######################
ps <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      !is.na(Phylum)
  )



## Specific Preprocessing #############################################################################################
## ps_rarefied: rarefy ps to the same sequencing depth #########################
# rarefaction is not actually recommended!
ps_rarefied <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),
                                 rngseed = 42, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
saveRDS(ps_rarefied, paste0(proc_data_dir, '/ps_rarefied.rds'))

## ps_noRare: remove single hit ASVs from ps ###################################
#prevalence, here defined as the number of samples in which a taxon appears at least once.

prev = apply(X = otu_table(ps),
             MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prev_df = data.frame(Prevalence = prev,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

# filter out single hit taxa 
keepTaxa = rownames(prev_df)[(prev_df$Prevalence > 1)]
ps_noRare = prune_taxa(keepTaxa, ps)

saveRDS(ps_noRare, paste0(proc_data_dir, '/ps_norare.rds'))

## ps_noCorr: collapse ASVs that are only one base away ########################

# Find all ASV sequences ##
incl_asvs = as_tibble(as.data.frame(tax_table(ps)), rownames = 'ASV')
all_seq_df = as_tibble(taxa, rownames = 'seq')
all_seq_df['ASV'] = paste0('ASV', 1:dim(all_seq_df)[1])

# Compute all hamming distances between ASVs 
# simple hamming distance returns infinite sometimes because the sequences are not of equal length
distmat <- stringdistmatrix(all_seq_df$seq, all_seq_df$seq, method = 'hamming')
colnames(distmat) <- all_seq_df$ASV
rownames(distmat) <- all_seq_df$ASV
distmat[lower.tri(distmat, diag = T)] <- NA

# restructure the distance matrix to a long format
seq_dist_df <- as_tibble(distmat, rownames = 'ASV_left') %>%
  pivot_longer(-ASV_left, names_to = 'ASV_right', values_to = 'dist') %>%
  filter(!is.na(dist)) %>%
  left_join(all_seq_df[,c('Family', 'Genus','Species', 'ASV')], by = join_by('ASV_left' == 'ASV')) %>%
  left_join(all_seq_df[,c('Family', 'Genus','Species', 'ASV')], by = join_by('ASV_right' == 'ASV'), suffix = c('_left', '_right'))

# there are a number of sequences just 1 nucleotide away
pairs_dist1 <- seq_dist_df %>%
  filter(dist == 1) 

# create a df with the ASV abundance information of ps
ps_df <- psmelt(ps)

# List all pairs of ASVs with distance of 1 nucleotide + high abundance correlation ##
## this piece of code takes a while to run (> min)
ASV_allcor = data.frame()
for (ASV_left_choice in pairs_dist1 %>% pull(ASV_left) %>% unique()){
  ASV_l = ps_df %>%
    ungroup() %>%
    select(OTU, Abundance, rep, day, inoc, salt, temp) %>%
    filter(OTU == ASV_left_choice)
  
  ASV_r = ps_df %>%
    ungroup() %>%
    select(OTU, Abundance, rep, day, inoc, salt, temp) %>%
    filter(OTU %in% (pairs_dist1 %>% filter(ASV_left == ASV_left_choice) %>% pull(ASV_right))) %>%
    right_join(ASV_l, by = c('rep', 'day', 'inoc', 'salt', 'temp'), suffix = c('_right', '_left'))
  
  ASV_allcor = bind_rows(ASV_allcor, ASV_r)
}

strong_cor = ASV_allcor %>%
  group_by(OTU_left, OTU_right) %>%
  summarise(ASV_cor = cor(Abundance_right, Abundance_left, method = 'pearson'))  %>%
  arrange(-ASV_cor) %>%
  filter(ASV_cor >= 0.9) %>%
  left_join(ASV_allcor, by = c('OTU_left', 'OTU_right'))

#saveRDS(strong_cor, paste0(proc_data_dir, '/strong_cor.rds'))
strong_cor <- readRDS(paste0(proc_data_dir, '/strong_cor.rds'))

ggplot(strong_cor) +
  geom_point(aes(x = rel_abund_right, y = rel_abund_left, color = inoc )) +
  facet_wrap(vars(OTU_right), scales = 'free')

# make tibble with a list per asv of other asvs it is strongly correlated with
# compare this to the list per asv of blast hist to the isolates

foo = strong_cor %>%
  select(OTU_left, OTU_right, ASV_cor) %>%
  distinct()

ps_noCorr = ps
for (i in 1:dim(foo)[1]){
  otu_table(ps_noCorr)[, foo[[i, 'OTU_left']] ] <- otu_table(ps)[, foo[[i, 'OTU_right']] ] + otu_table(ps)[, foo[[i, 'OTU_left']] ]  
  otu_table(ps_noCorr)[, foo[[i, 'OTU_right']] ] <- 0
}

saveRDS(ps_noCorr, paste0(proc_data_dir, '/ps_nocorr.rds'))

## Create integrated data objects per ps (diversity, mcn) ############################################################
## Existing ps: ps, ps_rarefied, ps_noRare, ps_noCorr ####

ps_list = c(ps, ps_rarefied, ps_noRare, ps_noCorr)
#ps_list = c(ps_noRare)
ps_names = c('raw', 'rarefied', 'norare', 'nocorr')
#ps_names = c('norare')

## If starting from here, just load the ps dataframes into the ps_list
ps_list = list(ps)
for (i in 2:length(ps_names)){
  ps_name = ps_names[[i]]
  ps_list[[i]] <- readRDS(paste0(proc_data_dir, '/ps_', ps_name, '.rds'))
}



## Calculate diversity measures ###########################
#D0 = richness

diversity_df_list = list()
for (i in 1:length(ps_list)){
  ps_partition <- ps_list[[i]]
  ps_name <- ps_names[[i]]
  
  D0 <- estimate_richness(ps_partition,split=TRUE, measures=c("Observed", "Shannon"))

  diversity_df_list[[i]] = data.frame(sample = row.names(D0), D0) %>%
    separate(sample, into = c('plate', 'well'), sep = '_') %>%
    left_join(metadata, by = c('plate', 'well')) %>%
    pivot_longer(c('Observed', "Shannon"), names_to = 'diversity_metric', values_to = 'div_value')
  #print(diversity_df_list[[i]])
  
  write_csv(diversity_df_list[[i]], paste0(proc_data_dir, '/diversity_df_', ps_name, '.csv'))
}
names(diversity_df_list) <- ps_names

## Calculate rel abundance & mean copynumber ###############

rel_abund_df_list = list()
mcn_df_list = list()
df_names <- c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'GenusSpecies')

for (i in 1:length(ps_list)){
  ps_partition <- ps_list[[i]]
  ps_name <- ps_names[[i]]

  # turn ps object into a dataframe
  df_all <- psmelt(ps_partition) %>% 
    mutate(GenusSpecies = paste0(Genus, ' ', Species))
  
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
  print(ps_name)
  print('Check one, should read TRUE: ', dim(df_rrndb)[1] == dim(df_all)[1])
  print('Check two, should read TRUE: ', dim(df_rrndb %>% filter(is.na(copy_number_NCBI)))[1] == 0 )
  
  rel_abund_df_list[[i]] <- df_rrndb %>%
    mutate(weigh_abund = Abundance/copy_number_NCBI) %>%
    group_by(Sample) %>%
    mutate(rel_abund = weigh_abund / sum(weigh_abund, na.rm = T))
  write_csv(rel_abund_df_list[[i]], paste0(proc_data_dir, '/rel_abund_df_', ps_name, '.csv'))
  
  mcn_df_list[[i]] <- rel_abund_df_list[[i]] %>%
    group_by(rep, day, inoc, salt, temp) %>%
    summarise(mcn = weighted.mean(copy_number_NCBI, rel_abund, na.rm = T) )
  write_csv(mcn_df_list[[i]], paste0(proc_data_dir, '/mcn_df_', ps_name, '.csv'))
  
}

