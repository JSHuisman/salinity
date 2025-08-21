################################################################################
## Process isolate data 
## Author: JS Huisman
################################################################################

library(tidyverse)
library(readxl)
library(seqinr)
library(stringr)
library(dada2) #bioconductor
library(phyloseq)
library(rBLAST) #bioconductor

################################################################################
## Assign Taxonomy ####
# proc_data_dir = paste0('../data/1.2_isolates')
# 
# all_seq_data <- read.fasta(paste0(proc_data_dir, '/all_isolates.fasta'), as.string = T, forceDNAtolower = F)
# all_seqs <- sapply(unlist(all_seq_data), function(x){toupper(x)})
# all_seq_names <- names(all_seq_data)
# 
# id_blocks = list(1:20, 21:41, 42:63, 64:85, 86:107, 108:129, 130:153)
# 
# # assign taxa 
# all_taxa <- data.frame()
# for(i in 1:7){
#   
# taxa <- assignTaxonomy(all_seqs[id_blocks[[i]]], "../sequencing/dilution_experiment_Illumina/Jana_Martina_data/silva_nr99_v138.1_train_set.fa.gz", 
#                        multithread=TRUE)
# 
# pick_longest_nonN = function(x){
#   pieces = str_split_1(x, 'N')
#   lengths = sapply(pieces, str_length)
#   return(pieces[which.max(lengths)])
# }
# 
# rownames(taxa) <- sapply(rownames(taxa), pick_longest_nonN)
# 
# # add exact species assignment 
# taxa <- addSpecies(taxa, "../sequencing/dilution_experiment_Illumina/Jana_Martina_data/silva_species_assignment_v138.1.fa.gz", 
#                    allowMultiple = T)
# 
# all_taxa <- bind_rows(all_taxa, as_tibble(as.matrix(taxa)))
# print(i)
# }
# all_taxa['id'] = all_seq_names
# write_csv(all_taxa, file = paste0(proc_data_dir, "/taxa_isolates.csv"))

################################################################################
### Assign copynumber to the isolates - up to Family level ######

taxa_isolates <- read_csv(file = paste0(proc_data_dir, "/taxa_isolates.csv"))

rrndb_NCBI <- read_delim(paste0(proc_data_dir, "/rrnDB-5.8_pantaxa_stats_NCBI.tsv"), delim = '\t') %>%
  select (rank, name, copy_number_NCBI = mean)
rrndb_names <- c("superkingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species')

df_names <- c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'GenusSpecies')

df_rrndb <- taxa_isolates %>%
  rowwise() %>%
  mutate(GenusSpecies = list(paste0(Genus, ' ', unlist(str_split(Species, '/'))) )) %>%
  unnest_longer(GenusSpecies) %>%
  dplyr::rename(name = df_names[7] ) %>%
  left_join(rrndb_NCBI %>% filter(rank == rrndb_names[7]), by = 'name') %>%
  dplyr::rename(!!df_names[7] := name) # this just renames the column back to the taxon level

# Assign copynumber up to Family level
for(level in 6:5){
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

isolate_rrndb <- df_rrndb %>% 
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, id) %>%
  summarize(copy_number_NCBI = mean(copy_number_NCBI),
            #rank = rrndb_names[max(which(rank == rrndb_names))]
            rank = list(rank),
            .groups = 'drop') %>%
  left_join(taxa_isolates, by = c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', "id"))


write_csv(isolate_rrndb, file = paste0(proc_data_dir, "/isolate_rrndb.csv"))
################################################################################
#### Assign isolates to ASVs ####

proc_data_dir = '../data/1_community_experiment/'

## To run BLAST from R:
## use your conda environment in Rstudio by navigating to your project directory in the terminal. 
## there, type "open project.Rproj" and Rstudio will launch using your environment. 

## Note: no spaces are allowed in the file path of your blast database

#conda activate default

####### Make an isolate database ################################
#makeblastdb -in all_isolates.fasta -parse_seqids -dbtype nucl -out isolatedb

# SAVE in location blast_db_loc

####### Read in ASVs ################################
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

######### Query all ASVs against isolate db  #############

# load isolate database
bl <- blast(db="blast_db_loc")

# load asv seqs
seq <- refseq(ps)

# Blast ASVs against isolate database
cl <- predict(bl, refseq(ps), BLAST_args = "-perc_identity 97")

clean_blast_results = cl %>%
  rename(sseqid = 'isolate_id',
         qseqid = 'QueryID',
         pident = 'Perc.Ident',
         length = 'Alignment.Length',
         bitscore = 'Bits')

# Some ASVs match multiple isolates equally well
length(unique(clean_blast_results$QueryID))

write_csv(clean_blast_results, paste0(proc_data_dir, "/isolate_blast.csv"))

#### Reformat BLAST results for future use  ################################
#clean_blast_results <- read_csv(paste0(proc_data_dir, "/isolate_blast.csv"))

# for each ASV determine the best possible matching percentage/bitscore
# store all isolates that match the ASV with that score in the matching_isolates column

asv_isolate_map <- clean_blast_results %>%
  group_by(QueryID, Perc.Ident, Bits) %>%
  summarize(matching_isolates = paste(isolate_id, collapse = "/") ) %>%
  group_by(QueryID) %>%
  arrange(-Bits) %>%
  filter(row_number()==1) %>%
  mutate(ASV_number = as.numeric(gsub('ASV', '', QueryID))) %>%
  arrange(ASV_number)

write_csv(asv_isolate_map, paste0(proc_data_dir, "/asv_isolate_map.csv"))


################################################################################
# Load growth data
all_results_df <- read_csv('../data/1.2_isolates/all_growthdata.csv')

S1_table <- isolate_info[, c("isolate_id", "plate_date", "inoc", "plate_temperature","plate_salinity","morphology" , "color" )] %>%
  left_join(isolate_rrndb, by = join_by(isolate_id == id)) %>%
  select(-rank) %>%
  left_join(n_growth_meas, by = 'isolate_id') %>%
  filter(n >= 3) %>%
  rename(n = 'n_growth_meas')

write_csv(S1_table, paste0(proc_data_dir, '../S1_table.csv'))

################################################################################
#### Stats on the number of isolates #####
isolate_info <- read_csv('../data/1.2_isolates/paper_isolate_ids.csv')
all_taxa <- read_csv(file ="../data/1.2_isolates/taxa_isolates.csv")

integrated_isolate_info <- isolate_info %>%
  select(-class, -order, -family, -genus, -species, -morphology, -color, -plate_storage) %>%
  left_join(all_taxa %>% dplyr::rename("isolate_id" = id), by = 'isolate_id') %>%
  filter(plate_date != 'Cordero',
         sanger == 'x')

# Number of sanger sequenced isolates
length(unique(integrated_isolate_info$isolate_id))
# Number of unique Families
length(unique(integrated_isolate_info$Family))
# Number of unique Genus-Species assignments
length(unique(integrated_isolate_info$Genus))
# Number of unique Genus-Species assignments
length(unique(integrated_isolate_info %>% mutate(GS = paste0(Genus, ' ', Species)) %>% pull(GS)))

# Number of isolates with at least 3 growth measurements
all_results_df <- read_csv('../data/1.2_isolates/all_growthdata.csv')
length(unique(average_growth_df %>% filter(n_growth_meas >= 3) %>% pull(isolate_id)))