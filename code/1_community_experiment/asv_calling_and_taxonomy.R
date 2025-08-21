###########################################################
# Extracting ASVs from our 16S amplicon sequencing
# based on dada2 tutorial https://benjjneb.github.io/dada2/tutorial.html
# Author: JS Huisman
###########################################################

library(tidyverse)
library(ggplot2)
library(dada2)

###########################################################
## Read data ####
main_data_dir <- "../sequencing/dilution_experiment_Illumina/Jana_Martina_data"
proc_data_dir <- '../processed_data/daily_dilutions/sequencing'

# collect dirs with fastq files
data_dirs <- grep('RawData[1-9]/P', list.dirs(path = main_data_dir, recursive = T), value = T)

# Forward and reverse fastq filenames have format: WELL_1.fastq and WELL_2.fastq.gz
fnFs <- sort(list.files(data_dirs, pattern="[0-9]_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(data_dirs, pattern="[0-9]_2.fastq.gz", full.names = TRUE))
# Extract well names, assuming filenames have format: WELL_XXX.fastq
sample.names <- gsub("_[1-2].fastq.gz", "", basename(fnFs))


###########################################################
## Inspect read quality profiles ####

plotQualityProfile(fnFs[291:302], n = 10000, aggregate = T)

plotQualityProfile(fnRs[301:302], n = 10000)

## Surprisingly good quality? No real characteristic dropoff...

###########################################################
## Filter and trim ####
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(main_data_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(main_data_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,190),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose = T) 
head(out)

###########################################################
## Learn error rates ####

errF <- learnErrors(filtFs, multithread=TRUE)
#101763270 total bases in 484587 reads from 4 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#116144530 total bases in 611287 reads from 5 samples will be used for learning the error rates.

plotErrors(errF, nominalQ=TRUE)

###########################################################
## Sample inference! ####

## Tried adding pool=TRUE, but this takes too long
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Just checking the number of ASVs for the first sample
dadaFs[[1]]
dadaRs[[1]]

###########################################################
## Merge paired reads ####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###########################################################
## Construct sequence table ####

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 338 60034

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###########################################################
## Remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#338 16083
saveRDS(seqtab.nochim, file = paste0(proc_data_dir, "/seqtab.nochim.rds"))

sum(seqtab.nochim)/sum(seqtab)
#0.9306844

###########################################################
## Track pipeline ####

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## Might need to revisit primer removal??

#################################################################################################################################################################################
## Assign Taxonomy ####

# assign taxa
taxa <- assignTaxonomy(seqtab.nochim, paste0(main_data_dir, "/silva_nr99_v138.1_train_set.fa.gz"), multithread=TRUE)

# add exact species assignment 
taxa <- addSpecies(taxa, paste0(main_data_dir, "/silva_species_assignment_v138.1.fa.gz"), allowMultiple = T)
saveRDS(taxa, file = paste0(proc_data_dir, "/taxa.rds"))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

###########################################################
## Load metadata ####

library(readxl)

raw_meta_P14 = read_xlsx(paste0(main_data_dir, '/../JSH_MDB_sequencing_map_modified.xlsx'), range = c('A3:M11')) %>%
  rename(row_letter = '...1') %>%
  pivot_longer(cols = !row_letter, names_to = 'col_number', values_to = 'sample_description')

raw_meta_P15 = read_xlsx(paste0(main_data_dir, '/../JSH_MDB_sequencing_map_modified.xlsx'), range = c('A14:M22')) %>%
  rename(row_letter = '...1') %>%
  pivot_longer(cols = !row_letter, names_to = 'col_number', values_to = 'sample_description')

raw_meta_P16 = read_xlsx(paste0(main_data_dir, '/../JSH_MDB_sequencing_map_modified.xlsx'), range = c('A25:M33')) %>%
  rename(row_letter = '...1') %>%
  pivot_longer(cols = !row_letter, names_to = 'col_number', values_to = 'sample_description')

raw_meta_P17 = read_xlsx(paste0(main_data_dir, '/../JSH_MDB_sequencing_map_modified.xlsx'), range = c('A36:M44')) %>%
  rename(row_letter = '...1') %>%
  pivot_longer(cols = !row_letter, names_to = 'col_number', values_to = 'sample_description')

metadata <- bind_rows(P14 = raw_meta_P14, P15 = raw_meta_P15, P16 = raw_meta_P16, P17 = raw_meta_P17, .id = 'plate') %>%
  mutate(well = paste0(row_letter, col_number)) %>%
  separate_wider_delim(sample_description, delim = " ", names = c("rep", "day", "inoc", "salt", "temp"), too_few = "align_start") %>%
  filter(rep != 'empty')

write_csv(metadata, file = paste0(proc_data_dir, '/metadata.csv'))







