###########################################################
## Extracting and preprocessing the data from Beaufort sea
## Author: JS Huisman
###########################################################

library(tidyverse)
library(readxl)
library(dada2)

main_data_dir <- "../sequencing/dilution_experiment_Illumina/Jana_Martina_data"

data_dir = '../sequencing/environmental_data/Kellogg_Beaufort'

###########################################################
acc_list = read_csv(paste0(data_dir, '/SRR_Acc_List.txt'), col_names = 'acc_no')
run_meta_df = read_csv(paste0(data_dir, '/SraRunTable.txt'))

metadata <- run_meta_df %>%
  select(sample = Run, sample_name = `Sample Name`, date = Collection_Date, 
         depth = Depth, lat_lon, geo_loc_name, site_name, env_medium, Organism) 

write_csv(metadata, file = paste0(data_dir, '/processed_data/metadata.csv'))

######################################################################################################################

## Forward and reverse fastq filenames ####
fnFs <- sort(list.files(paste0(data_dir, '/raw_data'), pattern="_1.fastq.gz", full.names = T))
fnRs <- sort(list.files(paste0(data_dir, '/raw_data'), pattern="_2.fastq.gz", full.names = T))
## Extract sample names ####
sample.names <- sort(acc_list$acc_no)

###########################################################
## Inspect read quality profiles ####

plotQualityProfile(fnFs[1:5], n = 5e+05)
plotQualityProfile(fnRs[1:10], n = 5e+05)

###########################################################
## Filter and trim ####
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(data_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(data_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft = 20,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose = T)
head(out)

###########################################################
## Learn error rates ####

errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errR, nominalQ=TRUE)

###########################################################
## Sample inference! ####

## Tried adding pool=TRUE, but this takes too long
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Just checking the number of ASVs for the first sample
dadaFs[[1]]
dadaRs[[1]]
# 299 sequence variants were inferred from 18549 input unique sequences.

###########################################################
## Merge paired reads ####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###########################################################
## Construct sequence table ####

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 66 245068

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


###########################################################
## Remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#66 30146
saveRDS(seqtab.nochim, file = paste0(data_dir, "/processed_data/seqtab.nochim.rds"))

sum(seqtab.nochim)/sum(seqtab)
#0.4422776

###########################################################
## Track pipeline ####

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#################################################################################################################################################################################
## Assign Taxonomy ####

# assign taxa
taxa <- assignTaxonomy(seqtab.nochim, paste0(main_data_dir, "/silva_nr99_v138.1_train_set.fa.gz"), multithread=TRUE)

# add exact species assignment 
taxa <- addSpecies(taxa, paste0(main_data_dir, "/silva_species_assignment_v138.1.fa.gz"), allowMultiple = T)
saveRDS(taxa, file = paste0(data_dir, "/processed_data/taxa.rds"))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print, n = 15)
