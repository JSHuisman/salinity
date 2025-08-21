###########################################################
## Extracting and preprocessing the data from the Chesapeake bay
## Author: JS Huisman
###########################################################

library(tidyverse)
library(readxl)
library(dada2)

main_data_dir <- "../sequencing/dilution_experiment_Illumina/Jana_Martina_data"
data_dir = '../sequencing/environmental_data/Cram_Chesapeake'

###########################################################
acc_list = read_csv(paste0(data_dir, '/SRR_Acc_List.txt'), col_names = 'acc_no')
run_meta_df = read_csv(paste0(data_dir, '/SraRunTable.txt'))
# meta_df = read_xlsx(paste0(data_dir, '/metadata.xlsx')) %>%
#   mutate(name = tolower(gsub('_', '', name)))

download_list = list.files(paste0(data_dir, '/Paired-end_data'), full.names = F)

#write_csv(combi_metadata, file = paste0(data_dir, '/processed_data/metadata.csv'))

######################################################################################################################

## collect dirs with fastq files ####
data_dirs <- list.dirs(path = paste0(data_dir, '/Paired-end_data'))
## Forward and reverse fastq filenames ####
fnFs <- sort(list.files(data_dirs, pattern="forward.fastqsanger.gz", full.names = TRUE))
fnRs <- sort(list.files(data_dirs, pattern="reverse.fastqsanger.gz", full.names = TRUE))
## Extract sample names ####
sample.names <- sort(grep('S', list.dirs(path = paste0(data_dir, '/Paired-end_data'), full.names = F), value = T))

###########################################################
## Inspect read quality profiles ####

plotQualityProfile(fnFs, n = 10000)

#Reverse basically doesn't exist!
plotQualityProfile(fnRs, n = 10000)

###########################################################
## Filter and trim ####
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(data_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(data_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Did you remove primers from your raw reads? This needs to be done before running dada2, 
#or can be performed with the trimLeft parameter of the filterAndTrim function.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220, 220),
                     trimLeft = 30,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose = T) 
head(out)

###########################################################
## Learn error rates ####

errF <- learnErrors(filtFs, multithread=TRUE)
#125227740 total bases in 569217 reads from 6 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)
#125227740 total bases in 569217 reads from 6 samples will be used for learning the error rates.
plotErrors(errR, nominalQ=TRUE)

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
# 110 140730

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


###########################################################
## Remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 93311 bimeras out of 140730 input sequences.
dim(seqtab.nochim)
#110 23513
saveRDS(seqtab.nochim, file = paste0(data_dir, "/processed_data/seqtab.nochim.rds"))

sum(seqtab.nochim)/sum(seqtab)
#0.8393171

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
head(taxa.print)
