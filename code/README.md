# Code used for preprocessing, analysis and plotting for the manuscript "Microbial communities demonstrate robustness in deteriorating environments due to predictable composition shifts"

Scripts are organized in order of appearance in the manuscript. 

Some scripts (clearly marked) will not run due to absence of the input data, but were included to be completely transparent on how more basic/low-level files were combined and pre-processed prior to the primary data analysis. These low-level files were not included because they are either numerous or too large for github, but are available from Jana S. Huisman upon request. 


Software requirements:
- R version 4.4.0 (or greater)
- R packages: tidyverse (v2.0.0), dada2 (v1.32.0), patchwork, readxl, segmented (v2.1-4), seqinr, ape, ggtree, deSolve, mgcv (v1.9-3), performance (v0.14.0), rnaturalearth (v1.0.1), sf (v1.0-21)
- Optional R packages: phyloseq, stringdist, stringr (v1.5.1), rBLAST
- Additional packages used in preprocessing: gcplyr (v1.11.0), xml2, openxlsx

Tested on Mac OS Ventura 13.6.6.

Suggestions for use: download the full github repository, or at least the folders containing data and code. Open the code file of interest in a development environment, e.g. Rstudio, and run the code line by line. Sections that are responsible for reading in data, processing it, or producing specific figures have been demarcated and commented. Chunks of code that are commented out typically take longer than 1-2 minutes to run, and a line of code to read in the produced file is included below that section. The figures folder includes all expected output figures of the code in this repository.


## 1. Community serial dilution experiment (Figs 1, S1-4, S6-7)
- plot_community_od_ph.R : analyses and plots the OD data of the serial passaging experiment, corresponding to Fig S4

- process_asvs.R : the main script to do 16S data preprocessing (cleaning ; rrn assignment). generates files ps_*.rds, diversity_df_*.rds, rel_abund_*.rds, mcn_df_*.rds . The "noRare" preprocessing is used in the rest of the analyses.
- plot_processed_asvs.R : plot diversity over time and at the end of the experiment
- plot_CCI_CGR_community_experiment.R : analyze and plot most panels of Fig 1


Additional files (included for completeness, will not run because individual datafiles are not included):
- preprocess_community_od_ph.R : combines OD and pH measurement files of the 7 serial dilution cycles into files all_od_data and all_ph_data
- asv_calling_and_taxonomy.R : runs dada2 on 16S sequences of community serial propagation experiment. Produces the files seqtab.nochim.rds , taxa.rds , metadata.csv


## 1.2 Isolate information and growth rates (Figs 1, S5, S8, S9, S15)
- process_isolate_data.R : assign taxonomy and rrn to isolates, assign isolates to asvs

- plot_isolate_growth.R : Plot growth related figures 1B, S8, S15 
- plot_phylogeny.R : plot phylogeny Fig S5, and the family level statistics shown in Fig 1C
- fit_growthrates.R : fit a 2 segment model to the measured salinity performance curves - produces fig S9

Additional files (included for completeness):
- read_od_xml.R : Transform xml files from the Tecan platereader into an interpretable csv file format.
- estimate_all_growthrates.R : the script used to estimate growthrates from raw OD data. Generates the files all_growthdata.csv and average_growthrates.csv . 


## 2. Modeling (Figs 2, S10-13)
- Fig2.R : script that generates and plots all simulation results
- model_pairwise.nb : mathematica file detailing the analytical derivations reported in the supplement

## 3. Pairwise competition (Figs 3, S14)
Main file:
- Fig3.R : analyses data and plots all panels of Fig 3 + S14

Additional files:
- preprocess_competition_data.R : included for completeness, will not run. Combines individual files with experimental count data into the file competition_plating_data.csv and the OD data into od_data.csv
- plot_OD_competition_data.R : plots the OD over time for the propagated pairwise competitions and corresponding monocultures. Data not shown in the manuscript.

## 4. Environmental data analysis (Fig 4)
Main analysis file:
- Fig4.R : computes GAM models for the 4 environmental datasets, plots all components of Fig. 4

Included for completeness. These scripts were used to preprocess the publicly available data, but will not run as provided because the raw data is too large to add to this github repository + we did not update file-paths. 
- preprocess_Kellogg_Beaufort.R
- analyze_Kellogg_Beaufort.R
- preprocess_Cram_Chesapeake.R
- analyze_Cram_Chesapeake.R
- analyze_Latz_Baltic.R