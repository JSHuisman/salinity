################################################################################
## Plot phylogenetic tree with some additional isolate information
## Author: JS Huisman
################################################################################

library(tidyverse)
# library(readxl)
library(seqinr)
# library(stringr)
# library(dada2) #bioconductor
# library(phyloseq)
# library(rBLAST) #bioconductor
library(ape)
#BiocManager::install("msa")
#library(msa)
library(ggtree)
library(patchwork)

#raw_data_dir = '../raw_data/isolates/'
theme_set(theme_minimal() + 
            theme(text = element_text(size = 20)))

##### Growth data #####################################################

proc_data_dir = '../data/1.2_isolates/'
fig_dir = '../figures/1.2_isolates/'

average_growth_df <- read.csv('../data/1.2_isolates/average_growthrates.csv')

isolates_growth <- average_growth_df %>%
  filter(n_growth_meas >= 3) %>%
  pull(isolate_id) %>% unique()



##### MAFFT / RAXML tree #####################################################
# Load data
fasta_16S = '../data/1.2_isolates/all_isolates.fasta'
all_16Sseqs <- read.fasta(fasta_16S)
best_tree = ape::read.tree("../data/1.2_isolates/RAxML_bestTree.all_isolates_raxml")

Tree = best_tree

ggtree(best_tree, ladderize = F)  + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) +
  geom_treescale() 

######
# growth data for tree
growth_Tree <- ape::drop.tip(Tree, setdiff(names(all_16Sseqs), isolates_growth))
growth_Tree = root(growth_Tree, node=98, resolve.root=T)


growth_forplot <- average_growth_df %>%
  filter(n_growth_meas >= 3) %>%
  pivot_wider(id_cols = c(isolate_id, class, order, family, genus, species), names_from = well_salt, names_prefix = 'salt', values_from = growthrate) %>%
  arrange(class, order, family, genus, species)

growth_df_forplot <- data.frame(growth_forplot[, 7:18])
rownames(growth_df_forplot) <- growth_forplot$isolate_id
colnames(growth_df_forplot) <- gsub('salt', '', colnames(growth_df_forplot))

# number of isolates per family
growth_forplot %>%
  group_by(family) %>%
  summarize(count = n())

#######
## display phylogenetic tree
options(ignore.negative.edge=TRUE)

# Identify internal nodes that are the MRCA of large families
pseudo_ids = growth_forplot %>% filter(family == 'Pseudoalteromonadaceae') %>% pull(isolate_id) %>% unique()
pseudo_mrca = MRCA(growth_Tree, pseudo_ids)

shewa_ids = growth_forplot %>% filter(family == 'Shewanellaceae') %>% pull(isolate_id) %>% unique()
shewa_mrca = MRCA(growth_Tree, shewa_ids)

flavo_ids = growth_forplot %>% filter(family == 'Flavobacteriaceae') %>% pull(isolate_id) %>% unique()
flavo_mrca = MRCA(growth_Tree, flavo_ids)

morax_ids = growth_forplot %>% filter(family == 'Moraxellaceae') %>% pull(isolate_id) %>% unique()
morax_mrca = MRCA(growth_Tree, morax_ids)

pseudom_ids = growth_forplot %>% filter(family == 'Pseudomonadaceae') %>% pull(isolate_id) %>% unique()
pseudom_mrca = MRCA(growth_Tree, pseudom_ids)

ocean_ids = growth_forplot %>% filter(family == 'Oceanospirillaceae') %>% pull(isolate_id) %>% unique()
ocean_mrca = MRCA(growth_Tree, ocean_ids)

rhodo_ids = growth_forplot %>% filter(family == 'Rhodobacteraceae') %>% pull(isolate_id) %>% unique()
rhodo_mrca = MRCA(growth_Tree, rhodo_ids)

baci_ids = growth_forplot %>% filter(family == 'Bacillaceae') %>% pull(isolate_id) %>% unique()
baci_mrca = MRCA(growth_Tree, baci_ids)

coll_ids = growth_forplot %>% filter(family == 'Colwelliaceae') %>% pull(isolate_id) %>% unique()
coll_mrca = MRCA(growth_Tree, coll_ids)


color_list_families = c("Pseudoalteromonadaceae" = '#6E969B', 
                        "Moraxellaceae"  = '#C0E49B', "Pseudomonadaceae"  = '#82E2F3', 
                        "Flavobacteriaceae" = '#F5BD23',"Shewanellaceae"= '#741B7B',
                        "Oceanospirillaceae" = '#004D07',
                        "Rhodobacteraceae"  = '#A72626',
                        "Bacillaceae" = '#B77406',  "Colwelliaceae" = '#356FEA')


p <- ggtree(growth_Tree, ladderize = F)  + 
  #geom_tiplab(size=2, align=TRUE, linesize=.5) +
  geom_treescale() + #geom_text(aes(label=node), hjust=-1) +
  geom_cladelabel(node=pseudo_mrca, label="Pseudo-\nalteromonadaceae", fill="#6E969B", align = T) +
  geom_hilight(node=pseudo_mrca, fill="#6E969B") +
  geom_cladelabel(node=flavo_mrca, label="Flavobacteriaceae", color="#F5BD23", align = T) +
  geom_hilight(node=flavo_mrca, fill="#F5BD23") +
  geom_cladelabel(node=shewa_mrca, label="Shewanellaceae", color='#741B7B', align = T) +
  geom_hilight(node=shewa_mrca, fill='#741B7B') +
  geom_cladelabel(node=pseudom_mrca, label="Pseudomonadaceae", color='#82E2F3', align = T) +
  geom_hilight(node=pseudom_mrca, fill='#82E2F3') +
  geom_cladelabel(node=morax_mrca, label="Moraxellaceae", color='#C0E49B', align = T) +
  geom_hilight(node=morax_mrca, fill='#C0E49B') +
  geom_cladelabel(node=ocean_mrca, label="Oceanospirillaceae", color='#004D07', align = T) +
  geom_hilight(node=ocean_mrca, fill='#004D07')+
  geom_cladelabel(node=rhodo_mrca, label="Rhodobacteraceae", color='#A72626', align = T) +
  geom_hilight(node=rhodo_mrca, fill='#A72626')+
  geom_cladelabel(node=baci_mrca, label="Bacillaceae", color='#B77406', align = T) +
  geom_hilight(node=baci_mrca, fill='#B77406')+
  geom_cladelabel(node=coll_mrca, label="Colwelliaceae", color='#356FEA', align = T) +
  geom_hilight(node=coll_mrca, fill='#356FEA')
  
p

gheatmap(p, growth_df_forplot, offset = 0.5, width=1, 
         colnames=TRUE, colnames_offset_y = -0.4) +
  theme(legend.position = 'right') +
  labs(fill = 'Growth \nrate (1/h)') +
  scale_fill_viridis_c() 


ggsave('../figures/1.2_isolates/FigS5_RaxML_phylogeny_growth.pdf', height = 6, width = 10)

###########
## Plot sopt and rmax for families
max_growth_salinity_df <- average_growth_df %>%
  filter(n_growth_meas >= 3,
         !is.na(family)) %>%
  group_by(class, order, family, genus, species, isolate_id) %>%
  summarize(max_gr = max(growthrate, na.rm = T),
            sopt = well_salt[growthrate == max_gr],
            rdiff_3045 = growthrate[well_salt == 45] - growthrate[well_salt == 30]) %>%
  group_by(family) %>%
  summarize(n_isolate = n(),
            mean_max_gr = mean(max_gr, na.rm = T),
            sd_max_gr = sd(max_gr, na.rm = T),
            mean_sopt = mean(sopt, na.rm = T),
            sd_sopt = sd(sopt, na.rm = T),
            mean_rdiff3045 = mean(rdiff_3045, na.rm = T),
            sd_rdiff3045 = sd(rdiff_3045, na.rm = T)) %>%
  filter(n_isolate >= 3) %>%
  mutate(ordered_family_id = factor(family, levels = c("Pseudomonadaceae" , "Moraxellaceae", "Oceanospirillaceae" ,
                                                       "Shewanellaceae","Pseudoalteromonadaceae" ,"Colwelliaceae", 
                                                       "Flavobacteriaceae", "Bacillaceae" ,  "Rhodobacteraceae"  )))




max_growth_plot <- ggplot(max_growth_salinity_df) +
  geom_point(aes(y = ordered_family_id, x = mean_max_gr, color = ordered_family_id), size = 6, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = ordered_family_id, xmin = pmax(mean_max_gr - sd_max_gr, 0), xmax = mean_max_gr + sd_max_gr, color = ordered_family_id), 
                linewidth = 2, width = 0, position = position_dodge(width = 1)) +
  labs(y = '', x = latex2exp::TeX('$r_{max}$ (1/h)'), color = 'Family') +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = color_list_families) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'none')

max_growth_plot
ggsave('../figures/1.2_isolates/Fig1C_family_growthrate.pdf', width = 8, height = 5)


sopt_plot <- ggplot(max_growth_salinity_df) +
  geom_point(aes(y = ordered_family_id, x = mean_sopt, color = ordered_family_id), size = 6, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = ordered_family_id, xmin = pmax(mean_sopt - sd_sopt, 0), xmax = mean_sopt + sd_sopt, color = ordered_family_id), 
                linewidth = 2, width = 0, position = position_dodge(width = 1)) +
  geom_vline(xintercept = 35, linetype = 'dashed') +
  scale_y_discrete(limits=rev) +
  labs(y = '', x = latex2exp::TeX('$s_{opt}$ (g/L)'), color = 'Family') +
  scale_color_manual(values = color_list_families) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'none')

sopt_plot
ggsave('../figures/1.2_isolates/Fig1C_family_sopt.pdf', width = 8, height = 5)


rdiff3045_plot <- ggplot(max_growth_salinity_df) +
  geom_point(aes(y = ordered_family_id, x = mean_rdiff3045, color = ordered_family_id), size = 6, position = position_dodge(width = 1)) +
  geom_errorbar(aes(y = ordered_family_id, xmin = mean_rdiff3045 + sd_rdiff3045, xmax = mean_rdiff3045 - sd_rdiff3045, color = ordered_family_id), 
                linewidth = 2, width = 0, position = position_dodge(width = 1)) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks = c(0, -0.15, -0.3)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(y = '', x = latex2exp::TeX('$Delta r_{45 - 30}$ (1/h)'), color = 'Family') +
  scale_color_manual(values = color_list_families) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'none')

rdiff3045_plot
ggsave('../figures/1.2_isolates/Fig1C_family_rdiff3045.pdf', width = 8, height = 5)


(sopt_plot + theme(axis.text.y = element_text(size = 15))) +
  (max_growth_plot + theme(axis.text.y = element_blank())) +
  (rdiff3045_plot + theme(axis.text.y = element_blank())) &
  theme(axis.text = element_text(size = 15))
ggsave('../figures/1.2_isolates/Fig1C_family_sopt_maxgrowth.pdf', width = 9, height = 4)

