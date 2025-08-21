######################################################################################################################
## Simulating community dynamics as a function of salinity
## This code generates panels for Figures 2, S10-13 of the paper
## Author: JS Huisman
######################################################################################################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(deSolve)

fig_dir = '../figures/2_model/'

theme_set(theme_minimal() + theme(text = element_text(size = 20)))


###### Define Model Functions ################################################################

# pops: named list of numerics; initial population sizes pop.Xi
# parms: named list of numerics; (r_vec, alpha_m, death)
model_many_species <- function(time, pops, parms){
  with(as.list(parms), {
    pops[pops < 0] <- 0
    
    deriv.pops <- r_vec*pops*(1 - alpha_m%*%pops) - death*pops
    list(deriv.pops);
  })
}

integrate_model <- function(pops, t.vec, parms, sum_stat = F){
  
  model_out <- as.data.frame(deSolve::lsoda(y = pops, times = t.vec,
                                            func = model_many_species, parms = parms))
  
  # only return "summary statistic", namely the last value
  if(sum_stat){
    model_out = model_out[length(t.vec), 2:(length(pops)+1)]
  }
  return(model_out)
}


run_scenario <- function(salt_val, r_vec, r_salt, alpha_m, death, t.vec = 1:5000, n_species = 50, existence_threshold = 1e-6){
  pops = rep(1/n_species, n_species)
  r_rank_orig = order(r_vec, decreasing = T)
  r_rank_new = order(r_salt, decreasing = T)
  
  parms <- list(r_vec = r_salt, alpha_m = alpha_m, death = death)
  out = integrate_model(pops, t.vec, parms, sum_stat = T)
  
  abundance_vec = unlist(out, use.names = FALSE)
  abundance_vec[abundance_vec < existence_threshold] = 0
  
  result_df <- data.frame(salt = salt_val, r_id = 1:n_species, orig_r = r_vec, orig_rank = factor(r_rank_orig),
                          salt_r = r_salt, salt_rank = factor(r_rank_new),
                          abundance = abundance_vec,
                          mean_salt_r = max(weighted.mean(r_salt, abundance_vec), 0, na.rm = T),
                          mean_r = max(weighted.mean(r_vec, abundance_vec), 0, na.rm =T),
                          survival_frac = sum(abundance_vec > existence_threshold)/n_species,
                          tot_biomass = sum(abundance_vec) )
  
  return(result_df)
}


analytic_2species_ratio <- function(salt_val, r_salt, alpha_m, death, slopes, r_orig, type = 'parallel', smax = 40){

  if(r_orig[1] == 0 | r_orig[2] == 0 | r_orig[2] >= r_orig[1]){
    n1n2_ratio = NA
  } else if (type == 'parallel'){
    n1n2_ratio = ( (alpha_m[1,2]*alpha_m[2,1] - alpha_m[1,1]*alpha_m[2,2])*death*slopes*
      (r_orig[1]-r_orig[2])*(death - r_salt[1] - r_salt[2]) ) / (alpha_m[2,1]*r_salt[2]*
        (death - r_salt[1]) - alpha_m[1,1]*r_salt[1]*(death - r_salt[2]) )^2
  } else if (type == 'convergent'){
    n1n2_ratio = ( -(alpha_m[1,2]*alpha_m[2,1] - alpha_m[1,1]*alpha_m[2,2])*death*r_orig[1]*r_orig[2]*
                     (r_orig[1]-r_orig[2]) ) / (smax*(alpha_m[2,1]*r_orig[2]*
                      (death - r_salt[1]) - alpha_m[1,1]*r_orig[1]*(death - r_salt[2]) )^2)
  } else if (type == 'divergent'){
    n1n2_ratio = ( (alpha_m[1,2]*alpha_m[2,1] - alpha_m[1,1]*alpha_m[2,2])*death*(slopes[1] - slopes[2])*
                     (r_orig[1]^2 - r_orig[1]*death - slopes[1]*slopes[2]*salt_val^2)) / (alpha_m[2,1]*r_salt[2]*
                   (death - r_salt[1]) - alpha_m[1,1]*r_salt[1]*(death - r_salt[2]) )^2
    }
  return(n1n2_ratio)
}


###### 2-species example ##################################################################

r_vec <- c(0.7, 0.4)
death = 0.2
alpha_m <- matrix(c(1, -0.8, 0.8, 1), nrow = 2)

#### Scenario 1: #####
# growth rates decrease linearly at same rate from 
# a sampled rmax, until they hit 0
slope = -0.01
smax = 40

n2_result_df <- data.frame()
analytic_n1n2 <- vector(mode = 'numeric', length = smax+1)

for (salt_val in 0:smax){
  r_salt = pmax(r_vec + slope*salt_val, 0)
  result_df <- run_scenario(salt_val, r_vec, r_salt, alpha_m, death, t.vec = 1:5000, n_species = 2, existence_threshold = 1e-6)

  n2_result_df = bind_rows(n2_result_df, result_df)
  
  analytic_n1n2[salt_val+1] <- analytic_2species_ratio(salt_val, r_salt, alpha_m, death, -slope, 
                                                       r_orig = r_vec, type = 'parallel', smax)
}

crash2_salt = n2_result_df %>% filter(orig_rank == 2, abundance == 0) %>% pull(salt) %>% min()
analytic_n1n2 <- data.frame(ratio = analytic_n1n2, salt = 0:smax)

# Growthrate
commgrowth_plot <- ggplot(n2_result_df) +
  geom_vline(xintercept = crash2_salt, linetype = 'dashed') +
  geom_hline(yintercept = death, linetype = 'dashed') +
  geom_line(aes(x = salt, y = salt_r, group = orig_rank, color = orig_rank), 
            linewidth = 2, show.legend = F) +
  geom_line(aes(x = salt, y = mean_salt_r),
            linewidth = 3, show.legend = F) +
  scale_color_manual(values = c('#A72626', '#F5BD23')) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 0.8)) +
  labs(x = 'Salinity (g/L)', y = 'Growth rate (1/h)') +
  theme_minimal() + 
  theme(text = element_text(size = 20))

commgrowth_plot

ggsave(paste0(fig_dir, 'Fig2B.pdf'), width = 6, height = 4)

# Relative Abundance
relabund_plot <- ggplot(n2_result_df) +
  geom_vline(xintercept = crash2_salt, linetype = 'dashed') +
  geom_line(aes(x = salt, y = abundance / tot_biomass, colour = orig_rank, 
                group = orig_rank), linewidth = 2, show.legend = F) +
  scale_color_manual(values = c('#A72626', '#F5BD23')) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 1)) +
  labs(x = 'Salinity (g/L)', y = 'Relative abundance') +
  theme_minimal() +
  theme(text = element_text(size = 20),
  )

relabund_plot

ggsave(paste0(fig_dir, 'Fig2C.pdf'), width = 6, height = 4)


#### Use analytical solution to scan scenario 1 across different parameters #####

alpha_m <- matrix(c(1, -0.8, 0.8, 1), nrow = 2)

r2 = 0.6
analytic_n1n2 <- crossing(r1 = seq(r2, 1.2, 0.1), death_val = seq(0.1, 0.5, 0.1), slope_val = seq(0.005, 0.04, 0.001)) %>%
  mutate(ratio_0 = NA, ratio_10 = NA, ratio_20 = NA, ratio_30 = NA, ratio_40 = NA)


for (row_id in 1:dim(analytic_n1n2)[1]){

    r_vec <- c(analytic_n1n2[[row_id, 'r1']], r2)

    for (salt_val in c(0, 10, 20, 30, 40)){
      analytic_n1n2[row_id, paste0('ratio_', salt_val)] <- analytic_2species_ratio(salt_val, r_salt = pmax(r_vec + analytic_n1n2[[row_id, 'slope_val']]*salt_val, 0), 
                                                                   alpha_m, analytic_n1n2[[row_id, 'death_val']], analytic_n1n2[[row_id, 'slope_val']], 
                                                                   r_orig = r_vec, type = 'parallel')  
    }
    
}


long_df <- analytic_n1n2 %>% 
  pivot_longer(cols = starts_with('ratio_'), names_prefix = 'ratio_', names_to = 'salt', values_to = 'deriv')

# Plot only at salinity = 0, color slopes
ggplot( long_df %>% filter(salt ==0, death_val != 0.1)) +
  geom_line(aes(x = r1, color = slope_val, y = deriv, group = interaction(slope_val, salt)) , alpha = 0.7) +
  facet_wrap(vars(death_val)) +
  labs(x = 'Growth rate Species 1 (1/h)', color = 'Slope', y = 'd/ds N1/N2') +
  scale_color_viridis_c() +
  theme(legend.position = 'bottom',
        panel.spacing.x = unit(2, 'lines'),
        legend.key.width = unit(2, 'cm'))

ggsave(paste0(fig_dir, 'FigS11.pdf'), width = 10, height = 6)


####### 50-species examples: 3 scenarios #####################################################################

n_species = 50
r_vec <- rnorm(n_species, mean = 0.8, sd = 0.4)
death = 0.3

alpha_m <- matrix(runif(n_species^2, min = 0, max = 0.5), ncol=n_species)
diag(alpha_m) <- 1

#### Scenario 1: #####
# growth rates decrease linearly at same rate from 
# a sampled rmax, until they hit 0
slope = -0.02
smax = 40

scen1_result_df <- data.frame()
for (salt_val in 0:smax){
  r_salt = pmax(r_vec + slope*salt_val, 0)
  result_df <- run_scenario(salt_val, r_vec, r_salt, alpha_m, death, t.vec = 1:5000, n_species = 50, existence_threshold = 1e-6)
  
  scen1_result_df = bind_rows(scen1_result_df, result_df)
}


#### Scenario 2: #####
# growth rates decrease linearly from a sampled rmax down to smax
# rmax/smax = slope
smax = 40
slopes = - (r_vec-death) /(smax-5)

scen2_result_df <- data.frame()
for (salt_val in 0:smax){
  r_salt = r_vec + slopes*salt_val

  result_df <- run_scenario(salt_val, r_vec, r_salt, alpha_m, death, t.vec = 1:5000, 
                            n_species = 50, existence_threshold = 1e-6)

  scen2_result_df = bind_rows(scen2_result_df, result_df)
}

#### Scenario 3: #####
# growth rates decrease linearly from the same rmax down to a sampled smax
# rmax/smax = slope

smax = runif(n_species, min = 20, max = 60)
slopes = - max(r_vec) /smax

scen3_result_df <- data.frame()
for (salt_val in 0:40){
  r_salt = pmax(max(r_vec) + slopes*salt_val, 0)

  result_df <- run_scenario(salt_val, r_vec, r_salt, alpha_m, death, 
                            t.vec = 1:5000, n_species = 50, existence_threshold = 1e-6)
  
  scen3_result_df = bind_rows(scen3_result_df, result_df)
}

##### Plot the results for the 50-species #######################

results_list <- list(scen1_result_df, scen2_result_df, scen3_result_df)

for (scen_id in 1:3){
  # Abundance / Biomass
  absabund_plot <- ggplot(results_list[[scen_id]]) +
    geom_line(aes(x = salt, y = abundance, group = orig_rank),  color = 'grey', linewidth = 0.5, show.legend = F) +
    geom_line(aes(x = salt, y = tot_biomass), color = 'black', linewidth = 3) +
    labs(x = 'Salinity (g/L)', y = 'Biomass') +
    coord_cartesian(xlim = c(0, 40)) +
    theme_minimal() +
    theme(text = element_text(size = 20))
  
  #absabund_plot
  #ggsave(paste0(fig_dir, 'abund_salt_scen', scen_id, '_n50.pdf'), width = 6, height = 4)
  
  # Relative Abundance
  relabund_plot <- ggplot(results_list[[scen_id]]) +
    geom_line(aes(x = salt, y = abundance / tot_biomass, group = orig_rank),  color = 'grey', linewidth = 0.5, show.legend = F) +
    labs(x = 'Salinity (g/L)', y = 'Relative abundance') +
    coord_cartesian(xlim = c(0, 40)) +
    theme_minimal() +
    theme(text = element_text(size = 20))
  
  #relabund_plot
  #ggsave(paste0(fig_dir, 'relabund_salt_scen', scen_id, '_n50.pdf'), width = 6, height = 4)
  
  # Growthrate
  commgrowth_plot <- ggplot(results_list[[scen_id]]) +
    geom_line(aes(x = salt, y = salt_r, group = orig_rank), color = 'grey', 
              linewidth = 0.5, show.legend = F) +
    geom_line(aes(x = salt, y = mean_salt_r),
              linewidth = 3, show.legend = F) +
    geom_hline(yintercept = death, linetype = 'dashed') +
    coord_cartesian(xlim = c(0, 40)) +
    labs(x = 'Salinity (g/L)', y = 'Growth rate (1/h)') +
    theme_minimal() + 
    theme(text = element_text(size = 20))
  
  #commgrowth_plot
  #ggsave(paste0(fig_dir, 'grate_salt_scen', scen_id, '_n50.pdf'), width = 6, height = 4)
  
  relabund_plot + commgrowth_plot + absabund_plot +
    plot_layout(ncol = 1, guides = 'collect') &
    theme(legend.position = "none")
  
  ggsave(paste0(fig_dir, 'FigS10_', scen_id, '_n50.pdf'), width = 6, height = 12)
  
}


### Statistics across many random simulations - Panel F ########################
# # Scenario 1 but across 100 sets of random parameters
# # growth rates decrease linearly at same rate from 
# # a sampled rmax, until they hit 0
# 
# # iterate over 4 slopes
# # and several values for the mean of the interaction matrix
# n_species = 50
# death_val = 0.3
# 
# result_df <- data.frame()
# for (rep_i in 1:100){
# 
#   r_vec <- rnorm(n_species, mean = 0.8, sd = 0.4)
#   r_rank_orig = order(r_vec, decreasing = T)
# 
#   for (slope_val in c(-0.01, -0.02, -0.04)){
#   #for (slope in c(-0.02)){
#     for (alpha_sd in c(0.1, 0.2, 0.4, 0.5)){
#     #for (alpha_sd in c(0.5)){
# 
#       alpha_m <- matrix(runif(n_species^2, min = 0, max = 2*alpha_sd), ncol=n_species)
#       diag(alpha_m) <- 1
# 
#       for (salt_val in c(0, 5, 10, 15, 20, 25, 30, 35, 40)){
#         r_salt = pmax(r_vec + slope_val*salt_val, 0)
# 
#          new_result_df <- run_scenario(salt_val, r_vec, r_salt, alpha_m, death_val, 
#                                        t.vec = 1:5000, n_species = 50, existence_threshold = 1e-6)
# 
#          new_result_df <- new_result_df %>%
#            mutate(rep = rep_i,
#                   slope = slope_val,
#                   alpha = alpha_sd,
#                   death = death_val)
# 
#         result_df = bind_rows(result_df, new_result_df)
#       }
#     }
#   }
# }
# write_csv(result_df, '../data/2_model/n50_simulations.csv')

result_df <- read_csv('../data/2_model/n50_simulations.csv')

salt_baselines <- result_df %>%
  group_by(slope, alpha) %>%
  summarize(gr = mean(mean_r[salt == 0]))


## Plot Community composition as a function of salinity - Panel 2F
slope_val = -0.02
prior_df <- crossing(r_prior = rnorm(1000, mean = 0.6, sd = 0.3), salt = seq(0, 40, 5)) %>%
  mutate(r_eff_prior = pmax(r_prior + slope_val*salt, 0))

ggplot(result_df %>% filter(slope == slope_val, alpha == 0.1)) +
  #geom_boxplot(data = prior_df, aes(x = salt, y = r_eff_prior, group = salt), color = 'grey', width = 2) + 
  geom_boxplot(aes(x = salt, y = mean_salt_r, group = salt), color = 'black', width = 2) + 
  geom_boxplot(aes(x = salt, y = mean_r, group = salt), color = '#A72626', width = 2) + 
  geom_hline(data = salt_baselines %>% filter(slope == slope_val, alpha == 0.2), 
             aes(yintercept = gr), alpha = 0.3, color = '#A72626', linewidth = 2) + 
  geom_abline(intercept = salt_baselines %>% filter(slope == slope_val, alpha == 0.2) %>% pull(gr),
              slope = slope_val, alpha = 0.3, color = 'black', linewidth = 2) + 
  #geom_abline(intercept = 0.6, slope = slope_val, alpha = 0.3, color = 'grey', linewidth = 2) + 
  labs(y = 'Growth rate (1/h)', x = 'Salinity (g/L)') +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 2)) +
  theme_minimal() + 
  theme(text = element_text(size = 20),
        legend.position = 'bottom')

ggsave(paste0(fig_dir, 'Fig2E_', slope_val, '.pdf'), width = 6, height = 4)


# # checking diversity at each salinity
# result_df %>%
#   filter(abundance  >0) %>%
#   group_by(salt, slope, alpha, rep) %>%
#   summarize(n_species = n()) %>% 
#   ggplot() + 
#   geom_boxplot(aes(x = salt, y = n_species, group = interaction(salt, slope, alpha))) +
#   facet_grid(vars(slope), vars(alpha)) +
#   labs(x = 'Salinity (g/L)', y = '# Species')
# 
# ggsave(paste0(fig_dir, 'diversity_salinity_parallel_all.pdf'), width = 12, height = 9)


# Plot Community composition as a function of salinity - all facets #####
# ggplot(result_df) +
#   geom_boxplot(aes(x = salt, y = mean_r, group = interaction(salt, slope, alpha)), color = '#A72626', width = 2) + 
#   geom_boxplot(aes(x = salt, y = mean_salt_r, group = interaction(salt, slope, alpha)), color = 'black', width = 2) + 
#   geom_hline(data = salt_baselines, 
#              aes(yintercept = gr), alpha = 0.3, color = '#A72626', linewidth = 2) + 
#   geom_abline(data = salt_baselines, aes(intercept =  gr, slope = slope), 
#               alpha = 0.3, color = 'black', linewidth = 2) + 
#   labs(y = 'Growth rate (1/h)', x = 'Salinity (g/L)') +
#   facet_grid(vars(slope), vars(alpha)) +
#   coord_cartesian(xlim = c(0, 40)) +
#   theme_bw() + 
#   theme(text = element_text(size = 20),
#         legend.position = 'bottom')
# 
# ggsave(paste0(fig_dir, 'growthrate_salinity_parallel_all.pdf'), width = 12, height = 9)

# Plot change relative to declining species growth rate
rel_result_df <- result_df %>%
  left_join(salt_baselines, by = c('slope', 'alpha')) %>%
  mutate(rel_CCI = mean_r - gr,
         rel_CGR = mean_salt_r - gr - slope*salt)

ggplot(rel_result_df %>% mutate(slope = factor(slope, levels = c('-0.01', '-0.02', '-0.04')))) +
  geom_boxplot(aes(x = salt, y = rel_CGR, group = interaction(salt, slope, alpha), color = factor(alpha)), width = 2) + 
  labs(y = 'CGR robustness (1/h)', x = 'Salinity (g/L)', color = 'Mean interaction strength') +
  facet_wrap(vars(slope)) +
  scale_color_manual(values = c('#6E969B', '#82E2F3', '#F5BD23',  '#A72626', '#B77406',  '#356FEA')) +
  theme(text = element_text(size = 20),
        legend.position = 'bottom')

ggsave(paste0(fig_dir, 'FigS12.pdf'), width = 10, height = 5)

## Pivot wider to allow comparison between the growth rate at 20 and 40 g/L ####
# summary_df <- result_df %>%
#   filter(r_id == 1 ) %>%
#   pivot_wider(id_cols = c('alpha', 'slope', 'rep'),
#               names_from = salt, names_prefix = 'mean_r_', values_from = mean_r) %>%
#   mutate(slope = factor(slope, levels = c(-0.01, -0.02, -0.04)),
#          mean_r_orig = mean_r_0) %>%
#   pivot_longer(c('mean_r_0', 'mean_r_20', 'mean_r_40'), names_to = 'salt',
#                names_prefix = 'mean_r_') 
# 
# # Mean growth rate
# ggplot(summary_df ) +
#   geom_line(data = data.frame(x = c(0, 2), y = c(0, 2)), aes(x = x, y= y),
#             linetype = 'dashed') +
#   geom_point(aes(x = value, y = mean_r_orig, color = salt), size = 3, alpha = 0.5) + 
#   labs(y = 'Mean growth rate (1/h) at 0 g/L',
#        x = 'Mean growth rate (1/h) after salinity increase',
#        color = 'Salinity increase (g/L)') +
#   coord_cartesian(xlim = c(0,2), ylim = c(0,2)) +
#   scale_color_manual(values = c('black', '#A72626', '#6E969B')) +
#   facet_grid(vars(slope), vars(alpha)) +
#   theme_bw() + 
#   theme(text = element_text(size = 20),
#         legend.position = 'bottom')
# 
# # save the complete faceted plots with rows/cols corresponding to slope/alpha
# ggsave(paste0(fig_dir, 'growthrate_shift_slope.pdf'))

### Simulating a mixture of species with different rmax, smax ########################

# 25 species with a maximum at 0, slopes of -0.05
# 25 species with a maximum at 20, slopes of -0.1 towards 0 g/L, -0.05 towards higher salinities
# drawn from the same distribution of rmax 

n_species = 50
r_vec <- rnorm(n_species, mean = 0.8, sd = 0.4)
death = 0.3

alpha_m <- matrix(runif(n_species^2, min = 0, max = 0.5), ncol=n_species)
diag(alpha_m) <- 1

slope = -0.01

mix_result_df <- data.frame()
for (salt_val in 0:40){
  r_salt = numeric(n_species)
  r_salt[1:(n_species/2)] = pmax(r_vec[1:(n_species/2)] + slope*salt_val, 0)
  if(salt_val > 20){
    r_salt[(n_species/2+1):n_species] =   pmax(r_vec[(n_species/2+1):n_species] + 2*slope*(salt_val - 20), 0)
  } else{
    r_salt[(n_species/2+1):n_species] = pmax(r_vec[(n_species/2+1):n_species] - 0.1*(20 - salt_val), 0)
  }

  result_df <- run_scenario(salt_val, r_vec, r_salt, alpha_m, death, t.vec = 1:5000, n_species = 50, existence_threshold = 1e-6)
  
  mix_result_df = bind_rows(mix_result_df, result_df)
}

mix_result_df <- mix_result_df %>%
  mutate(group_id = ifelse(r_id <= 25, "1", "2"))

# Growthrate
CGR_mix_plot <- ggplot(mix_result_df ) +
  geom_line(aes(x = salt, y = salt_r, group = orig_rank, color = group_id), alpha = 0.5, 
            linewidth = 0.5, show.legend = F) +
  geom_line(aes(x = salt, y = mean_salt_r),
            linewidth = 3, show.legend = F) +
  geom_hline(yintercept = death, linetype = 'dashed') +
  coord_cartesian(xlim = c(0, 40)) +
  scale_color_manual(values = c("#B77406", "#6E969B")) +
  labs(x = 'Salinity (g/L)', y = 'Growth rate (1/h)') +
  theme_minimal() + 
  theme(text = element_text(size = 20))

CGR_mix_plot
ggsave(paste0(fig_dir, 'FigS13A.pdf'), width = 6, height = 4)

# Abundance of environment groups
total_mix <- mix_result_df %>%
  group_by(salt) %>%
  summarise(total_abundance = sum(abundance),
            total_mean_salt_r = max(weighted.mean(salt_r, abundance), 0, na.rm = T),
            total_mean_r = max(weighted.mean(orig_r, abundance), 0, na.rm =T) )

test_mix_df <- mix_result_df %>%
  group_by(salt, group_id) %>%
  summarise(group_abundance = sum(abundance),
            group_mean_salt_r = max(weighted.mean(salt_r, abundance), 0, na.rm = T),
            group_mean_r = max(weighted.mean(orig_r, abundance), 0, na.rm =T) ) %>%
  left_join(total_mix, by = c('salt')) %>%
  mutate(rel_group_abund = group_abundance/total_abundance)

# CCI
CCI_mix_plot <- ggplot(mix_result_df ) +
  geom_point(data = test_mix_df, aes(x = salt, y = group_mean_r, color = group_id), alpha = 0.75, size = 3) +
  geom_line(data = test_mix_df, aes(x = salt, y = total_mean_r), size = 2, color = '#A72626') +
  
  geom_hline(yintercept = death, linetype = 'dashed') +
  coord_cartesian(xlim = c(0, 40)) +
  scale_color_manual(values = c('#B77406', '#6E969B')) +
  labs(x = 'Salinity (g/L)', y = 'Community composition index', color = 'Community') +
  theme_minimal() + 
  theme(text = element_text(size = 20),
        legend.position = 'bottom')


CGR_mix_plot + CCI_mix_plot +
  plot_annotation(tag_levels = 'A')

ggsave(paste0(fig_dir, 'FigS13.pdf'), width = 12, height = 6)





