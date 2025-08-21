######################################################################################################################
## Fit growth dynamics 
## Author: JS Huisman
######################################################################################################################
library(tidyverse)
library(patchwork)
library(segmented)
library(readxl)

proc_data_dir = '../data/1.2_isolates/'

theme_set(theme_minimal() + 
            theme(text = element_text(size = 20)))

fig_dir = '../figures/1.2_isolates/'

all_results_df <- read_csv(paste0(proc_data_dir, 'all_growthdata.csv'))
average_growth_df <- read_csv(paste0(proc_data_dir, 'average_growthrates.csv'))

isolate_info <- read_csv('../data/1.2_isolates/paper_isolate_ids.csv')

color_list <- c('Estuary' = '#A72626', 
                'Marine 2' = '#F5BD23', 
                'Marine 1' = '#356FEA', 
                'Brackish' = '#004D07')

######################################################################################################################

# all_exp_ids = all_results_df %>% pull(exp_id) %>% unique()
# map_exp_ids = setNames(1:length(all_exp_ids), all_exp_ids)
# 
# fit_results_df = data.frame()
# for(id in unique(all_results_df$isolate_id)){
#   isolate_df <- all_results_df %>%
#     filter(isolate_id == id) %>%
#     dplyr::select(well_salt, growthrate = max_percap, measurement_date, plate_id, exp_id) %>%
#     mutate(exp_id_number = map_exp_ids[exp_id])
#   
#   if(dim(isolate_df %>% filter(growthrate != 0))[1] == 0){next}
#   
#   n_reps = length(unique(isolate_df$exp_id))
#   n_dates = length(unique(isolate_df$measurement_date))
#   
#   #smax = average_growth_df %>% filter(isolate_id == id, growthrate == 0) 
#   
#   # create a linear model
#   my.lm <- lm(growthrate ~ well_salt , data = isolate_df)
#   coef_linear <- round(coef(my.lm), 4)
#   sd_coef_linear <- round(summary(my.lm)$coefficients[, 2], 4)
#   lin_sigma <- round(summary(my.lm)$sigma, 4)
#   lin_adj_rsquared <- round(summary(my.lm)$adj.r.squared, 4)
#   
#   # and a segmented version of this model
#   my.seg <- try(segmented(my.lm, seg.Z = ~ well_salt, psi = 30))
#   if(any(class(my.seg) == "try-error")){next}
#   
#   coef_seg <- round(coef(my.seg), 4)
#   sd_coef_seg <- round(summary(my.seg)$coefficients[, 2], 4)
#   seg_sigma <- round(summary(my.seg)$sigma, 4)
#   seg_adj_rsquared <- round(summary(my.seg)$adj.r.squared, 4)
#   
#   ### get the slopes ####
#   # At the breakpoint the segments b and c intersect: b0 + b1*x = c0 + c1*x
#   # Important: the coefficients are the differences in slope in comparison to the previous slope
#   #Solve for c0 (intercept of second segment):
#   changepoint = round(my.seg$psi[[2]], 4)
#   slope_2 <- round(coef_seg[[2]] + coef_seg[[3]], 4)
#   intercept_2 <- round(coef_seg[[1]] + coef_seg[[2]] * changepoint - slope_2 * changepoint, 4)
#   fit_smax = round(ifelse(slope_2 == 0, - coef_seg[[1]]/coef_seg[[2]], -intercept_2/slope_2), 4)
#   fit_gchange = round(coef_seg[[1]] + coef_seg[[2]]*changepoint, 4)
#   
#   # return all these points in a dataframe structure
#   isolate_fit = data.frame(isolate_id = id, nreps = n_reps, ndates = n_dates,
#                  lin_x0 = coef_linear[[1]], lin_slope = coef_linear[[2]],
#     lin_x0_sd = sd_coef_linear[[1]], lin_slope_sd = sd_coef_linear[[2]], lin_sigma,lin_adj_rsquared,
#     seg_x0_1 = coef_seg[[1]], seg_slope1 = coef_seg[[2]], changepoint = changepoint, seg_slope2 = slope_2,
#     seg_x0_2 = intercept_2, fit_smax = fit_smax, fit_gchange = fit_gchange, seg_sigma, seg_adj_rsquared,
#     fit_sopt = changepoint)
#   
#   fit_results_df = bind_rows(fit_results_df, isolate_fit)
#   
# }
# 
# write_csv(fit_results_df, paste0(proc_data_dir, 'segment_fits.csv'))
###########################################################
fit_results_df <- read_csv(paste0(proc_data_dir, 'segment_fits.csv'))


###########################################################
# Some initial plotting to assess goodness of fit and the like

# Histogram of the Goodness of fit for segmented models
ggplot(fit_results_df %>% filter(ndates >= 3) ) +
  geom_histogram(aes(x = lin_sigma), fill = 'darkblue', alpha = 0.3, position = 'identity') +
  geom_histogram(aes(x = seg_sigma), fill = 'tomato', alpha = 0.5, position = 'identity') +
  geom_vline(xintercept = 0.3) +
  labs(x = 'Sigma of fit', y = '# Isolates')

ggplot(fit_results_df %>% filter(ndates >= 3) ) +
  geom_histogram(aes(x = lin_adj_rsquared), fill = 'darkblue', alpha = 0.3, position = 'identity') +
  geom_histogram(aes(x = seg_adj_rsquared), fill = 'tomato', alpha = 0.5, position = 'identity') +
  geom_vline(xintercept = 0.2) +
  labs(x = 'Adj. r squared of fit', y = '# Isolates')

# All negative slopes for the second segment
ggplot(fit_results_df) +
  geom_histogram(aes(x = seg_slope2), color = 'black', alpha = 0.5, bins = 25) +
  geom_histogram(aes(x = lin_slope), fill = 'darkred', color = 'darkred', alpha = 0.5, bins = 25) +
  labs(x = "Slope", y = "# Isolates", fill = "Linear fit or second segment") +
  theme(legend.position = 'bottom')

# Plot the dynamics of a specific bug
test_id = 'C1'
test_df <- fit_results_df %>% filter(isolate_id == test_id)

ggplot(average_growth_df %>% filter(isolate_id == test_id)) +
  geom_point(aes(x = well_salt, y = growthrate), size = 2) +
  geom_errorbar(aes(x = well_salt, ymin = pmax(growthrate - gr_sd, 0), ymax = growthrate + gr_sd)) + 
  geom_abline(intercept = test_df$lin_x0, slope = test_df$lin_slope, colour = "darkblue") +
  geom_vline(xintercept = 35, linetype = 'dotted') +
  geom_abline(intercept =  test_df$seg_x0_1, slope = test_df$seg_slope1, colour = "tomato") +
  geom_abline(intercept =  test_df$seg_x0_2, slope = test_df$seg_slope2, colour = "tomato4") +
  geom_vline(xintercept = fit_results_df %>% filter(isolate_id == test_id) %>% pull(fit_sopt), linetype = "dashed") 


###########################################################
#### Filtering and organisation of the results ####
# Include only isolates with at least three replicates (on different days)
# Add information on the isolates
#  if the first of the two segmented slopes is negative, the sopt is assumed at 0
# If the linear fit is better (r^2 and sigma), pick those slopes
clean_fit_df <- fit_results_df %>% 
  filter(ndates >= 3) %>%
  left_join(isolate_info, by = "isolate_id") %>%
  mutate(y0 = case_when(
              lin_sigma < seg_sigma ~ lin_x0, 
              seg_slope2 == 0 ~ seg_x0_1,
              TRUE ~ seg_x0_2),
         slope = case_when(
           lin_sigma < seg_sigma ~ lin_slope, 
           seg_slope2 == 0 ~ seg_slope1,
           TRUE ~ seg_slope2),
         sopt = case_when(
           seg_slope1 <= 0 ~ 0, 
           seg_slope2 == 0 ~ 0,
           TRUE ~ fit_sopt),
         iso_day = ifelse(is.na(plate_salinity), 'D0', 'D7')) %>%
  dplyr::select(isolate_id, iso_day, y0, slope, fit_smax, fit_gchange, sopt, inoc, class, order, family, genus, species, plate_salinity) %>%
  mutate(inoc_map = case_when(
    inoc == 'charles' ~ 'Brackish',
    inoc == 'ica' ~ 'Estuary',
    inoc == 'nahant' ~ 'Marine 1',
    inoc == 'fucus' ~ 'Marine 2'
  )) %>%
  filter(!is.na(inoc_map))

###########################################################
# Most sopt are at 0 or around 25 g/l
ggplot(fit_results_df ) +
  geom_histogram(aes(x = fit_sopt), 
                 alpha = 0.5, position = 'identity', binwidth = 5) +
  labs(x = "Salinity at segment breakpoint", y = "# Isolates") +
  theme(legend.position = 'bottom')

sopt_plot <- ggplot(clean_fit_df ) +
  geom_histogram(aes(x = sopt, color = inoc_map, fill = inoc_map), alpha = 0.5, binwidth = 5, position = 'stack') +
  labs(x = "Segment changepoint (g/L)", y = "# Isolates", fill = "Inoculum", color = "Inoculum") +
  facet_wrap(vars(inoc_map)) +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  theme(legend.position = 'bottom')

sopt_plot


# Slopes by Inoculum
slope_plot <- ggplot(clean_fit_df ) +
  geom_histogram(aes(x = slope, color = inoc_map, fill = inoc_map), alpha = 0.5, bins = 25, position = 'identity') +
  geom_freqpoly(aes(x = slope, color = inoc_map), linewidth = 1,  bins = 10, position = 'identity') +
  labs(x = "Slope", y = "# Isolates", fill = "Inoculum", color = "Inoculum") +
  #facet_wrap(vars(inoc)) +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  theme(legend.position = 'bottom')

slope_plot

## Values of mean slope and sopt
clean_fit_df %>%
  group_by(inoc) %>%
  summarize(sopt_mean = mean(sopt, na.rm = T),
            sopt_median = median(sopt, na.rm = T),
            sopt_sd = sd(sopt, na.rm = T),
            slope_mean = mean(slope, na.rm = T),
            slope_sd = sd(slope, na.rm = T))

# Combine plots for paper figure

sopt_plot + slope_plot +
  plot_annotation(tag_levels = 'A') +
  plot_layout(nrow = 1)

ggsave(paste0(fig_dir, 'FigS9_Fit_results_by_inoc.png'), width = 14, height = 6)


