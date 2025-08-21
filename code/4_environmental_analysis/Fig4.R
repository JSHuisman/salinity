###########################################################
# Environmental data analysis
# 
# Author: Martina dal Bello, adapted by JS Huisman
###########################################################

library(tidyverse)
library(ggplot2)

# Statistics
library(mgcv) # to run gam
library(performance) # check_collinearity

# To plot maps
library(rnaturalearth) # ne_countries
library(sf) # st_centroid


theme_set(theme_minimal() + theme(text = element_text(size = 20)))

proc_data_dir = '../data/4_environmental_analysis'
fig_dir = '../figures/4_environmental_analysis'

color_list <- c('chesapeake' = '#A72626', 
                'beaufort' = '#F5BD23', 
                'baltic' = '#356FEA', 
                'pico' = '#004D07')

###########################################################
salinity_tb <- read_csv(paste0(proc_data_dir, "/mcn_df_all.csv")) %>%
  mutate(dataset = ifelse(dataset == 'beauford', 'beaufort', dataset))
salinity_tb

##### Baltic Sea dataset ######################################

baltic_tb <- salinity_tb %>%
  filter(dataset == "baltic",
         !is.na(temperature))

# Linear regression, just MCN vs salinity
baltic_lm <- lm(mcn ~ salinity, data = baltic_tb)
summary(baltic_lm)

# Run and check GAM model
gam_baltic_1 <- gam (mcn ~ ammonium + nitrate_nitrite + phosphate + salinity + temperature,  data=baltic_tb) 

summary(gam_baltic_1)
check_collinearity(gam_baltic_1)

gam.check(gam_baltic_1)

# Calculate standard deviation of environmental predictors 
env_stdvar_baltic_long <- baltic_tb %>%
  select(ammonium, nitrate_nitrite, phosphate, salinity, temperature) %>%
  pivot_longer(everything(), names_to = "env_var", values_to = "value")

env_stdvar_baltic <- env_stdvar_baltic_long %>%
  filter(!is.na(value)) %>%
  group_by(env_var) %>%
  summarize(stdvar = sd(value))

env_stdvar_baltic

### Extract gam coefficients for each environmental variable and weigh by their standard deviation

env_var <- c("intercept", "ammonium", "nitrate_nitrite", "phosphate", "salinity", "temperature")
sd_env_var <- c(1, env_stdvar_baltic$stdvar)
coef_baltic <- gam_baltic_1$coefficients[1:6]*sd_env_var

coef_baltic_tb <- tibble(dataset = "baltic", env_var, coef = coef_baltic)


##### Beaufort lagoons dataset ######################################

beaufort_tb <- salinity_tb %>%
  filter(dataset == "beaufort",
         !is.na(mcn))

# Linear regression, just MCN vs salinity
beaufort_lm <- lm(mcn ~ salinity, data = beaufort_tb)
summary(beaufort_lm)

### Run and check gam model
gam_beaufort_1 <- gam (mcn ~ ammonium + dissolved_oxygen + doC + poC + salinity + tdN + temperature,  data = beaufort_tb) 

summary(gam_beaufort_1)
check_collinearity(gam_beaufort_1)

gam.check(gam_beaufort_1)

### Calculate standard deviation of environmental predictors
env_stdvar_beaufort_long <- beaufort_tb %>%
  select(ammonium, dissolved_oxygen, doC, poC, salinity, tdN, temperature) %>%
  pivot_longer(everything(), names_to = "env_var", values_to = "value")

env_stdvar_beaufort <- env_stdvar_beaufort_long %>%
  filter(!is.na(value)) %>%
  group_by(env_var) %>%
  summarize(stdvar = sd(value))

### Extract gam coefficients for each environmental variable and weigh by their standard deviation
env_var <- c("intercept", "ammonium", "dissolved oxygen", "DOC", "POC", "salinity", "tdN", "temperature")
sd_env_var <- c(1, env_stdvar_beaufort$stdvar)
coef_beaufort <- gam_beaufort_1$coefficients[1:8]*sd_env_var

coef_beaufort_tb <- tibble(dataset = "beaufort", env_var, coef = coef_beaufort)
coef_beaufort_tb

##### Chesapeake bay dataset ######################################

chesa_tb <- salinity_tb %>%
  filter(dataset == "chesapeake")

# Linear regression, just MCN vs salinity
chesa_lm <- lm(mcn ~ salinity, data = chesa_tb)
summary(chesa_lm)

# Run and check GAM 
gam_chesa_1 <- gam (mcn ~ carbon_per_liter + nitrogen_per_liter + par + salinity + temperature,  data=chesa_tb) 

summary(gam_chesa_1)
check_collinearity(gam_chesa_1)

gam.check(gam_chesa_1)

### Calculate standard deviation of environmental predictors
env_stdvar_chesa_long <- chesa_tb %>%
  select(carbon_per_liter, nitrogen_per_liter, par, salinity, temperature) %>%
  pivot_longer(everything(), names_to = "env_var", values_to = "value")

env_stdvar_chesa <- env_stdvar_chesa_long %>%
  filter(!is.na(value)) %>%
  group_by(env_var) %>%
  summarize(stdvar = sd(value))

### Extract gam coefficients for each environmental variable and weigh by their standard deviation
env_var <- c("intercept","carbon_per_liter", "nitrogen_per_liter", "par", "salinity", "temperature")
sd_env_var <- c(1, env_stdvar_chesa$stdvar)
coef_chesa <- gam_chesa_1$coefficients[1:6]*sd_env_var

coef_chesa_tb <- tibble(dataset = "chesapeake", env_var, coef = coef_chesa)
coef_chesa_tb

#####  Pivers Island (Pico) dataset ######################################
#Note: in this dataset free living and particulate fractions are mixed

pivers_tb <- salinity_tb %>%
  filter(dataset == "pico") %>%
  mutate_at(vars(starts_with("Date")), list(year=year, month=month, day=day)) %>%
  filter(!is.na(mcn))

# Linear regression, just MCN vs salinity
pivers_lm <- lm(mcn ~ salinity, data = pivers_tb)
summary(pivers_lm)

## GAM
pivers_gam_1 <- gam (mcn ~ ammonium + diC + dissolved_oxygen + phosphate + nitrate_nitrite + salinity + temperature + s(month, bs="cs"),  data=pivers_tb)

summary(pivers_gam_1)
check_collinearity(pivers_gam_1)
gam.check(pivers_gam_1)


### Calculate standard deviation of environmental predictors
env_stdvar_pivers_long <- pivers_tb %>%
  select(ammonium, diC, dissolved_oxygen, phosphate, nitrate_nitrite, salinity, temperature) %>%
  pivot_longer(everything(), names_to = "env_var", values_to = "value")

env_stdvar_pivers <- env_stdvar_pivers_long %>%
  filter(!is.na(value)) %>%
  group_by(env_var) %>%
  summarize(stdvar = sd(value))

env_stdvar_pivers


### Extract gam coefficients for each environmental variable and weigh by their standard deviation
env_var <- c("intercept","ammonium", "DIC", "dissolved oxygen", "nitrate_nitrite", "phosphate", "salinity", "temperature")
sd_env_var <- c(1, env_stdvar_pivers$stdvar)
coef_pivers <- pivers_gam_1$coefficients[1:8]*sd_env_var

coef_pivers_tb <- tibble(dataset = "pivers", env_var, coef = coef_pivers)
coef_pivers_tb

######################################################################################################################
#####  Plot all results - Fig. 4 ######################################

## Panel C: GAM Coefficients
env_var_order = c('salinity', 'temperature', 'ammonium', 'nitrate + nitrite', 
                  'd. nitrogen', 'nitrogen/L',
                  'phosphate', 'd. oxygen', 'd. organic carbon',
                  'p. organic carbon', 'd. inorganic carbon', 'carbon/L',  'light')


coef_tb <- rbind(coef_baltic_tb, coef_beaufort_tb, coef_chesa_tb, coef_pivers_tb) %>%
  filter (env_var != "intercept") %>%
  mutate(dataset = factor(dataset, levels = c('beaufort', 'baltic', 'pivers', 'chesapeake'),
                          labels = c('Beaufort Sea', 'Baltic Sea', 'Pivers Island', 'Chesapeake Bay')),
         env_var = case_when(
           env_var == 'nitrate_nitrite' ~ "nitrate + nitrite",
           env_var == 'carbon_per_liter' ~ 'carbon/L',
           env_var == 'nitrogen_per_liter' ~ 'nitrogen/L',
           env_var == 'dissolved oxygen' ~ 'd. oxygen',
           env_var == 'DOC' ~ 'd. organic carbon',
           env_var == 'POC' ~ 'p. organic carbon',
           env_var == 'DIC' ~ 'd. inorganic carbon',
           env_var == 'tdN' ~ 'd. nitrogen',
           env_var == 'par' ~ 'light',
           T ~ env_var
         )) %>%
  mutate(env_var = factor(env_var, levels = env_var_order ))


ggplot(coef_tb ) + 
  geom_tile(aes(y = dataset, x = env_var, fill= coef), color = 'gray60')+
  #scale_fill_fermenter(palette = "RdBu") +
  scale_fill_gradient2(high = '#A72626', low = '#356FEA') +
  labs(y ="", x = "", fill = "Adjusted \nGAM \ncoefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(paste0(fig_dir, '/Fig4C.pdf'), width = 12, height = 5) 


## Panel B: MCN vs salinity ####
ylims = data.frame(dataset = c('beaufort', 'baltic', 'pico', 'chesapeake'),
                   ymin = c(1.8, 0.8, 1.2, 1.5),
                   ymax = c(3.2, 4.5, 2.25, 2.8))

data_list = list('beaufort' = beaufort_tb, 'baltic' = baltic_tb,
                 'pico' = pivers_tb, 'chesapeake' = chesa_tb)

for(data_id in c('beaufort', 'baltic', 'pico', 'chesapeake')){
  
  ggplot(data_list[[data_id]]) +
    geom_point(aes(salinity, mcn, fill = dataset), color = "white" , size = 5, alpha = 0.4, shape = 21, show.legend = F) +
    geom_smooth(aes(salinity, mcn, color = dataset), method = lm, linetype = "dashed", linewidth = 2, se = FALSE, show.legend = F ) +
    scale_color_manual(values = color_list)+
    scale_fill_manual(values = color_list) +
    coord_cartesian(ylim = c(ylims %>% filter(dataset == data_id) %>% pull(ymin),
                             ylims %>% filter(dataset == data_id) %>% pull(ymax))) +
    labs(y="MCN", x = "Salinity (g/L)")
  
  ggsave(paste0(fig_dir, '/Fig4B_', data_id, '.pdf'), width = 6, height = 4)
}


## Panel A: Maps ####

world <- ne_countries(scale = "medium", returnclass = "sf")
world_points<- st_centroid(world)

# extract labels
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

map_lims = data.frame(baltic = c(0.00, 40.00, 50.00, 71.00),
                      chesapeake = c(-73.00, -79.00, 41.00, 36.00),
                      beaufort = c(-140.50, -145.00, 69.5, 70.5))

for(data_id in c('beaufort', 'baltic', 'chesapeake')){
  
  ggplot(data = world) +
    geom_sf(fill= "gray90") +
    geom_point(data = data_list[[data_id]], aes(longitude, latitude, fill = salinity), color = "white", 
               size = 10, shape = 21, alpha = 0.8)+
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(map_lims[1, data_id], map_lims[2, data_id]), 
             ylim = c(map_lims[3, data_id], map_lims[4, data_id]), expand = FALSE) +
    #scale_fill_viridis_c(option = 'viridis') +
    scale_fill_gradient(high = '#A72626', low = '#356FEA', limits = c(0, 45)) +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = "white"))+
    theme(legend.position="none",
          axis.text = element_blank(),
          axis.title = element_blank())

  ggsave(paste0(fig_dir, '/Fig4A_', data_id, '.png'))
}


### Pivers Island timeseries
ggplot(pivers_tb %>%
         mutate(date = as.Date(date)), aes(x=date, y = salinity))+
  geom_smooth(method = "lm", formula = y ~ poly(x, 10), size = 0.5, color = "gray60") +
  geom_point(aes(fill = salinity), color = "white", size=3, shape =21, alpha = 0.8)+
  theme(legend.position="bottom", legend.key.width = unit(1, "cm"))+
  scale_fill_gradient(high = '#A72626', low = '#356FEA', limits = c(0, 45)) +
  labs(x="Time", y = "Salinity (g/L)", fill = "Salinity (g/L)") 

ggsave(paste0(fig_dir, '/Fig4A_pivers_timeseries.pdf'), width = 6, height = 4)
