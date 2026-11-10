###...................................###
# Author: Georgia G. Hernandez
# Date: 2026-11-05
# R code for reproducing analyses and figures
# Paper: Photosynthetic heat tolerance is shaped by trait coordination, acclimation, and evolutionary history
###...................................###

# Load packages 

library(tidyverse)
library(readxl)
library(ggcorrplot)
library(factoextra)

# Bayesian modeling and diagnostics
library(brms)
library(bayesplot)
library(posterior)
library(modelr)
library(tidybayes)

# Visualization
library(viridis)
library(patchwork)
library(cowplot)
library(ggh4x)
library(gridExtra)
library(ggrepel)

# Spatial data and mapping
library(measurements)
library(sf)
library(geodata)
library(terra)
library(stars)
library(usmap)
library(ggspatial)

# Phylogeny
library(ape)
library(ggtree)
library(caper)

# Import data

# Heat tolerance + traits + general info
data_cali <- read_excel("data_california.xlsx")

# Read tree 
tree <- read.tree("figure_1a.tre")

# Read weather data 
spp_wd <- read_excel("sum_temp_prism.xlsx")
prev_wd <- read_csv("temp_sum_prevdays.csv")

# ------------------------------------------------------------------------------
# Data handling 

# Weather data -----------------------------------------------------------------
# Weather data, summarized by years
# Match sites names of the weather sites to the heat tolerance data
spp_wd <- spp_wd |> 
  dplyr::select(!c(...1, TMAX_sd_years1981_2024, TMIN_sd_years1981_2024, 
                   MAT_sd_years1981_2024, TMAX_sd_years2015_2024,
                   TMIN_sd_years2015_2024, MAT_sd_years2015_2024, 
                   TMAX_sd_year2024, TMIN_sd_year2024, MAT_sd_year2024)) |> 
  filter(!site == "va") |> 
  mutate(site = ifelse(site == "snarl", "va", site)) # Change code of snarl to va

# Merge weather data with heat tolerance data 
data_cali <- left_join(data_cali, spp_wd, by = "site")

# Weather data, summarized by previous days
prev_wd <- prev_wd |> 
  dplyr::select(!c(year, month, mean_temp, max_temp, min_temp)) |> 
  filter(!site == "va") |> 
  mutate(site = ifelse(site == "snarl", "va", site)) 

# Merge weather data with heat tolerance data 
data_cali <- left_join(data_cali, prev_wd, by = c("site", "date"))

# Heat tolerance data ----------------------------------------------------------
# Exclude species that are non-native 
data_cali <- data_cali |> 
  filter(!origin %in% c("non-native", "non_native")) |> 
  dplyr::select(!c(origin, ind))

# Filter NAs to use only for T50
t50 <- data_cali |> 
  group_by(unique_id) |> 
  mutate(any_na = any(is.na(T50))) |> 
  ungroup() |> 
  filter(!any_na) |> 
  dplyr::select(-any_na)

# Change tau's name
t50 <- t50 |> 
  rename(TAU = tau,
         Elevation = elevation)

# Coordinates ------------------------------------------------------------------
# Convert coordinates from (degrees decimal minutes, DDM) to a Decimal Degrees
# (DD), suitable for spatial analysis

#Latitude
t50$coordinates_n <- measurements::conv_unit(t50$coordinates_n, 
                                             from = 'deg_dec_min', 
                                             to = 'dec_deg') 

# Longitude
t50$coordinates_w <- measurements::conv_unit(t50$coordinates_w, 
                                             from = 'deg_dec_min', 
                                             to = 'dec_deg') 

# Add the negative sign to longitude to state it is from the West
t50 <- t50 |> 
  mutate(neg = "-") |> 
  unite(coordinates_longitude, 
        c(neg, coordinates_w), 
        sep = "") |> 
  rename(coordinates_latitude = coordinates_n) 

# Make coordinates numeric 
t50$coordinates_latitude <- as.numeric(as.character(t50$coordinates_latitude))
t50$coordinates_longitude <- as.numeric(as.character(t50$coordinates_longitude))

# Isotopes ---------------------------------------------------------------------
# Turn the content of nitrogen and carbon into ratios
t50 <- t50 |> 
  mutate(d13C = (d13C/SampleWeight_mg),
         totalC = (totalC/SampleWeight_mg)/1000*100,
         d15N = (d15N/SampleWeight_mg),
         totalN = (totalN/SampleWeight_mg)/1000*100)

# Scale variables --------------------------------------------------------------
t50 <- t50 |> 
  mutate(tau = as.numeric(scale(log(TAU))),
         la = as.numeric(scale(log(leaf_area))),
         lw = as.numeric(scale(log(leaf_width))),
         lt = as.numeric(scale(log(leaf_thickness))),
         sla = as.numeric(scale(SLA)),
         lma = as.numeric(scale(log(LMA))),
         ldmc = as.numeric(scale(LDMC)),
         lwc = as.numeric(scale(LWC)),
         rwc = as.numeric(scale(RWC)),
         ewt = as.numeric(scale(log(EWT))), 
         elevation = as.numeric(scale(Elevation, scale = FALSE)) / 500, 
         TMAX_1981_2024 = as.numeric(scale(TMAX_years1981_2024)),
         TMIN_1981_2024 = as.numeric(scale(TMIN_years1981_2024)),
         MAT_1981_2024 = as.numeric(scale(MAT_years1981_2024)), 
         TMAX_2015_2024 = as.numeric(scale(TMAX_years2015_2024)), 
         TMIN_2015_2024 = as.numeric(scale(TMIN_years2015_2024)),
         MAT_2015_2024 = as.numeric(scale(MAT_years2015_2024)), 
         TMAX_2024 = as.numeric(scale(TMAX_year2024)),      
         TMIN_2024 = as.numeric(scale(TMIN_year2024)),       
         MAT_2024 = as.numeric(scale(MAT_year2024)),
         dc = as.numeric(scale(d13C)),
         t_c = as.numeric(scale(totalC)),
         dn = as.numeric(scale(d15N)),
         t_n = as.numeric(scale(totalN)))

# To change color for bayesian plots -------------------------------------------
orange_scheme <- c("#ffebcc", "#ffcc80", "#ffad33", 
                   "#e68a00", "#995c00", "#663d00")
color_scheme_set(orange_scheme)

# Construct a covariance matrix of species for all models ----------------------
A <- ape::vcv.phylo(tree)

# Turning covariance matrix into a correlation matrix (to make results more interpretative)
A_c <- diag(1 / sqrt(diag(A)),
            nrow = NROW(A),
            ncol = NCOL(A)) %*%
  A %*% 
  diag(1 / sqrt(diag(A)),
       nrow = NROW(A),
       ncol = NCOL(A))

rownames(A_c) <- rownames(A)
colnames(A_c) <- colnames(A)

# Summarize data for the phylogenetic analyses
t50_evo <- data_cali |> 
  group_by(unique_id) |> 
  mutate(any_na = any(is.na(T50))) |> 
  ungroup() |> 
  filter(!any_na) |> 
  dplyr::select(-any_na)

# Make summary of T50 by species for phylogenetic signal for all traits
traits_t50_all <- t50_evo |> 
  group_by(species_3) |> 
  summarise(t50 = mean(T50), 
            t15 = mean(T15),
            tau = mean(tau),
            LMA = mean(LMA),
            LDMC = mean(LDMC), 
            LWC = mean(LWC), 
            RWC = mean(RWC),
            EWT = mean(EWT), 
            lt = mean(leaf_thickness),
            la = mean(leaf_area),
            lw = mean(leaf_width),
            TMAX = mean(TMAX_year2024),
            MAT = mean(MAT_year2024),
            TMAX3 = mean(maxt_3day),
            MAT3 = mean(maxt_3day)) |> 
  ungroup() 

# I have to divide the traits so these 2 that have multiple NAs don't exclude species 
traits2 <- data_cali |> 
  group_by(species_3) |> 
  summarise(d13 = mean(d13C),
            tn = mean(totalN), 
            TMAX = mean(TMAX_year2024),
            MAT = mean(MAT_year2024)) |> 
  ungroup() |> 
  drop_na()

############################## Research Questions ##############################
# Q1. Are heat tolerance and other leaf traits constrained by evolutionary history? 

# Q2. How does heat tolerance vary in response to short (1 to 15-day) and long-term 
# temperature variation, and how is this response mediated by elevation?

# Q3. Which leaf traits (LMA, LDMC, tau, LWC, RWC, EWT) and plant life history 
# strategies are the strongest predictors of heat tolerance? 

# Q4.	Do traits related to carbon economics, energy balance, and heat tolerance
# show coordinated variation?

#  Analyses --------------------------------------------------------------------

# Q1. --------------------------------------------------------------------------
# Q1. Are heat tolerance and other leaf traits constrained by evolutionary history?

# Prevent a common caper bug by ensuring node labels are NULL
tree$node.label <- NULL

# Check order of tips in the tree and that the amount of species is the correct
tree$tip.label

# Transform trait data into a dataframe 
traits_t50_all <- as.data.frame(traits_t50_all)
trait_d13 <- as.data.frame(traits2)

# Combine tree and dataset into a comparative data object
comp_data_t50 <- comparative.data(phy = tree, 
                                  data = traits_t50_all, 
                                  names.col = "species_3", 
                                  vcv = TRUE)

comp_data_d13 <- comparative.data(phy = tree, 
                                  data = trait_d13, 
                                  names.col = "species_3", 
                                  vcv = TRUE)

# Phylogeny analyses for all the traits (Climate-controlled)
pgls_t50     <- pgls(t50 ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_tcrit   <- pgls(t15 ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_lma   <- pgls(LMA ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_ldmc  <- pgls(LDMC ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_lwc   <- pgls(LWC ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_rwc   <- pgls(RWC ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_ewt   <- pgls(EWT ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_lt    <- pgls(lt ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_la    <- pgls(la ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_lw    <- pgls(lw ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_tau <- pgls(tau ~ TMAX, data = comp_data_t50, lambda = "ML")
pgls_d13 <- pgls(d13 ~ TMAX, data = comp_data_d13, lambda = "ML")
pgls_tn  <- pgls(tn ~ TMAX, data = comp_data_d13, lambda = "ML")

# Phylogenetic signal
summary(pgls_t50)
summary(pgls_tcrit)
summary(pgls_lma)
summary(pgls_ldmc)
summary(pgls_lwc)
summary(pgls_rwc)
summary(pgls_ewt)
summary(pgls_lt)
summary(pgls_la)
summary(pgls_lw)
summary(pgls_tau)
summary(pgls_d13)
summary(pgls_tn)

# Q2. --------------------------------------------------------------------------
# Q2. How does heat tolerance vary in response to short (1 to 15-day) and long-term 
# temperature variation, and how is this response mediated by elevation?

### Short term temperature data (1,3, 15 days prior to collection) models
# TCrit 

# Best model: Max temp from the last 15 days x Mean temp from the last 15 days
mod_temp2_x15m15 <- brm(T15 ~
                          (1|gr(species_3, cov = A_c)) + 
                          s(coordinates_longitude, coordinates_latitude) +
                          elevation + 
                          s(maxt_15day, meant_15day), 
                        data = t50,
                        family = gaussian(),
                        data2 = list(A_c = A_c),
                        cores = 4,
                        iter = 9000,
                        control = list(adapt_delta = 0.89,
                                       max_treedepth = 10),
                        seed = 2025)
summary(mod_temp2_x15m15) 

# To asses convergence of MCMC
mod_temp2_x15m15_MCMC <- mcmc_combo(mod_temp2_x15m15, 
                                    combo = c("dens_overlay","trace"),
                                    regex_pars = "(^sd|^Intercept)")

# Second best model: Max temp from the last day x Mean temp from the last 15 days 
mod_temp2_x1m15 <- brm(T15 ~
                         (1|gr(species_3, cov = A_c)) + 
                         s(coordinates_longitude, coordinates_latitude) +
                         elevation + 
                         s(maxt_1day, meant_15day), 
                       data = t50,
                       family = gaussian(),
                       data2 = list(A_c = A_c),
                       cores = 4,
                       iter = 9000,
                       control = list(adapt_delta = 0.89,
                                      max_treedepth = 13),
                       seed = 2025)
summary(mod_temp2_x1m15)

mod_temp2_x1m15_MCMC <- mcmc_combo(mod_temp2_x1m15, 
                                   combo = c("dens_overlay","trace"),
                                   regex_pars = "(^sd|^Intercept)")

# T50
# Best model: Max temp from the last 15 days x Mean temp from the last 15 days
mod_temp_x15m15 <- brm(T50 ~
                         (1|gr(species_3, cov = A_c)) + 
                         s(coordinates_longitude, coordinates_latitude) +
                         elevation + 
                         s(maxt_15day, meant_15day), 
                       data = t50,
                       family = gaussian(),
                       data2 = list(A_c = A_c),
                       cores = 4,
                       iter = 9000,
                       control = list(adapt_delta = 0.89,
                                      max_treedepth = 10),
                       seed = 2025)

summary(mod_temp_x15m15)

# To asses convergence of MCMC
mod_temp_x15m15_MCMC <- mcmc_combo(mod_temp_x15m15, 
                                   combo = c("dens_overlay","trace"),
                                   regex_pars = "(^sd|^Intercept)")

# Second best model: Max temp from the last day x Mean temp from the last 15 days 
mod_temp_x1m15 <- brm(T50 ~
                        (1|gr(species_3, cov = A_c)) + 
                        s(coordinates_longitude, coordinates_latitude) +
                        elevation + 
                        s(maxt_1day, meant_15day), 
                      data = t50,
                      family = gaussian(),
                      data2 = list(A_c = A_c),
                      cores = 4,
                      iter = 9000,
                      control = list(adapt_delta = 0.89,
                                     max_treedepth = 13),
                      seed = 2025)
summary(mod_temp_x1m15)

mod_temp_x1m15_MCMC <- mcmc_combo(mod_temp_x1m15, 
                                  combo = c("dens_overlay","trace"),
                                  regex_pars = "(^sd|^Intercept)")

### Long term temperature data (1984-2024) models
# T50
mod_temp30_t50 <- brm(T50 ~
                        (1|gr(species_3, cov = A_c)) + 
                        s(coordinates_longitude, coordinates_latitude) +
                        elevation + MAT_1981_2024 + TMAX_1981_2024 + TMIN_1981_2024, 
                      data = t50,
                      family = gaussian(),
                      data2 = list(A_c = A_c),
                      cores = 4,
                      iter = 9000,
                      control = list(adapt_delta = 0.89,
                                     max_treedepth = 13),
                      seed = 2025)
summary(mod_temp30_t50)

# To asses convergence of MCMC
temp30_t50_MCMC <- mcmc_combo(mod_temp30_t50, 
                              combo = c("dens_overlay","trace"),
                              regex_pars = "(^sd|^Intercept)")

# TCrit
mod_temp30_t15 <- brm(T15 ~
                        (1|gr(species_3, cov = A_c)) + 
                        s(coordinates_longitude, coordinates_latitude) +
                        elevation + MAT_1981_2024 + TMAX_1981_2024 + TMIN_1981_2024, 
                      data = t50,
                      family = gaussian(),
                      data2 = list(A_c = A_c),
                      cores = 4,
                      iter = 5000,
                      control = list(adapt_delta = 0.8,
                                     max_treedepth = 13),
                      seed = 2025)
summary(mod_temp30_t15)

temp30_t15_MCMC <- mcmc_combo(mod_temp30_t15, 
                              combo = c("dens_overlay","trace"),
                              regex_pars = "(^sd|^Intercept)")

# Q3. --------------------------------------------------------------------------
# Q3. Which leaf traits (tau, LMA, LDMC, LWC, RWC, EWT) and plant life history 
# strategies are the strongest predictors of heat tolerance? 

# TCrit ----
# tau
mod_t15_tau <- brm(T15 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     tau, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 5000,
                   control = list(adapt_delta = 0.8),
                   seed = 2025)
summary(mod_t15_tau)

# To asses convergence of MCMC
mcmc_combo(mod_t15_tau, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LMA
mod_t15_lma <- brm(T15 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     lma, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 5000,
                   control = list(adapt_delta = 0.85),
                   seed = 2025)
summary(mod_t15_lma)

# To asses convergence of MCMC
mcmc_combo(mod_t15_lma, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LDMC
mod_t15_ldmc <- brm(T15 ~
                      (1|gr(species_3, cov = A_c)) + 
                      s(coordinates_longitude, coordinates_latitude) +
                      ldmc, 
                    data = t50,
                    family = gaussian(),
                    data2 = list(A_c = A_c),
                    cores = 4,
                    iter = 8000,
                    control = list(adapt_delta = 0.85),
                    seed = 2025)
summary(mod_t15_ldmc)

# To asses convergence of MCMC
mcmc_combo(mod_t15_ldmc, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LWC
mod_t15_lwc <- brm(T15 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     lwc, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.88),
                   seed = 2025)
summary(mod_t15_lwc)

# To asses convergence of MCMC
mcmc_combo(mod_t15_lwc, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# RWC
mod_t15_rwc <- brm(T15 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     rwc, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 5000,
                   control = list(adapt_delta = 0.8),
                   seed = 2025)
summary(mod_t15_rwc)

# To asses convergence of MCMC
mcmc_combo(mod_t15_rwc, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# EWT
mod_t15_ewt <- brm(T15 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     ewt, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 8000,
                   control = list(adapt_delta = 0.85),
                   seed = 2025)
summary(mod_t15_ewt)

# To asses convergence of MCMC
mcmc_combo(mod_t15_ewt, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LA
mod_t15_la <- brm(T15 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    la, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 8000,
                  control = list(adapt_delta = 0.86),
                  seed = 2025)
summary(mod_t15_la)

# To asses convergence of MCMC
mcmc_combo(mod_t15_la, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LW
mod_t15_lw <- brm(T15 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    lw, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 8000,
                  control = list(adapt_delta = 0.86),
                  seed = 2025)
summary(mod_t15_lw)

# To asses convergence of MCMC
mcmc_combo(mod_t15_lw, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Leaf thickness
mod_t15_lt <- brm(T15 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    lt, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 8000,
                  control = list(adapt_delta = 0.86),
                  seed = 2025)
summary(mod_t15_lt)

# To asses convergence of MCMC
mcmc_combo(mod_t15_lt, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# d13C
mod_t15_d13C <- brm(T15 ~
                      (1|gr(species_3, cov = A_c)) + 
                      s(coordinates_longitude, coordinates_latitude) +
                      dc, 
                    data = t50,
                    family = gaussian(),
                    data2 = list(A_c = A_c),
                    cores = 4,
                    iter = 8000,
                    control = list(adapt_delta = 0.85),
                    seed = 2025)
summary(mod_t15_d13C)

# total N
mod_t15_tn <- brm(T15 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    t_n, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 8000,
                  control = list(adapt_delta = 0.85),
                  seed = 2025)
summary(mod_t15_tn)

# Vegetation type
mod_t15_veg <- brm(T15 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     (1| veg_type), 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t15_veg) 

# To asses convergence of MCMC
mcmc_combo(mod_t15_veg, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Leaf habit
mod_t15_lh <- brm(T15 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    leaf_habit_spec, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 8000,
                  control = list(adapt_delta = 0.85),
                  seed = 2025)
summary(mod_t15_lh) 

# To asses convergence of MCMC
mcmc_combo(mod_t15_lh, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Seed plant
mod_t15_sp <- brm(T15 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    seed_plant, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 8000,
                  control = list(adapt_delta = 0.85),
                  seed = 2025)
summary(mod_t15_sp) 

# To asses convergence of MCMC
mcmc_combo(mod_t15_sp, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# T50 ----
# tau
mod_t50_tau <- brm(T50 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     tau, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t50_tau)

# To asses convergence of MCMC
mcmc_combo(mod_t50_tau, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LMA
mod_t50_lma <- brm(T50 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     lma, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t50_lma)

# To asses convergence of MCMC
mcmc_combo(mod_t50_lma, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LDMC
mod_t50_ldmc <- brm(T50 ~
                      (1|gr(species_3, cov = A_c)) + 
                      s(coordinates_longitude, coordinates_latitude) +
                      ldmc, 
                    data = t50,
                    family = gaussian(),
                    data2 = list(A_c = A_c),
                    cores = 4,
                    iter = 9000,
                    control = list(adapt_delta = 0.89),
                    seed = 2025)
summary(mod_t50_ldmc)

# To asses convergence of MCMC
mcmc_combo(mod_t50_ldmc, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LWC
mod_t50_lwc <- brm(T50 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     lwc, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t50_lwc)

# To asses convergence of MCMC
mcmc_combo(mod_t50_lwc, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# RWC
mod_t50_rwc <- brm(T50 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     rwc, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t50_rwc)

# To asses convergence of MCMC
mcmc_combo(mod_t50_rwc, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# EWT
mod_t50_ewt <- brm(T50 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     ewt, 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t50_ewt)

# To asses convergence of MCMC
mcmc_combo(mod_t50_ewt, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LA
mod_t50_la <- brm(T50 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    la, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 9000,
                  control = list(adapt_delta = 0.89),
                  seed = 2025)
summary(mod_t50_la)

# To asses convergence of MCMC
mcmc_combo(mod_t50_la, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# LW
mod_t50_lw <- brm(T50 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    lw, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 9000,
                  control = list(adapt_delta = 0.89),
                  seed = 2025)
summary(mod_t50_lw)

# To asses convergence of MCMC
mcmc_combo(mod_t50_lw, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Leaf thickness
mod_t50_lt <- brm(T50 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    lt, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 9000,
                  control = list(adapt_delta = 0.89),
                  seed = 2025)
summary(mod_t50_lt)

# To asses convergence of MCMC
mcmc_combo(mod_t50_lt, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# d13C
mod_t50_d13C <- brm(T50 ~
                      (1|gr(species_3, cov = A_c)) + 
                      s(coordinates_longitude, coordinates_latitude) +
                      dc, 
                    data = t50,
                    family = gaussian(),
                    data2 = list(A_c = A_c),
                    cores = 4,
                    iter = 9000,
                    control = list(adapt_delta = 0.89),
                    seed = 2025)
summary(mod_t50_d13C)

# To asses convergence of MCMC
mcmc_combo(mod_t50_d13C, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# total N
mod_t50_tn <- brm(T50 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    t_n, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 9000,
                  control = list(adapt_delta = 0.89),
                  seed = 2025)

summary(mod_t50_tn)

# To asses convergence of MCMC
mcmc_combo(mod_t50_tn, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Vegetation type
mod_t50_veg <- brm(T50 ~
                     (1|gr(species_3, cov = A_c)) + 
                     s(coordinates_longitude, coordinates_latitude) +
                     (1| veg_type), 
                   data = t50,
                   family = gaussian(),
                   data2 = list(A_c = A_c),
                   cores = 4,
                   iter = 9000,
                   control = list(adapt_delta = 0.89),
                   seed = 2025)
summary(mod_t50_veg) 
ranef(mod_t50_veg)$veg_type

# To asses convergence of MCMC
mcmc_combo(mod_t50_veg, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Leaf habit
mod_t50_lh <- brm(T50 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    leaf_habit_spec, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 9000,
                  control = list(adapt_delta = 0.89),
                  seed = 2025)
summary(mod_t50_lh) 

# To asses convergence of MCMC
mcmc_combo(mod_t50_lh, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Seed plant
mod_t50_sp <- brm(T50 ~
                    (1|gr(species_3, cov = A_c)) + 
                    s(coordinates_longitude, coordinates_latitude) +
                    seed_plant, 
                  data = t50,
                  family = gaussian(),
                  data2 = list(A_c = A_c),
                  cores = 4,
                  iter = 9000,
                  control = list(adapt_delta = 0.89),
                  seed = 2025)
summary(mod_t50_sp) 

# To asses convergence of MCMC
mcmc_combo(mod_t50_sp, 
           combo = c("dens_overlay","trace"),
           regex_pars = "(^sd|^Intercept)")

# Q4.	--------------------------------------------------------------------------
# Q4.	Do traits related to carbon economics, energy balance, and heat tolerance
# show coordinated variation?

# Fit model T50
mod_ceeb_t50 <- brm(T50 ~
                      (1|gr(species_3, cov = A_c)) +
                      s(coordinates_longitude, coordinates_latitude) +
                      tau + lma + ldmc + la + lw,
                    data = t50,
                    family = gaussian(),
                    data2 = list(A_c = A_c),
                    cores = 4,
                    iter = 9000, 
                    seed = 2025,
                    control = list(adapt_delta = 0.89))

# Summary of the model
summary(mod_ceeb_t50)

# Fit model T15
mod_ceeb_t15 <- brm(T15 ~
                      (1|gr(species_3, cov = A_c)) + 
                      s(coordinates_longitude, coordinates_latitude) +
                      tau + lma + ldmc + la + lw,
                    data = t50,
                    family = gaussian(),
                    data2 = list(A_c = A_c),
                    cores = 4,
                    iter = 9000, 
                    seed = 2025,
                    control = list(adapt_delta = 0.89))

# Summary of the model
summary(mod_ceeb_t15)

# Figures ----------------------------------------------------------------------

# Figure 1 ----

# Map of California 

# Create this data frame with actual site names and desired display names
site_name_map <- data.frame(
  stringsAsFactors = FALSE, # Good practice
  original_site_name = c("ac", "ab", "bdc", "cc", "ha", "jsj", "la", "ml", 
                         "ov", "pc", "sr", "sg", "va", "yo"), 
  display_name       = c("Angelo Coast", "Anza Borrego", "Boyd Deep Canyon",
                         "Crooked Creek", "Hastings", "James San Jacinto", 
                         "Lassen", "McLaughlin", "Owens Valley", "Putah Creek",
                         "Stunt Ranch", "Sweeney Granite", "SNARL", "Yosemite"))

# Select desired variables to work with
geo_data <- t50 |> 
  dplyr::select(site, coordinates_latitude, coordinates_longitude) |> 
  distinct(site, .keep_all = TRUE) |> 
  as.data.frame()

# Convert coordinates to an sf object
points <- geo_data |> 
  st_as_sf(coords = c("coordinates_longitude",
                      "coordinates_latitude"), 
           crs = 4326, 
           remove = FALSE) 

# Join the name map with data points
points_for_labels <- points |> 
  left_join(site_name_map, 
            by = c("site" = "original_site_name")) |> 
  mutate(map_display_label = ifelse(is.na(display_name), 
                                    site, 
                                    display_name))

# Get state map geometry
cali <- us_map(regions = "counties",
               include = "CA")

states <- us_map(regions = "states",
                 include = c("CA"))

# Add a layer with temperature to the map using BioClim: mean maximum daily 
# temperature of the warmest month:
# Download WorldClim data
bio <- worldclim_global(var = "bio", 
                        res = 2.5, 
                        path = "data/worldclim")

# Change names of the climatic variables 
names(bio) <- paste0("bio", 1:19)
plot(bio[[5]])

# Select only bio5 (Max Temperature of Warmest Month)
bio <- bio$bio5

# Transform map of California state to be able to crop bioclim map
states_trans <- st_transform(states, crs = crs(bio))

# Crop BioClim data
cali_cropped <- crop(bio, states_trans)
plot(cali_cropped)

cali_masked <- mask(cali_cropped, states_trans)

# Transform cropped map
cali_final_stars <- cali_masked |> 
  st_as_stars() |> 
  st_transform(st_crs(states))

# Plot
fig1 <- ggplot() +
  geom_stars(data = cali_final_stars, color = NA) +
  geom_sf(data = states, fill = NA, color = "black", linewidth = 0.5) +
  geom_sf(data = cali, fill = NA, color = "gray18") +
  geom_sf(data = points, color = "black", size = 2.5, alpha = 1) +

  geom_text_repel(data = points_for_labels, 
                  aes(label = map_display_label, geometry = geometry), 
                  stat = "sf_coordinates",
                  size = 4, 
                  family = "sans",
                  box.padding = 0.4,       
                  point.padding = 0.5,     
                  min.segment.length = 0,
                  segment.color = "grey50",
                  seed = 2026) +
  scale_fill_viridis_c(option = "B", 
                       name = "Max Temperature (°C)") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "gray92"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 10, family = "sans"), 
        legend.title = element_text(size = 14, face = "bold", family = "sans"),
        axis.text = element_text(size = 12)) + 
  ggspatial::annotation_scale(location = "tr", 
                              bar_cols = c("grey60", "white"),
                              text_cex = 1) + 
  ggspatial::annotation_north_arrow(location = "tr", 
                                    which_north = "true",
                                    pad_x = unit(0.9, "cm"), 
                                    pad_y = unit(0.9, "cm"),
                                    style = ggspatial::north_arrow_nautical(
                                      text_size = 12,
                                      fill = c("grey40", "white"),
                                      line_col = "grey20"))

fig1 

# Figure 2 ----

# New species names to be showed in all trees
name_map_df <- data.frame(
  stringsAsFactors = FALSE,
  original_label = c(
    "Bebbia_juncea", "Ambrosia_dumosa", "Encelia_actoni",              
    "Encelia_farinosa", "Hazardia_squarrosa", "Ericameria_nauseosa",         
    "Ericameria_discoidea", "Baccharis_pilularis", "Artemisia_tridentata",        
    "Artemisia_californica", "Senecio_flaccidus", "Tetradymia_stenolepis",       
    "Lonicera_subspicata", "Lonicera_interrupta", "Monardella_odoratissima",     
    "Salvia_dorrii", "Salvia_leucophylla", "Condea_emoryi",    
    "Justicia_californica", "Chilopsis_linearis", "Fraxinus_latifolia",          
    "Eriodictyon_californicum", "Lycium_cooperi", "Rhododendron_occidentale",    
    "Arctostaphylos_glauca", "Arctostaphylos_viscida", "Arctostaphylos_manzanita",    
    "Arctostaphylos_patula", "Arctostaphylos_pringlei", "Arctostaphylos_glandulosa",   
    "Arbutus_menziesii", "Fouquieria_splendens", "Cornus_nuttallii",  
    "Atriplex_torreyi", "Simmondsia_chinensis", "Eriogonum_umbellatum",        
    "Eriogonum_fasciculatum", "Neltuma_odorata", "Acmispon_glaber",       
    "Olneya_tesota", "Psorothamnus_schottii", "Lupinus_albifrons",           
    "Senegalia_greggii", "Parkinsonia_florida", "Rosa_woodsii",          
    "Rosa_californica", "Rubus_parviflorus", "Cercocarpus_betuloides",      
    "Purshia_tridentata", "Heteromeles_arbutifolia", "Holodiscus_discolor",     
    "Chamaebatiaria_millefolium", "Ceanothus_integerrimus", "Ceanothus_velutinus",         
    "Ceanothus_oliganthus", "Ceanothus_incanus", "Ceanothus_foliosus",          
    "Ceanothus_cordulatus", "Rhamnus_ilicifolia", "Frangula_californica",        
    "Corylus_cornuta", "Betula_occidentalis", "Alnus_incana",      
    "Alnus_rhombifolia", "Juglans_hindsii", "Notholithocarpus_densiflorus",
    "Quercus_durata", "Quercus_douglasii", "Quercus_garryana",
    "Quercus_lobata", "Quercus_berberidifolia", "Quercus_kelloggii",           
    "Quercus_agrifolia", "Quercus_wislizeni", "Quercus_chrysolepis",         
    "Salix_gooddingii", "Salix_laevigata", "Salix_exigua",        
    "Populus_trichocarpa", "Populus_fremontii", "Populus_tremuloides",         
    "Larrea_tridentata", "Cleomella_arborea", "Acer_macrophyllum",        
    "Aesculus_californica", "Malosma_laurina", "Vitis_californica",           
    "Ribes_sericeum", "Ribes_cereum", "Platanus_racemosa",           
    "Hesperoyucca_whipplei", "Yucca_schidigera", "Washingtonia_filifera",       
    "Umbellularia_californica", "Hesperocyparis_macnabiana", "Hesperocyparis_sargentii",    
    "Calocedrus_decurrens", "Sequoiadendron_giganteum", "Sequoia_sempervirens",  
    "Pinus_ponderosa", "Pinus_sabiniana", "Pinus_coulteri",       
    "Pinus_flexilis", "Pinus_lambertiana", "Pinus_longaeva",              
    "Pseudotsuga_menziesii", "Abies_concolor", "Abies_bracteata",             
    "Abies_magnifica"), 
  display_label = c(
    "Bebbia juncea", "Ambrosia dumosa", "Encelia actoni",              
    "Encelia farinosa", "Hazardia squarrosa", "Ericameria nauseosa",         
    "Ericameria discoidea", "Baccharis pilularis", "Artemisia tridentata",        
    "Artemisia californica", "Senecio flaccidus", "Tetradymia stenolepis",       
    "Lonicera subspicata", "Lonicera interrupta", "Monardella odoratissima",     
    "Salvia dorrii", "Salvia leucophylla", "Condea emoryi",    
    "Justicia californica", "Chilopsis linearis", "Fraxinus latifolia",          
    "Eriodictyon californicum", "Lycium cooperi", "Rhododendron occidentale",    
    "Arctostaphylos glauca", "Arctostaphylos viscida", "Arctostaphylos manzanita",    
    "Arctostaphylos patula", "Arctostaphylos pringlei", "Arctostaphylos glandulosa",   
    "Arbutus menziesii", "Fouquieria splendens", "Cornus nuttallii",  
    "Atriplex torreyi", "Simmondsia chinensis", "Eriogonum umbellatum",        
    "Eriogonum fasciculatum", "Neltuma odorata", "Acmispon glaber",       
    "Olneya tesota", "Psorothamnus schottii", "Lupinus albifrons",           
    "Senegalia greggii", "Parkinsonia florida", "Rosa woodsii",          
    "Rosa californica", "Rubus parviflorus", "Cercocarpus betuloides",      
    "Purshia tridentata", "Heteromeles arbutifolia", "Holodiscus discolor",     
    "Chamaebatiaria millefolium", "Ceanothus integerrimus", "Ceanothus velutinus",         
    "Ceanothus oliganthus", "Ceanothus incanus", "Ceanothus foliosus",          
    "Ceanothus cordulatus", "Rhamnus ilicifolia", "Frangula californica",        
    "Corylus cornuta", "Betula occidentalis", "Alnus incana",      
    "Alnus rhombifolia", "Juglans hindsii", "Notholithocarpus densiflorus",
    "Quercus durata", "Quercus douglasii", "Quercus garryana",
    "Quercus lobata", "Quercus berberidifolia", "Quercus kelloggii",           
    "Quercus agrifolia", "Quercus wislizeni", "Quercus chrysolepis",         
    "Salix gooddingii", "Salix laevigata", "Salix exigua",        
    "Populus trichocarpa", "Populus fremontii", "Populus tremuloides",         
    "Larrea tridentata", "Cleomella arborea", "Acer macrophyllum",        
    "Aesculus californica", "Malosma laurina", "Vitis californica",           
    "Ribes sericeum", "Ribes cereum", "Platanus racemosa",           
    "Hesperoyucca whipplei", "Yucca schidigera", "Washingtonia filifera",       
    "Umbellularia californica", "Hesperocyparis macnabiana", "Hesperocyparis sargentii",    
    "Calocedrus decurrens", "Sequoiadendron giganteum", "Sequoia sempervirens",  
    "Pinus ponderosa", "Pinus sabiniana", "Pinus coulteri",       
    "Pinus flexilis", "Pinus lambertiana", "Pinus longaeva",              
    "Pseudotsuga menziesii", "Abies concolor", "Abies bracteata",             
    "Abies magnifica"
  ))

# Map the base tree
p <- ggtree(tree, layout = "circular", ladderize = FALSE) 

# T50
traits_t50_b <- t50 |> 
  group_by(species_3) |> 
  summarise(t50 = mean(T50)) |> 
  ungroup()

# Map the trait data to the tree
tree_t50 <- p %<+% traits_t50_b 

# Prepare Data for Custom Tip Labels
tip_t50<- tree_t50$data %>%
  filter(isTip == TRUE) %>%
  dplyr::select(label, x, y, node, angle, isTip) 

tiplab_t50 <- tip_t50 %>%
  left_join(name_map_df, by = c("label" = "original_label")) %>%
  mutate(display_label_final = ifelse(is.na(display_label), 
                                      label, display_label))
# Aesthetics
fig2 <- tree_t50 +
  geom_tree(aes(color = t50),
            linewidth = 1.5) +
  scale_color_viridis(option = "B", 
                      name =  expression("T"[50]~"(°C)")) +
  geom_tiplab(linewidth = 2.8, 
              fontface = "italic",
              offset = 15, 
              data = tiplab_t50, 
              mapping = aes(label = display_label_final),
              align = TRUE,
              linesize = 0.00, 
              linetype = "blank", 
              inherit.aes = FALSE) +  
  xlim(0, 420) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.8, "cm"), 
        legend.text = element_text(size = 12), 
        title = element_text(size = 14)) +
  geom_treescale(fontsize = 3,
                 linesize = 1, 
                 offset = 0.8)
fig2 

# Figure 3 ----
# TCrit
# Best model: Max temp from the last 15 days x Mean temp from the last 15 days
px15m15_2 <- plot(conditional_smooths(mod_temp2_x15m15, 
                                      smooths = "s(maxt_15day, meant_15day)"), 
                  ask = FALSE,
                  stype = "raster")[[1]]

px15m15_2$data$estimate__[px15m15_2$data$meant_15day > px15m15_2$data$maxt_15day] <- 
  NA_real_

px15m15_2 <- px15m15_2 +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent") +
  theme_classic(base_size = 13) +
  labs(x = "Maximum temperature 15D (°C)", 
       y = "Mean temperature 15D (°C)", 
       fill = expression(Delta*T[Crit]*" (°C)")) +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18)) 

# Second best model: Max temp from the last day x Mean temp from the last 15 days 
px1m15_2 <- plot(conditional_smooths(mod_temp2_x1m15, 
                                     smooths = "s(maxt_1day, meant_15day)"), 
                 ask = FALSE,
                 stype = "raster")[[1]] 

px1m15_2$data$estimate__[px1m15_2$data$meant_15day > px1m15_2$data$maxt_1day] <- 
  NA_real_

px1m15_2 <- px1m15_2 +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent") +
  theme_classic(base_size = 13) +
  labs(x = "Maximum temperature 1D (°C)", 
       y = "Mean temperature 15D (°C)",
       fill = expression(Delta*T[Crit]*" (°C)")) +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18))

# T50
# Best model: Max temp from the last 15 days x Mean temp from the last 15 days
px15m15 <- plot(conditional_smooths(mod_temp_x15m15, 
                                    smooths = "s(maxt_15day, meant_15day)"), 
                ask = FALSE,
                stype = "raster")[[1]] 

px15m15$data$estimate__[px15m15$data$meant_15day > px15m15$data$maxt_15day] <- 
  NA_real_

px15m15 <- px15m15 +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent") +
  theme_classic(base_size = 13) +
  labs(x = "Maximum temperature 15D (°C)", 
       y = "Mean temperature 15D (°C)",
       fill = expression(Delta*T[50]*" (°C)")) +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18))

# Second best model: Max temp from the last day x Mean temp from the last 15 days 
px1m15 <- plot(conditional_smooths(mod_temp_x1m15, 
                                   smooths = "s(maxt_1day, meant_15day)"), 
               ask = FALSE,
               stype = "raster")[[1]] 

px1m15$data$estimate__[px1m15$data$meant_15day > px1m15$data$maxt_1day] <- 
  NA_real_

px1m15 <- px1m15 +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent") +
  theme_classic(base_size = 13) +
  labs(x = "Maximum temperature 1D (°C)", 
       y = "Mean temperature 15D (°C)",
       fill = expression(Delta*T[50]*" (°C)")) +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18))

# Taking out axis titles, ticks and numbers to make final figure less convoluted
px15m15_2 <- px15m15_2 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

px1m15_2 <- px1m15_2 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

px1m15 <- px1m15 +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Create the title (TCrit) for the top row as a separate object
title_tcrit <- ggdraw() + 
  draw_label(expression("T"[Crit]~"(°C)"),
             fontface = 'bold',
             size = 17,
             x = 0.1,  
             hjust = 0.5)

# Combine plots (a) and (b) side-by-side (Tcrit graphs)
row1_plots <- plot_grid(px15m15_2, px1m15_2, 
                        ncol = 2, 
                        labels = c("(a)", "(b)"),
                        label_x = 0.89,
                        hjust = 1,
                        label_fontface = "plain") 

# Combine title and Tcrit graphs
final_row1 <- plot_grid(title_tcrit, 
                        row1_plots, 
                        ncol = 1, 
                        rel_heights = c(0.1, 0.9))

# Create the title (T50) for the bottom row
title_t50 <- ggdraw() + 
  draw_label(expression("T"[50]~"(°C)"),
             fontface = 'bold',
             size = 17,
             x = 0.1,
             hjust = 0.5)

# Combine plots (c) and (d) side-by-side (T50 graps)
row2_plots <- plot_grid(px15m15, px1m15, 
                        ncol = 2, 
                        labels = c("(c)", "(d)"),
                        label_x = 0.89,
                        hjust = 1,
                        label_fontface = "plain")

# Stack the T50 title on top of its plot row
final_row2 <- plot_grid(title_t50, 
                        row2_plots, 
                        ncol = 1, 
                        rel_heights = c(0.1, 0.9))

# Combine the two final rows vertically
fig3 <- plot_grid(final_row1, 
                  NULL, 
                  final_row2, 
                  ncol = 1, 
                  rel_heights = c(1, 0.13, 1))

fig3

# Figure 4 ----
# Function to extract the specific trait row from the models
get_brms_est <- function(model, trait_name, label, heattolerance) {
  summ <- fixef(model) 
  row_data <- as.data.frame(t(summ[trait_name, ])) 
  row_data$Trait <- label
  row_data$HT <- heattolerance
  return(row_data)
}

# Extract data from individual models
# T50
res_tau <- get_brms_est(mod_t50_tau, "tau", "tau", "T50")
res_lma <- get_brms_est(mod_t50_lma, "lma", "LMA", "T50")
res_ldmc <- get_brms_est(mod_t50_ldmc, "ldmc", "LDMC", "T50")
res_lwc <- get_brms_est(mod_t50_lwc, "lwc", "LWC", "T50")
res_la <- get_brms_est(mod_t50_la, "la", "Leaf area", "T50")
res_lt <- get_brms_est(mod_t50_lt, "lt", "Leaf thickness", "T50")
res_lw <- get_brms_est(mod_t50_lw, "lw", "Leaf width", "T50")
res_rwc <- get_brms_est(mod_t50_rwc, "rwc", "RWC", "T50")
res_ewt <- get_brms_est(mod_t50_ewt, "ewt", "EWT", "T50")
res_d13c <- get_brms_est(mod_t50_d13C, "dc", "delta 13 C", "T50")
res_tn <- get_brms_est(mod_t50_tn, "t_n", "Total N", "T50")

# TCrit
res_tau2 <- get_brms_est(mod_t15_tau, "tau", "tau", "T15")
res_lma2 <- get_brms_est(mod_t15_lma, "lma", "LMA", "T15")
res_ldmc2 <- get_brms_est(mod_t15_ldmc, "ldmc", "LDMC", "T15")
res_lwc2 <- get_brms_est(mod_t15_lwc, "lwc", "LWC", "T15")
res_la2 <- get_brms_est(mod_t15_la, "la", "Leaf area", "T15")
res_lt2 <- get_brms_est(mod_t15_lt, "lt", "Leaf thickness", "T15")
res_lw2 <- get_brms_est(mod_t15_lw, "lw", "Leaf width", "T15")
res_rwc2 <- get_brms_est(mod_t15_rwc, "rwc", "RWC", "T15")
res_ewt2 <- get_brms_est(mod_t15_ewt, "ewt", "EWT", "T15")
res_d13c2 <- get_brms_est(mod_t15_d13C, "dc", "delta 13 C", "T15")
res_tn2 <- get_brms_est(mod_t15_tn, "t_n", "Total N", "T15")

# Combine data from all models
plot_data <- bind_rows(res_lma, res_ldmc, res_tau,res_lwc, res_la, res_lt, 
                       res_lw, res_d13c,res_tn, 
                       res_lma2, res_ldmc2, res_tau2, res_lwc2, res_la2, res_lt2, 
                       res_lw2, res_rwc, res_ewt, res_rwc2, res_ewt2,
                       res_d13c2, res_tn2) 

colnames(plot_data) <- c("Estimate", "Error", "Lower_CI", 
                         "Upper_CI", "Trait", "Heat_tolerance")

# Plot
fig4 <- ggplot(plot_data, 
               aes(x = Estimate, 
                   y = reorder(Trait, Estimate),
                   color =  Heat_tolerance)) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "gray50") +
  geom_errorbarh(aes(xmin = Lower_CI, 
                     xmax = Upper_CI), 
                 height = 0.2, 
                 linewidth = 1,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  labs(x = "Standardized Effect Size",
       y = NULL,
       color = NULL) +
  theme_classic(base_size = 26) +
  theme(panel.grid.major.y = element_line(color = "gray95",
                                          linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray95",
                                          linewidth = 0.5),
        legend.position = "inside",
        legend.position.inside = c(.15, 1),
        legend.background = element_blank(), 
        legend.direction = "horizontal") +
  scale_color_manual(values = c("T50" =  "#201174", "T15" = "#BCBDDC"),
                     breaks = c("T50", "T15"),
                     labels = c(expression(T[50]), expression(T[Crit]))) +
  scale_y_discrete(labels = c(
    "delta 13 C"   = expression(delta^13*C),
    "tau" = expression(tau),
    "Total N" = "Leaf N"))

fig4

# Figure 5 ----

# Select data for PCA
data_pca <- t50 |> 
  dplyr::select(T50, T15, TAU, LMA, LDMC, leaf_area, leaf_width) %>%
  mutate(LMA = log(LMA),
         TAU = log(TAU),
         leaf_area = log(leaf_area),
         leaf_width = log(leaf_width))

# PCA
pca_1 <- princomp(data_pca, cor = TRUE)

# Get the variance explained by each component
summary(pca_1)

# Get the loadings 
pca_loadings <- pca_1$loadings
print(pca_loadings)
print(unclass(pca_loadings))

# Get the original variable names from your PCA object
original_var_names <- names(pca_1$center)

# Create a data frame that maps original names to your new, formatted labels
name_map <- data.frame(
  original_name = original_var_names,
  new_label = c("T[50]", "T[Crit]", "tau", "LMA", "LDMC", "LA", "LW"))

#  Extract variable coordinates and join with new labels
var_coords <- as.data.frame(get_pca_var(pca_1)$coord)
var_coords$original_name <- rownames(var_coords)

# Join with your name map
var_coords <- var_coords |> 
  left_join(name_map, by = "original_name")

# Graph
fig5_base <- fviz_pca_biplot(pca_1,
                             col.var = "#ef3b2c",      
                             geom.ind = "point",
                             geom.var = "arrow",      
                             pointshape = 19,
                             arrowsize = 1,
                             pointsize = 2.5, 
                             col.ind = interaction(t50$leaf_habit_spec, 
                                                   t50$seed_plant),
                             ggtheme = theme_classic(base_size = 20),
                             title = NULL,
                             legend.title = " ", 
                             addEllipses = TRUE, 
                             label = "var")

fig5 <- fig5_base +
  geom_text_repel(data = var_coords,
                  aes(x = Dim.1 * 4.5, y = Dim.2 * 5, label = new_label), 
                  parse = TRUE,  
                  size = 5) +
  theme(legend.position = "top", 
        legend.title = element_text(""), 
        legend.text = element_text(size = 15)) +
  labs(x = "PC1 (38.1%)", 
       y = "PC2 (29.9%)") + 
  scale_color_manual(name = NULL,  
                     values = c("drought_deciduous.angiosperm" = "#4d4352",
                                "winter_deciduous.angiosperm" = "#b694c9",
                                "evergreen.angiosperm" = "#bdd7e7",
                                "evergreen.gymnosperm" = "#6e9948"),
                     labels = c("drought_deciduous.angiosperm" = "Drought-Deciduous Angiosperm",
                                "winter_deciduous.angiosperm" = "Winter-Deciduous Angiosperm",
                                "evergreen.angiosperm" = "Evergreen Angiosperm",
                                "evergreen.gymnosperm" = "Evergreen Gymnosperm"),
                     guide = guide_legend(nrow = 2)) +
  scale_fill_manual(values = c("drought_deciduous.angiosperm" = "#4d4352",
                               "winter_deciduous.angiosperm" = "#b694c9",
                               "evergreen.angiosperm" = "#bdd7e7",
                               "evergreen.gymnosperm" = "#6e9948"),
                    guide = "none")

fig5

# Supplementary Materials ------------------------------------------------------

# Figure S1 ----
# TCrit
# Bake titles directly into the grid objects
temp2_x15m15_MCMC <- arrangeGrob(mod_temp2_x15m15_MCMC, top = "Model Diagnostics (MCMC): 15-Day Mean and Max Temperature")
temp2_x1m15_MCMC <- arrangeGrob(mod_temp2_x1m15_MCMC, top = "Model Diagnostics (MCMC): 1-Day Mean and 15-Day Max Temperature")

figS1 <- wrap_plots(temp2_x15m15_MCMC, temp2_x1m15_MCMC, ncol = 1)
figS1

# Figure S2 ----
# T50
# Bake titles directly into the grid objects
x15m15_MCMC <- arrangeGrob(mod_temp_x15m15_MCMC, top = "Model Diagnostics (MCMC): 15-Day Mean and Max Temperature")
x1m15_MCMC <- arrangeGrob(mod_temp_x1m15_MCMC, top = "Model Diagnostics (MCMC): 1-Day Mean and 15-Day Max Temperature")

figS2 <- wrap_plots(x15m15_MCMC, x1m15_MCMC, ncol = 1)
figS2

# Figure S3 ----

# Change names of the plots
new_plot_names <- c("sd_species_3__Intercept" = "Species SD", 
                    "sds_scoordinates_longitudecoordinates_latitude_1" = "Spatial Smooth SD", 
                    "Intercept" = "Intercept")

# T50
trace_temp30_t50 <- mcmc_trace(mod_temp30_t50,
                               regex_pars = "(^sd|^Intercept)",
                               facet_args = list(
                                 labeller = as_labeller(new_plot_names)))

do_temp30_t50 <- mcmc_dens_overlay(mod_temp30_t50,
                                   regex_pars = "(^sd|^Intercept)", 
                                   facet_args = list(
                                     labeller = as_labeller(new_plot_names)))

# TCrit
trace_temp30_t15 <- mcmc_trace(mod_temp30_t15,
                               regex_pars = "(^sd|^Intercept)",
                               facet_args = list(
                                 labeller = as_labeller(new_plot_names)))

do_temp30_t15 <- mcmc_dens_overlay(mod_temp30_t15,
                                   regex_pars = "(^sd|^Intercept)",
                                   facet_args = list(
                                     labeller = as_labeller(new_plot_names)))

# Combine figures that asses convergence of MCMC
figS3 <- do_temp30_t15 +
  labs(tag = "(a)") +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        legend.position = "none", 
        text = element_text(size = 16)) +
  trace_temp30_t15 +
  theme(legend.position = "none", 
        text = element_text(size = 16)) +
  do_temp30_t50 +
  labs(tag = "(b)") +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none", 
        text = element_text(size = 16)) + 
  trace_temp30_t50 + 
  theme(legend.position = "none",
        text = element_text(size = 16)) +
  plot_layout(ncol = 1)

figS3

# Figure S4 ----
# Find the exact parameter names from the model
variables(mod_ceeb_t50)

# Change names of the plots
new_plot_names <- c("sd_species_3__Intercept" = "Species SD", 
                    "sds_scoordinates_longitudecoordinates_latitude_1" = "Spatial Smooth SD", 
                    "Intercept" = "Intercept")

# T50
trace_ceeb_t50 <- mcmc_trace(mod_ceeb_t50,
                             regex_pars = "(^sd|^Intercept)", 
                             facet_args = list(
                               labeller = as_labeller(new_plot_names)))

do_ceeb_t50 <- mcmc_dens_overlay(mod_ceeb_t50,
                                 regex_pars = "(^sd|^Intercept)", 
                                 facet_args = list(
                                   labeller = as_labeller(new_plot_names)))

# TCrit
trace_ceeb_t15 <- mcmc_trace(mod_ceeb_t15,
                             regex_pars = "(^sd|^Intercept)",
                             facet_args = list(
                               labeller = as_labeller(new_plot_names)))

do_ceeb_t15 <- mcmc_dens_overlay(mod_ceeb_t15,
                                 regex_pars = "(^sd|^Intercept)",
                                 facet_args = list(
                                   labeller = as_labeller(new_plot_names)))

# Combine figures that asses convergence of MCMC
figS4 <- do_ceeb_t15 +
  labs(tag = "(a)") +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        legend.position = "none", 
        text = element_text(size = 16)) +
  trace_ceeb_t15 +
  theme(legend.position = "none", 
        text = element_text(size = 16)) +
  do_ceeb_t50 +
  labs(tag = "(b)") +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 18), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none", 
        text = element_text(size = 16)) + 
  trace_ceeb_t50 + 
  theme(legend.position = "none",
        text = element_text(size = 16)) +
  plot_layout(ncol = 1)

figS4

# Figure S5 ----
traits_t15_b <- t50 |> 
  group_by(species_3) |> 
  summarise(t15 = mean(T15)) |> 
  ungroup() 

# Map the trait data to the tree
tree_t15 <- p %<+% traits_t15_b 

# Prepare Data for Custom Tip Labels
tip_t15 <- tree_t15$data %>%
  filter(isTip == TRUE) %>%
  dplyr::select(label, x, y, node, angle, isTip) 

tiplab_t15 <- tip_t15 %>%
  left_join(name_map_df, by = c("label" = "original_label")) %>%
  mutate(display_label_final = ifelse(is.na(display_label), 
                                      label, display_label))

# Aesthetics
figS5 <- tree_t15 +
  geom_tree(aes(color = t15),
            size = 1.5) +
  scale_color_viridis(option = "B", 
                      name = expression("T"[15]~"(°C)")) +
  geom_tiplab(size = 4, 
              fontface = "italic",
              offset = 15,
              data = tiplab_t15, 
              mapping = aes(label = display_label_final),
              align = TRUE,
              linesize = 0.00, 
              linetype = "blank",
              inherit.aes = FALSE) +  
  xlim(0, 420) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 18), 
        title = element_text(size = 22)) +
  geom_treescale(fontsize = 3,
                 linesize = 1, 
                 offset = 0.8)

figS5

# Figure S6 ----
traits_tau_b <- t50 |> 
  group_by(species_3) |> 
  summarise(tau = mean(tau)) |> 
  ungroup()

# Map the trait data to the tree
tree_tautrait <- p %<+% traits_tau_b

# Prepare Data for Custom Tip Labels
tip_tau <- tree_tautrait$data %>%
  filter(isTip == TRUE) %>%
  dplyr::select(label, x, y, node, angle, isTip) 

tiplab_tau <- tip_tau %>%
  left_join(name_map_df, by = c("label" = "original_label")) %>%
  mutate(display_label_final = ifelse(is.na(display_label), 
                                      label, display_label))
# Aesthetics
figS6 <- tree_tautrait +
  geom_tree(aes(color = log(tau)),
            size = 1.5) +
  scale_color_viridis(option = "B", 
                      name = "log(tau)") +
  geom_tiplab(size = 4,
              fontface = "italic",
              offset = 15,
              data = tiplab_tau, 
              mapping = aes(label = display_label_final),
              align = TRUE,
              linesize = 0.00, 
              linetype = "blank",
              inherit.aes = FALSE) + 
  xlim(0, 420) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.8, "cm"), 
        legend.text = element_text(size = 18), 
        title = element_text(size = 22)) +
  geom_treescale(fontsize = 3,
                 linesize = 1, 
                 offset = 0.8)

figS6

# Figure S7 ----
traits_lma_b <- t50 |> 
  group_by(species_3) |> 
  summarise(lma = mean(LMA)) |> 
  ungroup() 

# Map the trait data to the tree
tree_lma <- p %<+% traits_lma_b 

# Prepare Data for Custom Tip Labels
tip_lma <- tree_lma$data %>%
  filter(isTip == TRUE) %>%
  dplyr::select(label, x, y, node, angle, isTip) 

tiplab_lma <- tip_lma %>%
  left_join(name_map_df, by = c("label" = "original_label")) %>%
  mutate(display_label_final = ifelse(is.na(display_label), 
                                      label, display_label))
# Aesthetics
figS7 <-tree_lma +
  geom_tree(aes(color = log(lma)),
            size = 1.5) +
  scale_color_viridis(option = "B",
                      name = "log(LMA)") +
  geom_tiplab(size = 4,
              fontface = "italic",
              offset = 15, 
              data = tiplab_lma, 
              mapping = aes(label = display_label_final),
              align = TRUE,
              linesize = 0.00, 
              linetype = "blank",
              inherit.aes = FALSE) + 
  xlim(0, 420) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.8, "cm"), 
        legend.text = element_text(size = 18), 
        title = element_text(size = 22)) +
  geom_treescale(fontsize = 3,
                 linesize = 1, 
                 offset = 0.8)

figS7 

# Figure S8 ----
traits_ldmc_b <- t50 |> 
  group_by(species_3) |> 
  summarise(ldmc = mean(LDMC)) |> 
  ungroup() 

# Map the trait data to the tree
tree_ldmc <- p %<+% traits_ldmc_b 

# Prepare Data for Custom Tip Labels
tip_ldmc <- tree_ldmc$data %>%
  filter(isTip == TRUE) %>%
  dplyr::select(label, x, y, node, angle, isTip) 

tiplab_ldmc <- tip_ldmc %>%
  left_join(name_map_df, by = c("label" = "original_label")) %>%
  mutate(display_label_final = ifelse(is.na(display_label), 
                                      label, display_label))
# Aesthetics
figS8 <-tree_ldmc +
  geom_tree(aes(color = ldmc),
            size = 1.5) +
  scale_color_viridis(option = "B",
                      name = "LDMC") +
  geom_tiplab(size = 4,
              fontface = "italic",
              offset = 15,
              data = tiplab_ldmc, 
              mapping = aes(label = display_label_final),
              align = TRUE,
              linesize = 0.00, 
              linetype = "blank",
              inherit.aes = FALSE) + 
  xlim(0, 420) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.8, "cm"), 
        legend.text = element_text(size = 18), 
        title = element_text(size = 22)) +
  geom_treescale(fontsize = 3,
                 linesize = 1, 
                 offset = 0.8)

figS8

# Figure S9 ----
traits_lwc_b <- t50 |> 
  group_by(species_3) |> 
  summarise(lwc = mean(LWC)) |> 
  ungroup() 

# Map the trait data to the tree
tree_lwc <- p %<+% traits_lwc_b 

# Prepare Data for Custom Tip Labels
tip_lwc <- tree_lwc$data %>%
  filter(isTip == TRUE) %>%
  dplyr::select(label, x, y, node, angle, isTip) 

tiplab_lwc <- tip_lwc %>%
  left_join(name_map_df, by = c("label" = "original_label")) %>%
  mutate(display_label_final = ifelse(is.na(display_label), 
                                      label, display_label))
# Aesthetics
figS9 <- tree_lwc +
  geom_tree(aes(color = lwc),
            size = 1.5) +
  scale_color_viridis(option = "B", 
                      name = "LWC") +
  geom_tiplab(size = 4,
              fontface = "italic",
              offset = 15,
              data = tiplab_lwc, 
              mapping = aes(label = display_label_final),
              align = TRUE,
              linesize = 0.00, 
              linetype = "blank",
              inherit.aes = FALSE) + 
  xlim(0, 420) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.8, "cm"), 
        legend.text = element_text(size = 18), 
        title = element_text(size = 22)) +
  geom_treescale(fontsize = 3,
                 linesize = 1, 
                 offset = 0.8)

figS9

# Figure S10 ----
newdata_elevation <- t50 %>%
  data_grid(
    elevation = seq_range(elevation, n = 100), 
    maxt_15day = median(maxt_15day, na.rm = TRUE),
    meant_15day = median(meant_15day, na.rm = TRUE),
    coordinates_longitude = mean(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = mean(coordinates_latitude, na.rm = TRUE))

# T50
preds_elevation <- newdata_elevation %>%
  add_epred_draws(mod_temp_x15m15, 
                  re_formula = NA) 

# Back-transform variable
mean_elevation <- mean(t50$Elevation, na.rm = TRUE)

preds_elevation <- preds_elevation %>%
  mutate(elevation = (elevation * 500) + mean_elevation)

# Plot the estimates
pred_t50_ele <- ggplot(preds_elevation, aes(x = elevation, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#08519C") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = Elevation, y = T50),
             alpha = 0.2, color = "grey40") +
  labs(x = "Elevation (m asl)", 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 25) +
  labs(tag = "(b)") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 25))

# TCrit
# Generate posterior fitted draws for the new data
preds_elevation2 <- newdata_elevation %>%
  add_epred_draws(mod_temp2_x15m15, 
                  re_formula = NA) 

# Back-transform variable
preds_elevation2 <- preds_elevation2 %>%
  mutate(elevation = (elevation * 500) + mean_elevation)

# Plot the estimates
pred_t15_ele <- ggplot(preds_elevation2, aes(x = elevation, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#08519C") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = Elevation, y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = NULL, 
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 25) +
  labs(tag = "(a)") +
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 25))

# Plot
figS10 <- wrap_plots(pred_t15_ele, pred_t50_ele, ncol = 1)
figS10

# Figure S11 ----
pcoord <- plot(conditional_smooths(
  mod_temp_x15m15, 
  smooths = "s(coordinates_longitude, coordinates_latitude)"), 
  ask = FALSE,
  stype = "raster")[[1]] +
  theme_classic(base_size = 20) +
  labs(x = "Longitude", 
       y = "Latitude") +
  theme(legend.title = element_blank())

pcoord_2 <- plot(conditional_smooths(
  mod_temp2_x15m15, 
  smooths = "s(coordinates_longitude, coordinates_latitude)"), 
  ask = FALSE,
  stype = "raster")[[1]] +
  theme_classic(base_size = 20) +
  labs(x = "Longitude", 
       y = "Latitude") +
  theme(legend.title = element_blank())

pcoord <- pcoord +
  labs(tag = "(b)") +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 20))

pcoord_2 <- pcoord_2 +
  labs(tag = "(a)") +
  theme(plot.tag.position = "topright",
        plot.tag = element_text(size = 20))

figS11 <- wrap_plots(pcoord_2, pcoord, ncol = 2)
figS11

# Figure S12 ----
# tau 
# T50
newdata_tau <- t50 |> 
  data_grid(
    tau = seq_range(tau, n = 100), 
    coordinates_longitude = median(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = median(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1])

# Generate posterior fitted draws for the new data
predictions_tau <- newdata_tau |> 
  add_epred_draws(mod_t50_tau, 
                  re_formula = NULL) 

# Back-transform variable
mean_tau <- mean(log(t50$TAU))
sd_tau <- sd(log(t50$TAU))

predictions_tau <- 
  predictions_tau |> 
  ungroup() |> 
  mutate(tau_unscaled = (tau * sd_tau) + mean_tau,
         .before = tau)

# Plot the estimates
pred_t50_tau <- ggplot(predictions_tau, aes(x = tau_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(TAU), y = T50),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("Log-transformed" ~ tau ~ "(s)"), 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(m)") +
  theme(legend.position = "none",
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# TCrit
predictions_tau2 <- newdata_tau |> 
  add_epred_draws(mod_t15_tau, 
                  re_formula = NULL) 

# Back-transform variable
predictionstau2 <- 
  predictions_tau2 |> 
  ungroup() |> 
  mutate(tau_unscaled = (tau * sd_tau) + mean_tau,
         .before = tau)

# Plot the estimates
pred_t15_tau <- ggplot(predictions_tau2, aes(x = tau_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(TAU), y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("Log-transformed" ~ tau ~ "(s)"),
       y = c,
       fill = "Credible Interval") +
  theme_classic(base_size = 15) +
  labs(tag = "(a)") +
  theme(legend.position = "top",
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# LMA 
# Create a newdata grid to predict heat tolerance  across LMA
newdata_lma <- t50 |>
  data_grid(
    lma = seq_range(lma, n = 100),
    coordinates_longitude = median(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = median(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1])

# T50
summary(mod_t50_lma)
predictions_lma <- newdata_lma |>
  add_epred_draws(mod_t50_lma, 
                  re_formula = NULL) 

# Back-transform variable
mean_lma <- mean(log(t50$LMA))
sd_lma <- sd(log(t50$LMA))

predictions_lma <- predictions_lma |> 
  ungroup() |> 
  mutate(lma_unscaled = (lma * sd_lma) + mean_lma,
         .before = lma)

# Plot the estimates
pred_t50_lma <- ggplot(predictions_lma, aes(x = lma_unscaled, 
                                            y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(LMA), y = T50),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("Log-transformed LMA (kg m"^{-2}~")"), 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(g)") +
  theme(legend.position = "none", 
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# TCrit
predictions_lma2 <- newdata_lma |>
  add_epred_draws(mod_t15_lma, 
                  re_formula = NULL) 

# Unscaled varibles
predictions_lma2 <- 
  predictions_lma2 |>
  ungroup() |>
  mutate(lma_unscaled = lma * sd_lma + mean_lma,
         .before = lma)

# Plot the estimations
pred_t15_lma <- ggplot(predictions_lma2, aes(x = lma_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(LMA), y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("Log-transformed LMA (kg m"^{-2}~")"), 
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(a)") +
  theme(legend.position = "none",
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# LDMC 
# T50
newdata_ldmc <- t50 |>
  data_grid(
    ldmc = seq_range(ldmc, n = 100), 
    coordinates_longitude = median(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = median(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1]
  )

# Generate posterior fitted draws for the new data
predictions_ldmc <- newdata_ldmc |>
  add_epred_draws(mod_t50_ldmc, 
                  re_formula = NULL)

# Back-transform variable
mean_ldmc <- mean(t50$LDMC)
sd_ldmc <- sd(t50$LDMC)

predictions_ldmc <- 
  predictions_ldmc |>
  ungroup() |>
  mutate(ldmc_unscaled = (ldmc * sd_ldmc) + mean_ldmc,
         .before = ldmc)

# Plot the estimates
pred_t50_ldmc <- ggplot(predictions_ldmc, aes(x = ldmc_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = LDMC, y = T50),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("LDMC (kg kg"^{-1}~")"), 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(h)") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# TCrit
predictions_ldmc2 <- newdata_ldmc |>
  add_epred_draws(mod_t15_ldmc, 
                  re_formula = NULL) 

predictions_ldmc2 <- 
  predictions_ldmc2 |>
  ungroup() |>
  mutate(ldmc_unscaled = (ldmc * sd_ldmc) + mean_ldmc,
         .before = ldmc)

# Plot the estimates
pred_t15_ldmc <- ggplot(predictions_ldmc2, aes(x = ldmc_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = LDMC, y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("LDMC (kg kg"^{-1}~")"),
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(b)") +
  theme(legend.position = "none",
        plot.tag.position = "topright",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# LWC 
newdata_lwc <- t50 |>
  data_grid(
    lwc = seq_range(lwc, n = 100), 
    coordinates_longitude = median(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = median(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1]
  )

# T50
predictions_lwc <- newdata_lwc |>
  add_epred_draws(mod_t50_lwc, 
                  re_formula = NULL) 

# Back-transform variable
mean_lwc <- mean(t50$LWC)
sd_lwc <- sd(t50$LWC)

predictions_lwc <- 
  predictions_lwc |>
  ungroup() |>
  mutate(lwc_unscaled = (lwc * sd_lwc) + mean_lwc,
         .before = lwc)

# Plot the predictions for elevation
pred_t50_lwc <- ggplot(predictions_lwc, aes(x = lwc_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = LWC, y = T50), 
             alpha = 0.2, color = "grey40") +
  labs(x = "LWC (%)", 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(i)") +
  theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# TCrit
predictions_lwc2 <- newdata_lwc |>
  add_epred_draws(mod_t15_lwc, 
                  re_formula = NULL) 

# Back-transform variable
predictions_lwc2 <- 
  predictions_lwc2 |>
  ungroup() |>
  mutate(lwc_unscaled = (lwc * sd_lwc) + mean_lwc,
         .before = lwc)

# Plot the predictions for elevation
pred_t15_lwc <- ggplot(predictions_lwc2, aes(x = lwc_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = LWC, y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = "LWC (%)", 
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(c)") +
  theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# LA 
newdata_la <- t50 |>
  data_grid(
    la = seq_range(la, n = 100), 
    coordinates_longitude = mean(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = mean(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1]
  )

# T50
predictions_la <- newdata_la |>
  add_epred_draws(mod_t50_la, 
                  re_formula = NULL) 

# Unscaled varibles
mean_la <- mean(log(t50$leaf_area))
sd_la <- sd(log(t50$leaf_area))

predictions_la <- predictions_la |> 
  ungroup() |> 
  mutate(la_unscaled = la * sd_la + mean_la,
         .before = la)

# Plot the estimates
pred_t50_la <- ggplot(predictions_la, aes(x = la_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(leaf_area), y = T50),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("Log-transformed leaf area (m"^{2}~")"), 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(j)") +
  theme(legend.position = "none",
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# TCrit
predictions_la2 <- newdata_la |>
  add_epred_draws(mod_t15_la, 
                  re_formula = NULL) 

# Unscaled varibles
predictions_la2 <- predictions_la2 |> 
  ungroup() |> 
  mutate(la_unscaled = la * sd_la + mean_la,
         .before = la)

# Plot the estimates
pred_t15_la <- ggplot(predictions_la2, aes(x = la_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(leaf_area), y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = expression("Log-transformed leaf area (m"^{2}~")"), 
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(d)") +
  theme(legend.position = "none",
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# LW 
newdata_lw <- t50 |>
  data_grid(
    lw = seq_range(lw, n = 100), 
    coordinates_longitude = mean(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = mean(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1]
  )

# T50
predictions_lw <- newdata_lw |>
  add_epred_draws(mod_t50_lw, 
                  re_formula = NULL) 

# Unscaled varibles
mean_lw <- mean(log(t50$leaf_width))
sd_lw <- sd(log(t50$leaf_width))

predictions_lw <- predictions_lw |> 
  ungroup() |> 
  mutate(lw_unscaled = lw * sd_lw + mean_lw,
         .before = lw)

# Plot the estimates
pred_t50_lw <- ggplot(predictions_lw, aes(x = lw_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(leaf_width), y = T50),
             alpha = 0.2, color = "grey40") +
  labs(x = "Log-transformed leaf width (m)", 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(k)") +
  theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# TCrit
predictions_lw2 <- newdata_lw |>
  add_epred_draws(mod_t15_lw, 
                  re_formula = NULL) 

# Unscaled varibles
predictions_lw2 <- predictions_lw2 |> 
  ungroup() |> 
  mutate(lw_unscaled = lw * sd_lw + mean_lw,
         .before = lw)

# Plot the estimates
pred_t15_lw <- ggplot(predictions_lw2, aes(x = lw_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(leaf_width), y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = "Log-transformed leaf width (m)", 
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(e)") +
  theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# LT
newdata_lt <- t50 |>
  data_grid(
    lt = seq_range(lt, n = 100), 
    coordinates_longitude = mean(coordinates_longitude, na.rm = TRUE),
    coordinates_latitude = mean(coordinates_latitude, na.rm = TRUE),
    species_3 = t50$species_3[1]
  )

# T50
predictions_lt <- newdata_lt |>
  add_epred_draws(mod_t50_lt, 
                  re_formula = NULL)

# Unscaled varibles
mean_lt <- mean(log(t50$leaf_thickness))
sd_lt <- sd(log(t50$leaf_thickness))

predictions_lt <- predictions_lt |> 
  ungroup() |> 
  mutate(lt_unscaled = lt * sd_lt + mean_lt,
         .before = lt)

# Plot the estimates
pred_t50_lt <- ggplot(predictions_lt, aes(x = lt_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50, 
             aes(x = log(leaf_thickness), y = T50, color = leaf_habit),
             alpha = 0.2, 
             color = "grey40") +
  labs(x = "Log-transformed leaf thickness (m)", 
       y = expression("Estimated T"[50]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(l)") +
  theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3)) 

# TCrit
predictions_lt2 <- newdata_lt |>
  add_epred_draws(mod_t15_lt, 
                  re_formula = NULL) 

# Unscaled varibles
predictions_lt2 <- predictions_lt2 |> 
  ungroup() |> 
  mutate(lt_unscaled = lt * sd_lt + mean_lt,
         .before = lt)

# Plot the estimates
pred_t15_lt <- ggplot(predictions_lt2, aes(x = lt_unscaled, y = .epred)) +
  stat_lineribbon(.width = c(0.95, 0.5), 
                  color = "#201174") +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = t50,
             aes(x = log(leaf_thickness), y = T15),
             alpha = 0.2, color = "grey40") +
  labs(x = "Log-transformed leaf thickness (m)", 
       y = expression("Estimated T"[Crit]~"(°C)"),
       fill = "Credible Interval") +
  theme_classic(base_size = 13) +
  labs(tag = "(f)") +
  theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = "topright",
        plot.tag = element_text(size = 15),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.3))

# Plot 
figS12 <- wrap_plots(pred_t15_lma, pred_t15_ldmc, pred_t15_lwc, 
                     pred_t15_la, pred_t15_lw, pred_t15_lt, 
                     pred_t50_lma, pred_t50_ldmc, pred_t50_lwc, 
                     pred_t50_la, pred_t50_lw, pred_t50_lt,
                     ncol = 3)
figS12

# Figure S13 ----
t15_lma_lh <- ggplot() +
  geom_jitter(data = t50,
              aes(x = log(LMA),
                  y = T15,
                  colour = leaf_habit_spec),
              alpha = 0.6,
              size = 4) +
  theme_classic(base_size = 25) +
  annotate("text", 
           x = Inf, y = Inf, 
           label = "(a)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(x = expression("log(LMA) (kg m"^{-2}~")"), 
       y = expression("T"[Crit]~"(°C)")) +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  scale_colour_manual(values = c("#201174", "#bc96ff", "#ffb54b"),
                      breaks = c("drought_deciduous", "winter_deciduous", "evergreen"),
                      labels = c("Drought-deciduous", "Winter-deciduous","Evergreen")) +
  scale_y_continuous(limits = c(35, 60), breaks = seq(35, 60, by = 5)) 

t50_lma_lh <- ggplot() +
  geom_jitter(data = data_cali,
              aes(x = log(LMA),
                  y = T50,
                  colour = leaf_habit_spec),
              alpha = .6,
              size = 4) +
  theme_classic(base_size = 25) +
  annotate("text", 
           x = Inf, y = Inf, 
           label = "(b)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(x = expression("log(LMA) (kg m"^{-2}~")"), 
       y = expression("T"[50]~"(°C)")) +
  theme(legend.title = element_blank(),
        legend.position = "top", 
        legend.text = element_text(size = 15)) +
  scale_colour_manual(values = c("#201174", "#bc96ff", "#ffb54b"),
                      breaks = c("drought_deciduous", "winter_deciduous", "evergreen"),
                      labels = c("Drought-deciduous", "Winter-deciduous","Evergreen"),
                      guide = guide_legend(nrow = 2)) +
  scale_y_continuous(limits = c(35, 60), breaks = seq(35, 60, by = 5)) 

# Plot
figS13 <- wrap_plots(t15_lma_lh, t50_lma_lh, ncol = 2)
figS13

# Figure S14 ----

# Filtering out species observed at a single site
aux_species <- t50 |> 
  count(site, species_3) |> 
  group_by(species_3) |> 
  mutate(mult_site = n() > 1) |> 
  ungroup()

aux_species <- aux_species |> 
  filter(!mult_site) |> 
  pull(species_3)

# Pull repeated species in multiple sites
t50_rep_spp <- t50 |> 
  filter(!species_3 %in% aux_species) |> 
  dplyr::select(!c(genus, species_variation, date))

# Change species names 
spp_rep_names <- c("Acer_macrophyllum" = "Acer macrophyllum", 
                   "Arctostaphylos_glandulosa" = "Arctostaphylos glandulosa",
                   "Pseudotsuga_menziesii" = "Pseudotsuga menziesii",     
                   "Quercus_chrysolepis" = "Quercus chrysolepis",      
                   "Umbellularia_californica" = "Umbellularia californica",  
                   "Encelia_farinosa" = "Encelia farinosa",         
                   "Fouquieria_splendens" = "Fouquieria splendens",      
                   "Simmondsia_chinensis" = "Simmondsia chinensis",     
                   "Artemisia_tridentata" = "Artemisia tridentata",      
                   "Ribes_cereum" = "Ribes cereum",             
                   "Acmispon_glaber" = "Acmispon glaber",           
                   "Artemisia_californica" = "Artemisia californica",    
                   "Heteromeles_arbutifolia" = "Heteromeles arbutifolia",   
                   "Quercus_douglasii" = "Quercus douglasii",        
                   "Quercus_lobata" = "Quercus lobata",            
                   "Calocedrus_decurrens" = "Calocedrus decurrens",     
                   "Frangula_californica" = "Frangula californica",      
                   "Pinus_ponderosa" = "Pinus ponderosa",          
                   "Quercus_kelloggii" = "Quercus kelloggii",         
                   "Quercus_wislizeni" = "Quercus wislizeni",        
                   "Alnus_incana" = "Alnus incana",              
                   "Cornus_nuttallii" = "Cornus nuttallii",         
                   "Pinus_lambertiana" = "Pinus lambertiana",         
                   "Rubus_parviflorus" = "Rubus parviflorus",        
                   "Cercocarpus_betuloides" = "Cercocarpus betuloides",    
                   "Quercus_berberidifolia" = "Quercus berberidifolia",   
                   "Rosa_californica" = "Rosa californica",          
                   "Ericameria_nauseosa" = "Ericameria nauseosa",      
                   "Populus_fremontii" = "Populus fremontii",         
                   "Salix_exigua" = "Salix exigua",             
                   "Eriogonum_fasciculatum" = "Eriogonum fasciculatum",    
                   "Rhamnus_ilicifolia" = "Rhamnus ilicifolia")

# T50 of different populations based on the best temperature model 
# Model: mod_temp_x15m15 

# Get estimates graphs to have uncertainty quantification
dt_aux2 <- t50_rep_spp |>
  dplyr::select(site, species_3,
                coordinates_longitude, coordinates_latitude,
                elevation, maxt_15day, meant_15day) |>
  group_by(site, species_3) |>
  summarise(across(where(is.numeric), median)) |>
  ungroup() |>
  unique()

ppe2 <- posterior_epred(mod_temp_x15m15,
                        newdata = dt_aux2) |>
  posterior::summarise_draws("median", 
                             ~ posterior::quantile2(., probs = c(0.025, 0.975))) 

dt_aux2 <- dt_aux2 |>
  bind_cols(dplyr::select(ppe2, - variable))


# Identify species with non-overlapping intervals
significant_species <- dt_aux2 |>
  group_by(species_3) |>
  summarise(
    lowest_upper_bound = min(q97.5),
    highest_lower_bound = max(q2.5),
    is_distinct = highest_lower_bound > lowest_upper_bound) |>
  filter(is_distinct) |>
  pull(species_3)

# View the list of species
print(significant_species)

# Create a darker background for species that have different T50 among sites
# based on the previous model:

# Ensure significant_species is a character vector
significant_species <- as.character(significant_species)

# Create the highlight dataframe
highlight_df <- data.frame(species_3 = significant_species, 
                           start = -Inf, 
                           end = Inf, 
                           bottom = -Inf, 
                           top = Inf)

# Plot
figS14 <- ggplot() +
  geom_rect(data = highlight_df, 
            aes(xmin = start, xmax = end, 
                ymin = bottom, ymax = top), 
            fill = "grey80", 
            alpha = 0.5,
            inherit.aes = FALSE) +
  geom_linerange(aes(ymin = q2.5, ymax = q97.5, x = site), 
                 data = dt_aux2) +
  geom_point(aes(x = site, y = median), 
             data = dt_aux2) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.25, 
                                        linetype = 1)) +
  labs(x = expression("Sites"), 
       y = expression("Estimated T"[50]~"(°C)")) +
  facet_wrap(~ species_3, 
             scales = "free",
             ncol = 4,
             labeller = as_labeller(spp_rep_names)) +
  facetted_pos_scales( 
    y = list(
      !COL == 1 ~ scale_y_continuous(labels = NULL, 
                                     limits = c(45, 58), 
                                     breaks = seq(45, 55, by = 5)),
      COL == 1 ~ scale_y_continuous(limits = c(45, 58), 
                                    breaks = seq(45, 55, by = 5))))

figS14

# Figure S15 ----

# TCrit of different populations based on the best temperature model (max15D, mean15D) 
# Model: mod_temp2_x15m15 

# Get estimates graphs to have uncertainty quantification
dt_aux3 <- t50_rep_spp |>
  dplyr::select(site, species_3,
                coordinates_longitude, coordinates_latitude,
                elevation, maxt_15day, meant_15day) |>
  group_by(site, species_3) |>
  summarise(across(where(is.numeric), median)) |>
  ungroup() |>
  unique()

ppe3 <- posterior_epred(mod_temp2_x15m15,
                        newdata = dt_aux3) |>
  posterior::summarise_draws("median", 
                             ~ posterior::quantile2(., probs = c(0.025, 0.975))) 

dt_aux3 <- dt_aux3 |>
  bind_cols(dplyr::select(ppe3, - variable))

# Identify species with non-overlapping intervals
significant_species2 <- dt_aux3 |>
  group_by(species_3) |>
  summarise(
    lowest_upper_bound = min(q97.5),
    highest_lower_bound = max(q2.5),
    is_distinct = highest_lower_bound > lowest_upper_bound) |>
  filter(is_distinct) |>
  pull(species_3)

# View the list of species
print(significant_species2)


# Create a darker background for species that have different T50 among sites
# based on the previous model:

# Ensure significant_species is a character vector
significant_species2 <- as.character(significant_species2)

# Create the highlight dataframe
highlight_df2 <- data.frame(species_3 = significant_species2, 
                            start = -Inf, 
                            end = Inf, 
                            bottom = -Inf, 
                            top = Inf)

# Plot
figS15 <- ggplot(data = dt_aux3,
                 aes(x = site, y = median)) +
  geom_rect(data = highlight_df2, 
            aes(xmin = start, xmax = end, 
                ymin = bottom, ymax = top), 
            fill = "grey80", 
            alpha = 0.5,
            inherit.aes = FALSE) +
  geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.25, 
                                        linetype = 1)) +
  labs(x = expression("Sites"), 
       y = expression("Estimated T"[Crit]~"(°C)")) +
  facet_wrap(~ species_3, 
             ncol = 4,
             scales = "free",
             labeller = as_labeller(spp_rep_names)) +
  facetted_pos_scales( y = list(
    !COL == 1 ~ scale_y_continuous(labels = NULL, 
                                   limits = c(40, 55), 
                                   breaks = seq(40, 55, by = 5)),
    COL == 1 ~ scale_y_continuous(limits = c(40, 55), 
                                  breaks = seq(40, 55, by = 5))))

figS15

# END --------------------------------------------------------------------------
