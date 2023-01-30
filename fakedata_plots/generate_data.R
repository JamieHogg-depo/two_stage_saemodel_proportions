# Generate fake data for case study plots

rm(list = ls())
library(tidyverse)
library(scales)
library(sf)
library(MASS)
library(patchwork)
library(readr)
library(readxl)

# Load map
map_sa2 <- st_read("~/OneDrive - Queensland University of Technology/DataLab/Requested file load/2016/2016_SA2_Shape_min/2016_SA2_Shape_min.shp") %>%  
  filter(STATE_CODE %in% c("8", "1", "3", "2")) %>% 
  mutate(SA2 = as.numeric(SA2_MAIN16)) %>% 
  mutate(ps_area = 1:1695)

# Load remoteness
ra1 <- read_csv("~/OneDrive - Queensland University of Technology/DataLab/Requested file load/request 5/aust_ASGS2016_sa2_detls.csv") %>% dplyr::select(-1)
ra2 <- read_csv("~/OneDrive - Queensland University of Technology/DataLab/Requested file load/request 5/Australia_SA22016_geogdet.csv") %>% dplyr::select(-1)

ra <- ra2 %>% 
  mutate(ra_sa2 = factor(ra2$ra_cat, levels = c(0,1,2,3,4), labels = unique(ra1$ra_name)[c(2,1,3,5,4)])) %>% 
  rename(SA2 = sa2maincode) %>% 
  dplyr::select(SA2, ra_sa2) %>% 
  mutate(ra_sa2 = str_replace(ra_sa2, " Australia", ""),
         ra_sa2 = str_replace(ra_sa2, " of", "")) %>% 
  mutate(ra_sa2_3c = as.factor(case_when(
    ra_sa2 %in% c("Outer Regional", "Remote", "Very Remote") ~ "Outer regional to very remote",
    ra_sa2 == "Major Cities" ~ "Major Cities",
    ra_sa2 == "Inner Regional" ~ "Inner Regional" 
  ))) %>% 
  mutate(ra_sa2_3c = fct_relevel(ra_sa2_3c, "Major Cities"))
rm(ra1, ra2)

# Load IRSD values
irsd <- read_excel("C:/r_proj/two_stage_saemodel_proportions/fakedata_plots/irsd.xlsx") %>% 
  mutate(ABS_irsd_decile_nation = as.factor(ABS_irsd_decile_nation))

# Create aux
aux2 <- inner_join(ra, irsd, by = "SA2") %>% 
  mutate(ps_area = 1:nrow(.),
         N_persons = round(runif(nrow(.), 4500, 13500)))

## Generate summ_mu ## ---------------------------------------------------------
S = matrix(c(0.00387, 0.00371, 0.00131,
             0.00371, 0.00369, 0.00141,
             0.00131, 0.00141, 0.00315),nrow=3) #with covariance of 0.8
y = mvrnorm(1695,c(0.147, 0.165, 0.144),S)
y = ifelse(y < 0, 0, y)

summ_mu <- as.data.frame(y) %>% 
  mutate(ps_area = 1:1695) %>% 
  setNames(c("TSLN", "ELN", "LOG", "ps_area")) %>% 
  pivot_longer(-ps_area, names_to = "model", values_to = "median") %>% 
  mutate(sd = exp(-2.0267 + 0.71314 * log(median)),
         cisize = exp(-0.653 + 0.738 * log(median)),
         lower = median - (cisize/2),
         upper = median + (cisize/2))

sample_agg <- data.frame(ps_area = 1:1695,
                         HT = exp(0.2733 + 1.9373 * log(summ_mu[summ_mu$model == "TSLN",]$median))) %>% 
  mutate(HT = ifelse(ps_area %in% sample(1695, 0.6*1695), HT, NA))
direct_est <- sample_agg %>% 
  rename(median = HT) %>% 
  mutate(cisize = 2*1.96*exp(-0.80212 + 0.70103 * log(median)),
         model = "Direct")

## Generate summ_or ## ---------------------------------------------------------

S = matrix(c(0.48029, 0.40909, 0.09932, 
             0.40909, 0.35669, 0.09753, 
             0.09932, 0.09753, 0.24202),nrow=3) #with covariance of 0.8
y = mvrnorm(1695,c(1.092, 1.073, 1.010),S)
y = ifelse(y < 0, 0.01, y)

summ_or <- as.data.frame(y) %>% 
  mutate(ps_area = 1:1695) %>% 
  setNames(c("TSLN", "ELN", "LOG", "ps_area")) %>% 
  pivot_longer(-ps_area, names_to = "model", values_to = "median") %>% 
  mutate(sd = exp(-1.0338 + 1.11 * log(median)),
         cisize = exp(0.1799 + 0.9638 * log(median)),
         lower = median - (cisize/2),
         upper = median + (cisize/2),
         EP = 0.6573 + 0.1688 * sd,
         DPP = 2*abs(EP-0.5),
         DPP_sig = as.factor(ifelse(DPP > 0.6, 1, 0)))

## Generate sa4_summ ## -------------------------------------------------------

S = matrix(c(199.5, 184.7, 121.4, 56.6, 67.9, 82.6, 
             184.7, 174.5, 121.9, 54.5, 65.5, 84.5,
             121.4, 121.9, 128.7, 35.1, 46.8, 83.3,
             56.6, 54.4, 35.1, 51.2, 38.7, 55.0,
             67.9, 65.5, 46.8, 38.7, 45.9, 61.5, 
             82.6, 84.5, 83.3, 55.0, 61.5, 108.6),nrow=6)/100000
y = mvrnorm(65,c(0.15691957, 0.15769829, 0.15512497, 0.04294622, 0.06100198, 0.09384463),S)
y = ifelse(y < 0, 0.001, y)

summ_mu_sa4 <- as.data.frame(y) %>% 
  mutate(ps_sa4 = 1:65) %>% 
  setNames(c("median_TSLN", "median_ELN", "median_LOG", "cisize_TSLN", "cisize_ELN", "cisize_LOG", "ps_sa4")) %>% 
  pivot_longer(-ps_sa4, names_to = c("metric", "model"), names_sep = "_", values_to = "value") %>% 
  pivot_wider(names_from = "metric", values_from = "value") %>% 
  mutate(sd = (cisize/2)/1.96,
         lower = median - (cisize/2),
         upper = median + (cisize/2),
         HT = exp(0.30292 + 1.198 * log(median)),
         res = rnorm(65*3, 0.06736968, 0.02636784),
         HT_lower = HT - res, HT_upper = HT + res,
         RRMSE = exp(-2.15070 + 0.07933 * log(median)),
         ARB = exp(-1.3766 + 0.1908 * log(median))) %>% 
  dplyr::select(-res)
