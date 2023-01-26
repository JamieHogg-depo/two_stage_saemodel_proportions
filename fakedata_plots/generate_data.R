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
         