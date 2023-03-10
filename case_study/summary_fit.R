## -----------------------------------------------------------------------------
## summary_fit ## --------------------------------------------------------------
## -----------------------------------------------------------------------------

# Packages
library(tidyverse)
library(HDInterval)
rm(list = ls())

# Source functions
source('C:/r_proj/two_stage_saemodel_proportions/case_study/functions_ALL.R')

apa <- function(.data, .m_n = FALSE, model_name = NULL){
  if(.m_n){
    .data %>% 
      mutate(ps_area = 1:nrow(.),
             model = model_name)
  }else{
    .data %>% mutate(ps_area = 1:nrow(.))
  }
}

# Load data
raw_est <- readRDS("C:/r_proj/two_stage_saemodel_proportions/data/raw_est.rds")

# For loop over the benchmarked and non-benchmarked results
for(i in c("bench", "nobench")){
  mb_est <- readRDS(paste0("C:/r_proj/two_stage_saemodel_proportions/data/modelled_est_", i, ".rds"))

# Set HT and HT_mu
HT_mu <- 0.147
HT_odds <- 0.147/(1-0.147)

## Posterior draws ## ----

mb_est$draws_or <-
  list(TSLN = getORs(mb_est$draws_mu$TSLN, baseline_odds = HT_odds),
       ELN = getORs(mb_est$draws_mu$ELN, baseline_odds = HT_odds),
       LOG = getORs(mb_est$draws_mu$LOG, baseline_odds = HT_odds))

## Posterior summary ## ----

mb_est$summ_mu <- 
  bind_rows(
    list(TSLN = getMCMCsummary(mb_est$draws_mu$TSLN, model_name = "TSLN") %>% apa(),
         ELN = getMCMCsummary(mb_est$draws_mu$ELN, model_name = "ELN") %>% apa(),
         LOG = getMCMCsummary(mb_est$draws_mu$LOG, model_name = "LOG")%>% apa()))

mb_est$summ_or <- 
  bind_rows(
    list(TSLN = getMCMCsummary(mb_est$draws_or$TSLN, model_name = "TSLN") %>% apa(),
         ELN = getMCMCsummary(mb_est$draws_or$ELN, model_name = "ELN") %>% apa(),
         LOG = getMCMCsummary(mb_est$draws_or$LOG, model_name = "LOG") %>% apa()))

## Difference in posterior probability ## ----

mb_est$DPP_mu <-
  bind_rows(
    list(TSLN = bind_cols(getDPP(mb_est$draws_mu$TSLN, null_value = HT_mu)) %>% apa(.m_n = T, model_name = "TSLN"),
         ELN = bind_cols(getDPP(mb_est$draws_mu$ELN, null_value = HT_mu)) %>% apa(.m_n = T, model_name = "ELN"),
         LOG = bind_cols(getDPP(mb_est$draws_mu$LOG, null_value = HT_mu)) %>% apa(.m_n = T, model_name = "LOG")))

mb_est$DPP_or <-
  bind_rows(
    list(TSLN = bind_cols(getDPP(mb_est$draws_or$TSLN, null_value = 1)) %>% apa(.m_n = T, model_name = "TSLN"),
         ELN = bind_cols(getDPP(mb_est$draws_or$ELN, null_value = 1)) %>% apa(.m_n = T, model_name = "ELN"),
         LOG = bind_cols(getDPP(mb_est$draws_or$LOG, null_value = 1)) %>% apa(.m_n = T, model_name = "LOG")))

## Add concordance data ## ----
mb_est$area_concor <- raw_est$area_concor

## Add direct estimates ## ----
mb_est$sa2_direct <- raw_est$sa2_direct
mb_est$sa4_direct <- raw_est$sa4_direct

## Update sa4_
mb_est$summ_mu_sa4 <- mb_est$summ_mu_sa4 %>% 
  left_join(.,mb_est$sa4_direct, by = "ps_sa4")

## Save full object ## ----
mb_est$draws_mu <- NULL
mb_est$draws_mu_sa4 <- NULL
mb_est$draws_or <- NULL
saveRDS(mb_est, file = paste0("C:/r_proj/two_stage_saemodel_proportions/data/sf_list_", i, ".rds"))

}
    
## END SCRIPT ## ---------------------------------------------------------------
    