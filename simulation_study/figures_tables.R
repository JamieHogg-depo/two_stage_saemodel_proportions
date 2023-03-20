##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##            CREATE FIGURES AND TABLES FOR SIMULATION EXPERIMENT           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
library(tidyverse)
library(patchwork)
library(knitr)
library(readr)
library(ggtext)
library(latex2exp)

cur_date <- "20221129"

pr_all <- readRDS(paste0("ResultsHPC/", cur_date, "/pr_all.rds"))
sim_list <- readRDS(paste0("ResultsHPC/", cur_date, "/sim_list.rds"))
loc_plots <- paste0("ResultsHPC/", cur_date, "/plots")
loc_sheets <- paste0("ResultsHPC/", cur_date, "/sheets")

## Temporary renaming functions ## ---------------------------------------------

filter_w_tsln <- function(.data){filter(.data, model %in% c("BETA", "s2LN", 
                                                            "BIN", "LOG", 
                                                            "ELN", "HT_Direct"))}

rename_w_tsln <- function(.data){
  mutate(.data, 
         model = ifelse(model == "s2LN", "TSLN", model),
         model = ifelse(model == "HT_Direct", "HT", model))
}

## Table 4 - summary of simulated data ## --------------------------------------

# Create objects for table
sr <- bind_rows(pr_all$SR, .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(Med_SR = median(SR)) %>% 
  mutate_if(is.numeric, round, digits = 2)
wols <- data.frame(WOLSB = unlist(pr_all$WOLSB)) %>% 
  mutate(QaS = as.character(sort(rep(1:6, 100)))) %>% 
  group_by(QaS) %>% 
  summarise(Med_WOLSB = median(WOLSB)) %>% 
  mutate_if(is.numeric, round, digits = 2)
wcor <- data.frame(WCOR = unlist(pr_all$WCOR)) %>% 
  mutate(QaS = as.character(sort(rep(1:6, 100)))) %>% 
  group_by(QaS) %>% 
  summarise(Med_WCOR = median(WCOR)) %>% 
  mutate_if(is.numeric, round, digits = 2)
samp_var_sm <- bind_rows(pr_all$samp_var_sm, .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(Perc_increas_samp = median(100*(ratio-1), na.rm = T)) %>% 
  mutate_if(is.numeric, round, digits = 1)
bias_red <- bind_rows(pr_all$bias_red, .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(Med_bias_red = median(100*(1-ratio), na.rm = T)) %>% 
  mutate_if(is.numeric, round, digits = 1)
insta <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  filter(!is.na(prop_samp)) %>% 
  mutate(yes = ifelse(is.na(HT_el_VAR), 1, 0)) %>% 
  group_by(QaS, rep_counter) %>% 
  summarise(prop_unstable = 100*mean(yes),
            .groups = "drop") %>% 
  group_by(QaS) %>% 
  summarise(Mean_PropUnstable = median(prop_unstable)) %>% 
  mutate_if(is.numeric, round, digits = 1)
nbar <- data.frame(QaS = as.character(1:6), nbar = unlist(sim_list$nbarbar)) %>% 
  mutate_if(is.numeric, round, digits = 1)

# create full table and save as excel
list(nbar, insta, sr, samp_var_sm, bias_red, wcor, wols) %>% 
  reduce(left_join, by = "QaS") %>% 
  rename(Scenario = QaS,
         "Percent (%) of unstable areas" = Mean_PropUnstable,
         "Posterior medians of SR" = Med_SR, 
         "Weighted OLS" = Med_WOLSB,
         "Weighted Correlation" = Med_WCOR,
         "Percent (%) increase in sampling variance" = Perc_increas_samp,
         "Percent (%) reduction in bias" = Med_bias_red,
         "Area sample size" = nbar) #%>% 
  #write_excel_csv(., file = paste0(loc_sheets, "/sr.csv"))

## Table 5 - Bayesian MRRMSE and MARB, ci size and coverage ## -----------------

# Get median 95 credible interval sizes
post_sd <- bind_rows(pr_all$spm_pa, .id="QaS") %>% 
  dplyr::select(area, QaS, rep_counter, missing, contains("ci_size")) %>%
  rename_with(~gsub("_ci_size", "", .x)) %>% 
  pivot_longer(-c(area, QaS, missing, rep_counter)) %>% 
  group_by(QaS, name, missing) %>% 
  summarise(psd = median(value, na.rm = T), 
            .groups = "drop") %>% 
  rename(model = name) %>% 
  filter_w_tsln() %>% rename_w_tsln() %>% 
  pivot_wider(values_from = psd,
              names_from = missing) %>% 
  rename("CI_size_FALSE" = "FALSE",
         "CI_size_TRUE" = "TRUE") %>% 
  mutate_if(is.numeric, function(x){round(x, digits = 2)})

# frequentist coverage
freq_cov <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  dplyr::select(ends_with(c("upper", "lower")), area, prop, QaS, rep_counter, missing) %>% 
  pivot_longer(-c(area, prop, QaS, rep_counter, missing)) %>% 
  separate(name, c("model", "metric")) %>% 
  pivot_wider(names_from = metric,
              values_from = value) %>% 
  mutate(in_ci = ifelse(prop > lower & prop < upper, 1, 0)) %>% 
  group_by(QaS, model, missing) %>% 
  summarise(Coverage = round(mean(in_ci, na.rm = T), 2)) %>% 
  filter_w_tsln() %>% rename_w_tsln() 
freq_cov_wide <- freq_cov %>% 
  pivot_wider(values_from = Coverage,
              names_from = missing) %>% 
  rename("Coverage_FALSE" = "FALSE",
         "Coverage_TRUE" = "TRUE")

# Check what models give the closet coverage
bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  dplyr::select(ends_with(c("upper", "lower")), area, prop, QaS, rep_counter, missing) %>% 
  pivot_longer(-c(area, prop, QaS, rep_counter, missing)) %>% 
  separate(name, c("model", "metric")) %>% 
  pivot_wider(names_from = metric,
              values_from = value) %>% 
  mutate(in_ci = ifelse(prop > lower & prop < upper, 1, 0)) %>% 
  group_by(QaS, model, missing) %>% 
  summarise(Coverage = mean(in_ci, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(di = abs(Coverage - 0.95)) %>% 
  group_by(QaS, missing) %>% 
  slice_min(order_by = di, n = 1)

## table
bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter_w_tsln() %>% rename_w_tsln() %>% 
  group_by(QaS, model, missing) %>% 
  summarise(MRRMSE = round(median(RRMSE)*10, 2),
            MARB = round(median(ARB)*10, 2), 
            .groups = "drop") %>% 
  arrange(QaS, missing, MRRMSE) %>% 
  pivot_wider(values_from = c(MARB, MRRMSE),
              names_from = missing) %>% 
  left_join(.,post_sd, by = c("QaS", "model")) %>% 
  left_join(.,freq_cov_wide, by = c("QaS", "model")) %>% 
  relocate(QaS, model, MRRMSE_FALSE, MARB_FALSE, CI_size_FALSE, Coverage_FALSE, MRRMSE_TRUE, MARB_TRUE) %>%
  arrange(QaS, model) %>% 
  rename(Scenario = QaS, 
         Model = model) %>% 
  write_excel_csv(., file = paste0(loc_sheets, "/bay_metrics.csv"))

## Table 1 (supplementary) - Frequentist MSE ##--------------------------------

# create table
bind_rows(sim_list$fp_global, .id = "QaS") %>% 
  filter_w_tsln() %>% 
  arrange(QaS, missing, MSE) %>% 
  group_by(QaS, missing) %>% 
  mutate(id = 1:length(unique(model))) %>% 
  ungroup() %>% 
  mutate(Bias = 10000*(Bias^2),
         Variance = 10000*Variance,
         MSE = 10000*MSE) %>% 
  dplyr::select(missing, model, Bias, Variance, QaS, MSE) %>% 
  rename_w_tsln() %>% 
  arrange(QaS, missing, model) %>% 
  mutate(Bias = round(Bias, digits =2),
         Variance = round(Variance, digits =2),
         MSE = round(MSE, digits =2)) %>% 
  left_join(.,freq_cov, by = c("QaS", "missing", "model")) %>% 
  rename(Scenario = QaS) %>% 
  mutate(model = fct_relevel(as.factor(model), "HT")) %>% 
  pivot_wider(names_from = missing,
              values_from = c(Bias, Variance, MSE, Coverage)) %>% 
  relocate(Scenario, model, 
           Bias_FALSE, Variance_FALSE, MSE_FALSE, Coverage_FALSE) %>% 
  dplyr::select(-contains(c("Bias", "Variance"))) %>% 
  arrange(Scenario, model) #%>% 
  #write_excel_csv(., file = paste0(loc_sheets, "/freq_mse.csv"))

## Table 2 (supplementary) - Model parameters ## ------------------------------

bind_rows(sim_list$mpl, .id = "QaS") %>% 
  filter(!model %in% c("s1LN2", "s2LN2")) %>% 
  mutate(model = ifelse(model == "s2LN", "TSLN-S2", model),
         model = ifelse(model == "s1LN", "TSLN-S1", model)) %>% 
  mutate(median = round(median, 2)) %>% 
  pivot_wider(names_from = parameter,
              values_from = median) %>% 
  relocate(QaS, model, `beta_u[1]`, `beta_u[2]`,
           `beta_u[3]`, `beta_u[4]`, `sigma_e`) %>% 
  rename(Scenario = QaS, Model = model) #%>% 
  #write_excel_csv(., file = paste0(loc_sheets, "/coefs.csv"))

## Table 3 (NA) - Summaries of posterior medians of prevalences ## ------------

bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  dplyr::select(contains(".median"), QaS) %>% 
  pivot_longer(-QaS, names_to = "model") %>% 
  group_by(QaS, model) %>% 
  summarise(median = round(median(value), 2),
            min = sprintf(min(value), fmt = '%#.2f'),
            max = sprintf(max(value), fmt = '%#.2f')) %>% 
  mutate(model = str_replace(model, ".median", ""),
         out = paste0(median, " (", min, ", ", max, ")")) %>% 
  dplyr::select(-c(median, min, max)) %>% 
  filter(model %in% c("BETA", "BIN", "ELN", "LOG", "s2LN")) %>% rename_w_tsln() %>%
  write_excel_csv(., file = paste0(loc_sheets, "/summary_prev_ests.csv"))

## Figure 1 - MRRMSE ## --------------------------------------------------------

bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter_w_tsln() %>% rename_w_tsln() %>% 
  rename(Scenario = QaS) %>% 
  mutate(cl = as.factor(ifelse(model == "TSLN", 1, 0.8))) %>% 
  mutate(missing = as.factor(ifelse(missing == FALSE, "Sampled areas", "Nonsampled areas")),
         missing = fct_relevel(missing, "Sampled areas"),
         Scenario = paste0("Sc", Scenario)) %>% 
  ggplot(aes(y = RRMSE, x = model, fill = model, alpha = cl))+
  geom_boxplot()+theme_bw()+
  facet_grid(missing~Scenario, scales = "free_y")+
  labs(title = "",
       y = "Mean relative root mean square error (MRRMSE)",
       x = "")+
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = c("BETA", "BIN", "ELN", "LOG", "**TSLN**"))+
  scale_fill_discrete()+
  theme(text=element_text(size=18), 
        axis.text.x = element_markdown(),
        legend.position = "null")
#ggsave(paste0(loc_plots, "/MRRMSE.png"), width = 10, height = 8.35)

## Figure 2 - MARB ## ----------------------------------------------------------

bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter(model %in% c("BETA", "s2LN", "LOG", "ELN")) %>% 
  mutate(model = ifelse(model == "s2LN", "TSLN", model)) %>% 
  rename(Scenario = QaS) %>% 
  mutate(cl = as.factor(ifelse(model == "TSLN", 1, 0.8))) %>% 
  mutate(missing = as.factor(ifelse(missing == FALSE, "Sampled areas", "Nonsampled areas")),
         missing = fct_relevel(missing, "Sampled areas"),
         Scenario = paste0("Sc", Scenario)) %>%
  ggplot(aes(y = ARB, x = model, fill = model, alpha = cl))+
  geom_boxplot()+theme_bw()+
  facet_grid(missing~Scenario)+
  labs(title = "",
       y = "Mean absolute relative bias (MARB)",
       x = "")+
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = c("BETA", "ELN", "LOG", "**TSLN**"))+
  theme(text=element_text(size=18),
        axis.text.x = element_markdown(),
        legend.position = "null")
#ggsave(paste0(loc_plots, "/MARB_noBIN.png"), width = 10, height = 8.35)

