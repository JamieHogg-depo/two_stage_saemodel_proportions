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

setwd("C:/r_proj/two_stage_saemodel_proportions")

#cur_date <- "20230201"

# pr_all <- readRDS(paste0("ResultsHPC/", cur_date, "/pr_all.rds"))
# sim_list <- readRDS(paste0("ResultsHPC/", cur_date, "/sim_list.rds"))
# loc_plots <- paste0("ResultsHPC/", cur_date, "/plots")
# loc_sheets <- paste0("ResultsHPC/", cur_date, "/sheets")
pr_all <- readRDS("data/pr_all.rds")
sim_list <- readRDS("data/sim_list.rds")
loc_plots <- paste0("simulation_study/plots")
loc_sheets <- paste0("simulation_study/sheets")

# Load modeled results
spm_pa <- readRDS("data/20231016/spm_pa.rds")
spm_global <- readRDS("data/20231016/spm_global.rds")
fp_global <- readRDS("data/20231016/fp_global.rds")

# Load pr_all files
pr <- lapply(list.files("Z:/paper1/outputs/202310160", pattern = "pr_all", full.names = T), readRDS)
names(pr) <- list.files("Z:/paper1/outputs/202310160", pattern = "pr_all")

# Load sim_list files
sim_list <- lapply(list.files("Z:/paper1/outputs/202310160", pattern = "sim_list", full.names = T), readRDS)
names(sim_list) <- list.files("Z:/paper1/outputs/202310160", pattern = "sim_list")

## Temporary renaming functions ## ---------------------------------------------

filter_w_tsln <- function(.data){filter(.data, model %in% c("BETA", "s2LN", 
                                                            "BIN", "LOG", 
                                                            "ELN", "HT_Direct"))}

rename_w_tsln <- function(.data){
  mutate(.data, 
         model = ifelse(model == "s2LN", "TSLN", model),
         model = ifelse(model == "HT_Direct", "HT", model))
}

make_numeric_decimal <- function(.data){
  df <- .data
  cols_to_format <- unlist(lapply(df, is.numeric))
  df[,cols_to_format] <- bind_cols(lapply(df[,cols_to_format], sprintf, fmt = '%#.2f'))
  return(df)
}

getMedIQR_aschar <- function(x, digits = 2, na.rm = TRUE){
  fmt = paste0('%#.', digits, "f")
  if(digits == 0) fmt = '%0.0f'
  med <- sprintf(median(x, na.rm = na.rm), fmt = fmt)
  lower <- sprintf(unname(quantile(x, probs = 0.25, na.rm = na.rm)), fmt = fmt)
  upper <- sprintf(unname(quantile(x, probs = 0.75, na.rm = na.rm)), fmt = fmt)
  paste0(med, " (", lower, ", ", upper, ")")
}

## Table 4 - summary of simulated data ## --------------------------------------

# Create objects for table
wols <- bind_rows(lapply(1:5, FUN = function(x)data.frame(WOLSB = unlist(pr[[x]]$WOLSB)) %>% 
                           mutate(QaS = as.character(sort(rep(1:6, 100)))))) %>% 
  group_by(QaS) %>% 
  summarise(Med_WOLSB = getMedIQR_aschar(WOLSB))
samp_var_sm <- bind_rows(lapply(1:5, FUN = function(x)bind_rows(pr[[x]]$samp_var_sm, .id = "QaS"))) %>% 
  group_by(QaS) %>% 
  summarise(Perc_increas_samp = getMedIQR_aschar(100*(ratio-1), digits = 0))
bias_red <- bind_rows(lapply(1:5, FUN = function(x)bind_rows(pr[[x]]$bias_red, .id = "QaS"))) %>% 
  group_by(QaS) %>% 
  summarise(Med_bias_red = getMedIQR_aschar(100*(1-ratio), digits = 0))
insta <- spm_pa %>% 
  filter(!is.na(prop_samp)) %>% 
  mutate(yes = ifelse(is.na(HT_el_VAR), 1, 0)) %>% 
  group_by(QaS, rep_counter) %>% 
  summarise(prop_unstable = 100*mean(yes),
            .groups = "drop") %>% 
  group_by(QaS) %>% 
  summarise(Mean_PropUnstable = getMedIQR_aschar(prop_unstable, digits = 1))

# create full table and save as excel
list(insta, samp_var_sm, bias_red, wols) %>% 
  reduce(left_join, by = "QaS") %>% 
  dplyr::select(QaS, Mean_PropUnstable,
                Med_WOLSB, Perc_increas_samp, Med_bias_red) %>% 
  mutate(QaS = paste0("Sc", QaS)) %>% 
  rename(Scenario = QaS,
         "Percent (%) of unstable areas" = Mean_PropUnstable,
         "ALC" = Med_WOLSB,
         "Percent (%) increase in sampling variance" = Perc_increas_samp,
         "Percent (%) reduction in bias" = Med_bias_red) %>% 
  write_excel_csv(., file = paste0(loc_sheets, "/sr.csv"))

## Table 5 - Bayesian MRRMSE and MARB, ci size and coverage ## -----------------

# Get median 95 credible interval sizes
post_sd <- spm_pa %>% 
  #bind_rows(pr_all$spm_pa, .id="QaS") %>% 
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
freq_cov <- spm_pa %>% 
  #bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
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
spm_pa %>% 
  #bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
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
spm_global %>% 
  #bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter_w_tsln() %>% rename_w_tsln() %>% 
  group_by(QaS, model, missing) %>% 
  summarise(MRRMSE = round(median(RRMSE), 2),
            MARB = round(median(ARB), 2), 
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
  group_by(Scenario) %>% 
  mutate(MRF_rr = MRRMSE_FALSE / MRRMSE_FALSE[5],
         MAF_rr = MARB_FALSE / MARB_FALSE[5],
         MRT_rr = MRRMSE_TRUE / MRRMSE_TRUE[5],
         MAT_rr = MARB_TRUE / MARB_TRUE[5]) %>% 
  ungroup() %>% 
  mutate(MRRMSE_FALSE = paste0(sprintf(MRRMSE_FALSE, fmt = '%#.2f'), " (", sprintf(MRF_rr, fmt = '%#.2f'), ")"),
         MARB_FALSE = paste0(sprintf(MARB_FALSE, fmt = '%#.2f'), " (", sprintf(MAF_rr, fmt = '%#.2f'), ")"),
         MRRMSE_TRUE = paste0(sprintf(MRRMSE_TRUE, fmt = '%#.2f'), " (", sprintf(MRT_rr, fmt = '%#.2f'), ")"),
         MARB_TRUE = paste0(sprintf(MARB_TRUE, fmt = '%#.2f'), " (", sprintf(MAT_rr, fmt = '%#.2f'), ")"),
         CI_size_FALSE = sprintf(CI_size_FALSE, fmt = '%#.2f'),
         CI_size_TRUE = sprintf(CI_size_TRUE, fmt = '%#.2f'),
         Coverage_FALSE = sprintf(Coverage_FALSE, fmt = '%#.2f'),
         Coverage_TRUE = sprintf(Coverage_TRUE, fmt = '%#.2f'),
         Scenario = paste0("Sc", Scenario)) %>% 
  dplyr::select(-contains("_rr")) %>% 
  write_excel_csv(., file = paste0(loc_sheets, "/bay_metrics.csv"))

## Table 6 - Frequentist MSE ##-------------------------------------------------

# start table
start_table4 <-spm_global %>% 
  filter_w_tsln() %>% rename_w_tsln() %>% 
  group_by(QaS, model, missing) %>% 
  summarise(MRRMSE = median(RRMSE),
            MARB = median(ARB),
            Mcrpsprop = median(crpsprop),
            Mcrpslogit = median(crpslogit),
            .groups = "drop") %>% 
  arrange(QaS, missing, MRRMSE) %>% 
  pivot_wider(values_from = c(MARB, MRRMSE, Mcrpsprop, Mcrpslogit),
              names_from = missing) %>% 
  left_join(.,post_sd, by = c("QaS", "model")) %>% 
  left_join(.,freq_cov_wide, by = c("QaS", "model")) %>% 
  # put in correct order
  relocate(QaS, model, contains("_FALSE")) %>%
  #dplyr::select(-contains("logit")) %>% 
  dplyr::select(-contains("prop")) %>% 
  # arrange and rename 
  arrange(QaS, model) 

# start table
start_table6 <- fp_global %>% 
  filter_w_tsln() %>% 
  arrange(QaS, missing, MSE) %>% 
  group_by(QaS, missing) %>% 
  mutate(id = 1:length(unique(model))) %>% 
  ungroup() %>% 
  mutate(Bias = 100*(Bias^2),
         Variance = 100*Variance,
         MSE = 100*MSE) %>% 
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
  arrange(Scenario, model) %>% 
  left_join(.,dplyr::select(start_table4, QaS, model, Mcrpslogit_TRUE, Mcrpslogit_FALSE),
            by = c("model", "Scenario" = "QaS"))

# Get ratios
rankIT <- function(x){paste0("(", sprintf(x/x[6], fmt = '%#.2f'), ")")}
ratios <- start_table6 %>% 
  dplyr::select(-contains("Coverage")) %>% 
  group_by(Scenario) %>%
  mutate_at(vars(-c("Scenario", "model")), rankIT) %>% ungroup() %>% 
  setNames(c("Scenario", "model", paste0("Baseline", names(.)[-c(1,2)]))) %>% 
  dplyr::select(-c(Scenario, model))

# save to csv
start_table6 %>% 
  bind_cols(ratios) %>%
  dplyr::select(-contains(c("Bias", "Variance", "Coverage"))) %>%
  mutate(Scenario = paste0("Sc", Scenario)) %>% 
  dplyr::select(c("Scenario", "model",
                  "MSE_FALSE", "BaselineMSE_FALSE", "Mcrpslogit_FALSE", "BaselineMcrpslogit_FALSE",
                  "MSE_TRUE", "BaselineMSE_TRUE", "Mcrpslogit_TRUE", "BaselineMcrpslogit_TRUE")) %>% 
  mutate(model = factor(model, levels = c("HT", "BETA", "BIN", "ELN", "LOG", "TSLN")))  %>% 
  group_by(model) %>% mutate(id = cur_group_id()) %>% ungroup() %>% 
  mutate(Scenario = ifelse(id == 1, Scenario, "")) %>% dplyr::select(-id) %>% 
  make_numeric_decimal() %>% 
  setNames(str_remove_all(colnames(.), "_.*"))

## Table 2 (supplementary) - Model parameters ## ------------------------------

mpl <-  lapply(list.files("Z:/paper1/outputs/202310160/r/", pattern = "_mpl", full.names = T), readRDS)
names(mpl) <- list.files("Z:/paper1/outputs/202310160/r/", pattern = "_mpl")

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

spm_pa %>% 
  #bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
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

spm_global %>% 
  #bind_rows(pr_all$spm_global, .id = "QaS") %>% 
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
ggsave(paste0(loc_plots, "/MRRMSE.png"), width = 10, height = 8.35)

spm_global %>% 
  #bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter_w_tsln() %>% rename_w_tsln() %>% 
  rename(Scenario = QaS) %>% 
  mutate(cl = as.factor(ifelse(model == "TSLN", 1, 0.8))) %>% 
  mutate(missing = as.factor(ifelse(missing == FALSE, "Sampled areas", "Nonsampled areas")),
         missing = fct_relevel(missing, "Sampled areas"),
         Scenario = paste0("Sc", Scenario)) %>% 
  # filter out BIN and Sc2, Sc4 and Sc6
  filter(model != "BIN",
         !Scenario %in% c("Sc2", "Sc4", "Sc6")) %>% 
  mutate(Scenario = case_when(Scenario == "Sc1" ~ "50-50",
                              Scenario == "Sc3" ~ "Rare",
                              Scenario == "Sc5" ~ "Common")) %>% 
  # create plot
  ggplot(aes(y = RRMSE, x = model, fill = model, alpha = cl))+
  geom_boxplot()+theme_bw()+
  facet_grid(missing~Scenario, scales = "free_y")+
  labs(title = "",
       y = "Mean relative root mean square error (MRRMSE)",
       x = "")+
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = c("BETA", "ELN", "LOG", "**TSLN**"))+
  scale_fill_discrete()+
  theme(text=element_text(size=18), 
        axis.text.x = element_markdown(),
        legend.position = "null")
ggsave(paste0(loc_plots, "/sparse_MRRMSE.png"), width = 10, height = 8.35)

## Figure 2 - MARB ## ----------------------------------------------------------

spm_global %>% 
  #bind_rows(pr_all$spm_global, .id = "QaS") %>% 
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
ggsave(paste0(loc_plots, "/MARB_noBIN.png"), width = 10, height = 8.35)

spm_global %>% 
  #bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter(model %in% c("BETA", "s2LN", "LOG", "ELN")) %>% 
  mutate(model = ifelse(model == "s2LN", "TSLN", model)) %>% 
  rename(Scenario = QaS) %>% 
  mutate(cl = as.factor(ifelse(model == "TSLN", 1, 0.8))) %>% 
  mutate(missing = as.factor(ifelse(missing == FALSE, "Sampled areas", "Nonsampled areas")),
         missing = fct_relevel(missing, "Sampled areas"),
         Scenario = paste0("Sc", Scenario)) %>%
  # filter out BIN and Sc2, Sc4 and Sc6
  filter(model != "BIN",
         !Scenario %in% c("Sc2", "Sc4", "Sc6")) %>% 
  mutate(Scenario = case_when(Scenario == "Sc1" ~ "50-50",
                              Scenario == "Sc3" ~ "Rare",
                              Scenario == "Sc5" ~ "Common")) %>% 
  # create plot
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
ggsave(paste0(loc_plots, "/sparse_MARB_noBIN.png"), width = 10, height = 8.35)

