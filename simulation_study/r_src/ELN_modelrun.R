##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                        EMPIRICAL LOGIT NORMAL MODEL                      ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Model Code ## ---------------------------------------------------------------

MC <- "ELN"

## Prepare data for stan model ## ----------------------------------------------

tp_f <- filter(ss_list$true_prop, !missing)

# model design matrix
z_a <- model.matrix(~ k1, data = ss_list$census_agg)
z_a <- as.matrix(scale(z_a, scale = F)[,-1])  # no intercept

# gamma design matrix
T_gamma <- model.matrix(~ HT_Direct + log(area_ss), data = tp_f)

# for stan
ELN_data <- list(theta_o = tp_f$HT_el_Direct,
                  # vector of stable gamma's
                  sqrt_gamma_o = filter(tp_f, HT_stable)$HT_el_SD,
                  # number of areas with stable phi's
                  m_s = nrow(filter(tp_f, HT_stable)),
                  # model design matrix
                  q_a = 1,
                  z_a = z_a,
                  # gamma design matrix
                  q1_gvf = ncol(T_gamma),
                  T = T_gamma,
                  # other objects
                  m = ss_list$SPL$m,
                  M = ss_list$SPL$M,
				  id_s = with(tp_f, which(HT_stable)),
                  id_us = as.array(with(tp_f, which(!HT_stable))),
				  anyunstable = nrow(filter(tp_f, HT_stable)) != ss_list$SPL$m)

## Fit the stan Model ## -------------------------------------------------------

m_s <- Sys.time()
ELN_fit <- sampling(ELN_comp,
                     data = ELN_data,
                     chains = n.ch, 
                     iter = n.iter,
                     warmup = n.warm, 
					 refresh = 0)
ELN_loo <- loo(ELN_fit)
# ELN_loo$estimates[1,1]
ELN_d <- round(as.numeric(Sys.time() - m_s, units = "mins"), 2)

## Check HMC errors/performance ## ---------------------------------------------

# check HMC diagnostics
#check_hmc_diagnostics(ELN_fit)

# Check Rhat's
ELN_btess <- suppressMessages(get_bulktailESS(ELN_fit))
ELN_c <- cbind(as.data.frame(summary(ELN_fit)$summary),
                ELN_btess$df)
ELN_div <- length(which(get_divergent_iterations(ELN_fit)))
ELN_its <- rstan::extract(ELN_fit)

# check convergence for parameters of interest
ELN_c_ip <- ELN_c[str_detect(ELN_c$parameter, "mu"),]
# stop the current repetition if convergence has not been achieved
if(mean(as.numeric(ifelse(ELN_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE) != 1){
  message("Convergence was not achieved for model: ", MC, ". Proportion of important parameters with Rhat<1.01 is ", mean(as.numeric(ifelse(ELN_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE), 
          ". Max Rhat is ", max(ELN_c_ip$Rhat, na.rm = T), ".")
		  ELN_convergence <- c(ELN_convergence, 0)
  next
}
if(ELN_div > 0){
  message("Convergence was not achieved for model: ", MC, ". There were ", ELN_div, " divergent transition(s).")
		  ELN_convergence <- c(ELN_convergence, 0)
  next
}
ELN_convergence <- c(ELN_convergence, 1)

## Create "measures" object ## -------------------------------------------------

true_prop <- arrange(ss_list$true_prop, ps_area)

mu_ELN_its <- ELN_its$mu

assign(paste0(MC, "_measures"), 
       ELN_fit %>% 
         spread_draws(mu[ps_area]) %>% 
         # mean and highest density continuous interval
         median_hdci(.simple_names = F) %>% 
         # rename variables
         rename_with(~gsub("mu", "ELN", .x)) %>%
         rename(ELN.median = ELN) %>% 
         # remove useless columns
         dplyr::select(-.width, -.point, -.interval) %>% 
         # create performance metrics
         mutate(
           # Confidence interval size
           ELN_ci_size = ELN.upper - ELN.lower,
           # relative root mean square error
           ELN_RRMSE = getRRMSE(mu_ELN_its, true_prop$prop),
           # variance of posterior
           ELN_VAR = apply(mu_ELN_its, 2, var),
           # relative bias
           ELN_ARB = getARB(mu_ELN_its, true_prop$prop),
)
)


# add to draws_u list
ar$ELN[[JaS]] <- ELN_fit %>% 
  spread_draws(mu[ps_area]) %>% 
  rename(ELN = mu) %>% 
  left_join(.,dplyr::select(true_prop, ps_area, prop, area), by = "ps_area") %>% 
  ungroup() %>% 
  dplyr::select(-ps_area) %>%
  group_by(area) %>% 
  mutate(JaS = JaS)
  
# add to model parameters list - mpl
mpl$ELN[[JaS]] <- ELN_fit %>%
					gather_draws(lambda_a[i], sigma_v) %>%
					median_hdci() %>% 
					ungroup() %>% 
					mutate(model = MC)

# remove objects
rm(T_gamma, z_a, mu_ELN_its, tp_f)
