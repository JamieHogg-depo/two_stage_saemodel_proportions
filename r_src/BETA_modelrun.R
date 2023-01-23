##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                               Beta model                                 ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Model Code ## ---------------------------------------------------------------

MC <- "BETA"

## Prepare data for stan model ## ----------------------------------------------

tp_f <- filter(ss_list$true_prop, !missing)

# model design matrix
z_a <- model.matrix(~ k1, data = ss_list$census_agg)
z_a <- as.matrix(scale(z_a, scale = F)[,-1])  # no intercept

# for stan
BETA_data <- list(mu_d_o = tp_f$HT_Direct,
                  # vector of stable gamma's
                  sqrt_psi_o = filter(tp_f, HT_stable)$HT_SD,
                  # number of areas with stable phi's
                  m_s = nrow(filter(tp_f, HT_stable)),
                  q_a = 1,
                  z_a = z_a,
                  m = ss_list$SPL$m,
                  M = ss_list$SPL$M,
                  m_s = length(with(tp_f, which(HT_stable))),
                  logNi = scale(log(ss_list$census_agg$N_i))[,1],
                  logit_mu_bounds = jlogit(c(0.05, 0.95)),
                  id_s = with(tp_f, which(HT_stable)),
                  id_us_miss = with(ss_list$true_prop, which(HT_stable == FALSE | is.na(HT_stable))))

## Fit the stan Model ## -------------------------------------------------------

m_s <- Sys.time()
BETA_fit <- sampling(BETA_comp,
                     data = BETA_data,
                     chains = n.ch, 
                     iter = n.iter,
                     warmup = n.warm, 
                     control = list(adapt_delta = 0.9), 
                     init = 0,
					 refresh = 0)
BETA_loo <- loo(BETA_fit)
# BETA_loo$estimates[1,1]
BETA_d <- round(as.numeric(Sys.time() - m_s, units = "mins"), 2)

## Check HMC errors/performance ## ---------------------------------------------

# check HMC diagnostics
# check_hmc_diagnostics(BETA_fit)

# Check Rhat's
BETA_btess <- suppressMessages(get_bulktailESS(BETA_fit))
BETA_c <- cbind(as.data.frame(summary(BETA_fit)$summary),
                BETA_btess$df)
BETA_div <- length(which(get_divergent_iterations(BETA_fit)))
BETA_its <- rstan::extract(BETA_fit)

# check convergence for parameters of interest
BETA_c_ip <- BETA_c[str_detect(BETA_c$parameter, "mu"),]
# stop the current repetition if convergence has not been achieved
if(mean(as.numeric(ifelse(BETA_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE) != 1){
  message("Convergence was not achieved for model: ", MC, ". Proportion of important parameters with Rhat<1.01 is ", mean(as.numeric(ifelse(BETA_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE), 
          ". Max Rhat is ", max(BETA_c_ip$Rhat, na.rm = T), ".")
		  BETA_convergence <- c(BETA_convergence, 0)
  next
}
if(BETA_div > 0){
  message("Convergence was not achieved for model: ", MC, ". There were ", BETA_div, " divergent transition(s).")
		  BETA_convergence <- c(BETA_convergence, 0)
  next
}
BETA_convergence <- c(BETA_convergence, 1)

## Create "measures" object ## -------------------------------------------------

true_prop <- arrange(ss_list$true_prop, ps_area)

mu_BETA_its <- BETA_its$mu

assign(paste0(MC, "_measures"), 
       BETA_fit %>% 
         spread_draws(mu[ps_area]) %>% 
         # mean and highest density continuous interval
         median_hdci(.simple_names = F) %>% 
         # rename variables
         rename_with(~gsub("mu", "BETA", .x)) %>%
         rename(BETA.median = BETA) %>% 
         # remove useless columns
         dplyr::select(-.width, -.point, -.interval) %>% 
         # create performance metrics
         mutate(
           # Confidence interval size
           BETA_ci_size = BETA.upper - BETA.lower,
           # relative root mean square error
           BETA_RRMSE = getRRMSE(mu_BETA_its, true_prop$prop),
           # variance of posterior
           BETA_VAR = apply(mu_BETA_its, 2, var),
           # relative bias
           BETA_ARB = getARB(mu_BETA_its, true_prop$prop),
)
)


# add to draws_u list
ar$BETA[[JaS]] <- BETA_fit %>% 
  spread_draws(mu[ps_area]) %>% 
  rename(BETA = mu) %>% 
  left_join(.,dplyr::select(true_prop, ps_area, prop, area), by = "ps_area") %>% 
  ungroup() %>% 
  dplyr::select(-ps_area) %>%
  group_by(area) %>% 
  mutate(JaS = JaS)
  
# add to model parameters list - mpl
mpl$BETA[[JaS]] <- BETA_fit %>%
					gather_draws(lambda_a[i], sigma_v) %>%
					median_hdci() %>% 
					ungroup() %>% 
					mutate(model = MC)

# remove objects
rm(z_a, mu_BETA_its, tp_f)