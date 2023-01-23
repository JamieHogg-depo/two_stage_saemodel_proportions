##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                          Binomial Model                                  ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Model Code ## ---------------------------------------------------------------

MC <- "BIN"

## Prepare data for jags model ## ----------------------------------------------

  # model matrix
  z_a <- model.matrix(~ k1, data = ss_list$census_agg)
  z_a <- as.matrix(scale(z_a, scale = F)[,-1])
  
  # for jags
  BIN_data <- list(m = ss_list$SPL$m,
               M = ss_list$SPL$M,
               n_i = ss_list$true_prop$area_ss,
               N_i = ss_list$true_prop$area_pop,
               q_a = 1,
               z_a = z_a,
               sample_counts = ss_list$true_prop$y_count)
  rm(z_a)

## Create the JAGS Sampler ## --------------------------------------------------

  # Create sampler
  m_s <- Sys.time()
  BIN_comp <- jags.model(file = paste0(wd, "modelcode/BIN_jags.txt"), 
                         data = BIN_data, 
                         n.chains = n.ch,
                         n.adapt = n.warm*3,
						 quiet = TRUE)
  
  # sample from the model
  BIN_fit <- coda.samples(BIN_comp, 
                          variable.names = c("mu", "theta", "log_lik", 
                                             "sigma_v", "lambda_a"),
                          n.iter = n.iter*3 - n.warm*3)
  BIN_d <- round(as.numeric(Sys.time() - m_s, units = "mins"), 2)
  
## Check errors/performance ## -------------------------------------------------
  
  # Check Rhat's
  BIN_btess <- suppressMessages(get_bulktailESS(BIN_fit, Stan = FALSE))
  BIN_c <- cbind(MCMCsummary(BIN_fit),
                 BIN_btess$df)
  BIN_its <- as.data.frame(as.matrix(BIN_fit, chains=T))
  
  # check convergence for parameters of interest
  BIN_c_ip <- BIN_c[str_detect(BIN_c$parameter, c("mu|theta")),]
  # stop the current repetition if convergence has not been achieved
  if(mean(as.numeric(ifelse(BIN_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE) != 1){
    message("Convergence was not achieved for model: ", MC, ". Proportion of important parameters with Rhat<1.01 is ", mean(as.numeric(ifelse(BIN_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE), 
            ". Max Rhat is ", max(BIN_c_ip$Rhat, na.rm = T), ".")
    BIN_convergence <- c(BIN_convergence, 0)
    next
  }
  BIN_convergence <- c(BIN_convergence, 1)
  
  # get LOOCV
  BIN_log_lik <- as.matrix(BIN_its %>% dplyr::select(contains("log_lik")))
  BIN_reff <- relative_eff(exp(BIN_log_lik),
                           chain_id = c(rep(1,n.iter*3 - n.warm*3), 
                                        rep(2, n.iter*3 - n.warm*3),
                                        rep(3, n.iter*3 - n.warm*3),
                                        rep(4, n.iter*3 - n.warm*3)))
  BIN_loo <- loo(x = BIN_log_lik,
                 r_eff = BIN_reff)
  
## Create "measures" object ## -------------------------------------------------
  
true_prop <- arrange(ss_list$true_prop, ps_area)
  
# create separate iteration matrices
theta_BIN_its <- BIN_its %>% dplyr::select(contains("theta"))
mu_BIN_its <- BIN_its %>% dplyr::select(contains("mu"))
  
# create the measures object 
assign(paste0(MC, "_measures"), 
       BIN_fit %>% 
         spread_draws(mu[ps_area], # parameter of interest
                      theta[ps_area]) %>% 
         # mean and highest density continuous interval
         median_hdci(.simple_names = F) %>% 
         # rename variables
         rename_with(~gsub("theta", "BINpp", .x)) %>% 
         rename_with(~gsub("mu", "BIN", .x)) %>%
		 rename(BINpp.median = BINpp,
				BIN.median = BIN) %>%
         # remove useless columns
         dplyr::select(-.width, -.point, -.interval) %>% 
         mutate(
           # Confidence interval size
           BINpp_ci_size = BINpp.upper - BINpp.lower,
           BIN_ci_size = BIN.upper - BIN.lower,
           # relative root mean square error
           BIN_RRMSE = getRRMSE(mu_BIN_its, true_prop$prop),
           BINpp_RRMSE = getRRMSE(theta_BIN_its, true_prop$prop),
           # variance of posterior
           BIN_VAR = apply(mu_BIN_its, 2, var),
           BINpp_VAR = apply(theta_BIN_its, 2, var),
           # relative bias
           BIN_ARB = getARB(mu_BIN_its, true_prop$prop),
           BINpp_ARB = getARB(theta_BIN_its, true_prop$prop)
         ) 
)

# add to draws_u list
ar$BIN[[JaS]] <- BIN_fit %>% 
  spread_draws(mu[ps_area], # parameter of interest
               theta[ps_area]) %>% 
  rename(BINpp = theta,
         BIN = mu)%>% 
  left_join(.,dplyr::select(true_prop, ps_area, prop, area), by = "ps_area") %>% 
  ungroup() %>%
  dplyr::select(-ps_area) %>% 
  group_by(area) %>% 
  mutate(JaS = JaS)
  
# add to model parameters list - mpl
mpl$BIN[[JaS]] <- BIN_fit %>%
					gather_draws(lambda_a, sigma_v) %>%
					median_hdci() %>% 
					ungroup() %>% 
					mutate(model = MC)

# remove objects
rm(theta_BIN_its, mu_BIN_its)

