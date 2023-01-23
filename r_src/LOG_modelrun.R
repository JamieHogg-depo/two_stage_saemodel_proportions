##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                            Logistic model                                ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Model Code ## ---------------------------------------------------------------

MC <- "LOG"

## Prepare data for stan model ## ----------------------------------------------

# sample
x_u <- model.matrix(~ x2 + k1, 
                    data = ss_list$sample)
sample_colmeans <- colMeans(x_u)
x_u <- scale(x_u, scale = F)[,-1]

# census
x_u_census <- model.matrix(~ x2 + k1, 
                           data = ss_list$census)
x_u_census <- scale(x_u_census, center = sample_colmeans, scale = F)
x_u_census <- x_u_census[,-1]

# for stan
LOG_data <- list(n = nrow(ss_list$sample),
                  q_u = ncol(x_u),
                  x_u = x_u,
                  x_u_census = x_u_census,
                  y = ss_list$sample$y,
                  m = ss_list$SPL$m,
                  N = ss_list$SPL$N,
                  M = ss_list$SPL$M,
                  w_tilde = (ss_list$sample$w/sum(ss_list$sample$w)) * nrow(ss_list$sample),
                  ps_area = ss_list$sample$ps_area,
                  ps_area_census = ss_list$census$ps_area,
                  area_loc_id = first.changes(ss_list$census$ps_area), 
                  area_loc_size = rle(ss_list$census$ps_area)$lengths)

rm(x_u, x_u_census, sample_colmeans)

## Fit the stan Model ## -------------------------------------------------------

m_s <- Sys.time()
LOG_fit <- sampling(LOG_comp,
                    data = LOG_data,
                    chains = n.ch, 
                    iter = n.iter,
                    warmup = n.warm, 
					refresh = 0)
LOG_loo <- loo(LOG_fit)
# LOG_loo$estimates[1,1]
LOG_d <- round(as.numeric(Sys.time() - m_s, units = "mins"), 2) 

## Check HMC errors/performance ## ---------------------------------------------

# check HMC diagnostics
#check_hmc_diagnostics(LOG_fit)

# Check Rhat's
LOG_btess <- suppressMessages(get_bulktailESS(LOG_fit))
LOG_c <- cbind(as.data.frame(summary(LOG_fit)$summary),
               LOG_btess$df)
LOG_div <- length(which(get_divergent_iterations(LOG_fit)))
LOG_its <- rstan::extract(LOG_fit)

# check convergence for parameters of interest
LOG_c_ip <- LOG_c[str_detect(LOG_c$parameter, "mu"),]
# stop the current repetition if convergence has not been achieved
if(mean(as.numeric(ifelse(LOG_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE) != 1){
  message("Convergence was not achieved for model: ", MC, ". Proportion of important parameters with Rhat<1.01 is ", mean(as.numeric(ifelse(LOG_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE), 
          ". Max Rhat is ", max(LOG_c_ip$Rhat, na.rm = T), ".")
  LOG_convergence <- c(LOG_convergence, 0)
  next
}
if(LOG_div > 0){
  message("Convergence was not achieved for model: ", MC, ". There were ", LOG_div, " divergent transition(s).")
		  LOG_convergence <- c(LOG_convergence, 0)
  next
}
LOG_convergence <- c(LOG_convergence, 1)

## Create "measures" object ## -------------------------------------------------

true_prop <- arrange(ss_list$true_prop, ps_area)

assign(paste0(MC, "_measures"), 
       LOG_fit %>% 
         spread_draws(mu[ps_area]) %>%     
         # median and highest density continuous interval
         median_hdci(.simple_names = F) %>% 
         # rename variables
         rename_with(~gsub("mu", "LOG", .x)) %>% 
         rename(LOG.median = LOG) %>% 
         # remove useless columns
         dplyr::select(-.width, -.point, -.interval) %>% 
         right_join(.,dplyr::select(true_prop, ps_area), by = "ps_area") %>% 
         mutate(
           # Confidence interval size
           LOG_ci_size = LOG.upper - LOG.lower,
           # relative root mean square error
           LOG_RRMSE = getRRMSE(LOG_its$mu, true_prop$prop),
           # variance of posterior
           LOG_VAR = getVAR(LOG_its$mu, true_prop$prop),
           # relative bias
           LOG_ARB = getARB(LOG_its$mu, true_prop$prop),
         )
)

# add to model parameters list - mpl
mpl$LOG[[JaS]] <- LOG_fit %>%
					gather_draws(beta_u[i], sigma_e) %>%
					median_hdci() %>% 
					ungroup() %>% 
					mutate(model = MC)


# add to draws_u list
ar$LOG[[JaS]] <- LOG_fit %>% 
  spread_draws(mu[ps_area]) %>% 
  rename(LOG = mu) %>% 
  left_join(.,dplyr::select(true_prop, ps_area, prop, area), by = "ps_area") %>%
  ungroup() %>% 
  dplyr::select(-ps_area) %>%
  group_by(area) %>% 
  mutate(JaS = JaS)
