##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                              Normal model                                ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Model Code ## ---------------------------------------------------------------

MC <- "s2LN"

## Prepare data for stan model ## ----------------------------------------------

# randomly select direct estimate iterations from s1LN model
n_draws <- nrow(s1LN_its$theta_pd)
theta_pd_its <- s1LN_its$theta_pd[sample(n_draws, 
                                      round(n_draws * prop_draws)),]

# normal
z_a <- model.matrix(~ k1, 
                    data = ss_list$census_agg)
z_a <- as.matrix(scale(z_a, scale = F)[,-1])

# gamma design matrix
tp_f <- filter(ss_list$true_prop, !missing)
T_gamma <- model.matrix(~ log(area_ss), data = tp_f)

# for stan
s2LN_data <- list(q_a = 1,
                   z_a = z_a,
                  its = nrow(theta_pd_its),
                  theta_pd_its = theta_pd_its,
                  theta_pd_sd = apply(s1LN_its$theta_pd, 2, sd),
                  sqrt_gamma_input = apply(sqrt(s1LN_its$gamma_pd), 2, mean),
                  m = ss_list$SPL$m,
                  M = ss_list$SPL$M,
				  # gamma design matrix
                  q1_gvf = ncol(T_gamma),
                  T = T_gamma,
				  # other objects 
                  m_s = length(which(ss_list$true_prop$HT_stable == T)),
                  id_s = which(ss_list$true_prop$HT_stable == T),
                  id_us = as.array(which(ss_list$true_prop$HT_stable == F)),
				  anyunstable = nrow(filter(ss_list$true_prop, HT_stable)) != ss_list$SPL$m)

rm(z_a)

## Fit the stan Model ## -------------------------------------------------------

m_s <- Sys.time()
s2LN_fit <- sampling(NORM_comp,
                     data = s2LN_data,
                     chains = n.ch, 
                     iter = n.iter,
                     warmup = n.warm, 
                     refresh = 0)
s2LN_d <- round(as.numeric(Sys.time() - m_s, units = "mins"), 2) 

## Check HMC errors/performance ## ---------------------------------------------

# check HMC diagnostics
#check_hmc_diagnostics(s2LN_fit)

# Check Rhat's
s2LN_btess <- suppressMessages(get_bulktailESS(s2LN_fit))
s2LN_c <- cbind(as.data.frame(summary(s2LN_fit)$summary),
                s2LN_btess$df)
s2LN_div <- length(which(get_divergent_iterations(s2LN_fit)))
s2LN_its <- extract(s2LN_fit)

# check convergence for parameters of interest
s2LN_c_ip <- s2LN_c[str_detect(s2LN_c$parameter, "mu"),]
# stop the current repetition if convergence has not been achieved
if(mean(as.numeric(ifelse(s2LN_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE) != 1){
  message("Convergence was not achieved for model: ", MC, ". Proportion of important parameters with Rhat<1.02 is ", mean(as.numeric(ifelse(s2LN_c_ip$Rhat<1.02, 1, 0)), na.rm = TRUE), 
          ". Max Rhat is ", max(s2LN_c_ip$Rhat, na.rm = T), ".")
		  s2LN_convergence <- c(s2LN_convergence, 0)
  next
}
if(s2LN_div > 0){
  message("Convergence was not achieved for model: ", MC, ". There were ", s2LN_div, " divergent transition(s).")
		  s2LN_convergence <- c(s2LN_convergence, 0)
  next
}
s2LN_convergence <- c(s2LN_convergence, 1)

## Create "measures" object ## -------------------------------------------------

true_prop <- arrange(ss_list$true_prop, ps_area)

assign(paste0(MC, "_measures"), 
       s2LN_fit %>% 
         spread_draws(mu[ps_area]) %>%     
         # median and highest density continuous interval
         median_hdci(.simple_names = F) %>% 
         # rename variables
         rename_with(~gsub("mu", "s2LN", .x)) %>% 
         rename(s2LN.median = s2LN) %>% 
         # remove useless columns
         dplyr::select(-.width, -.point, -.interval) %>% 
         mutate(
           # Confidence interval size
           s2LN_ci_size = s2LN.upper - s2LN.lower,
           # relative root mean square error
           s2LN_RRMSE = getRRMSE(s2LN_its$mu, true_prop$prop),
           # variance of posterior
           s2LN_VAR = getVAR(s2LN_its$mu, true_prop$prop),
           # relative bias
           s2LN_ARB = getARB(s2LN_its$mu, true_prop$prop),
         )
)

# add to model parameters list - mpl
mpl$s2LN[[JaS]] <- s2LN_fit %>%
					gather_draws(lambda_a[i], sigma_v) %>%
					median_hdci() %>% 
					ungroup() %>% 
					mutate(model = MC)

# add to draws_u list
ar$s2LN[[JaS]] <- s2LN_fit %>% 
  spread_draws(mu[ps_area]) %>% 
  rename(s2LN = mu) %>% 
  left_join(.,dplyr::select(true_prop, ps_area, prop, area), by = "ps_area") %>%
  ungroup() %>% 
  dplyr::select(-ps_area) %>%
  group_by(area) %>% 
  mutate(JaS = JaS)
