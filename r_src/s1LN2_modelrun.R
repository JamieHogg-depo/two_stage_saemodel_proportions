##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                            Logistic model                                ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Model Code ## ---------------------------------------------------------------

MC <- "s1LN2"

## Prepare data for stan model ## ----------------------------------------------

# calculate overall weighted mean
HT_mu <- with(ss_list$sample, weighted.mean(y, w = w))
what_is_this <- ss_list$sample %>% 
  mutate(diff = w * (HT_mu - y)) %>% 
  group_by(ps_area) %>% 
  summarise(yes = abs(sum(diff)/sum(w))) %>% 
  summarise(bot_SR = sum(yes))
bot_SR <- what_is_this$bot_SR
rm(what_is_this, HT_mu)

# Logistic
x_u <- model.matrix(~ x1 + x2 + k1, 
                    data = ss_list$sample)
x_u <- scale(x_u, scale = F)
x_u <- x_u[,-1]

ni <- ss_list$sample_agg$area_ss
Ni <- ss_list$sample_agg$sum_w

# for stan
s1LN2_data <- list(n = nrow(ss_list$sample),
                 q_u = ncol(x_u),
                 x_u = x_u, 
                 y = ss_list$sample$y,
				 y2 = ss_list$sample$y,
                 w = ss_list$sample$w,
                 sum_w = ss_list$sample_agg$sum_w,
                 m = ss_list$SPL$m,
                 M = ss_list$SPL$M,
				 w_tilde = (ss_list$sample$w/sum(ss_list$sample$w)) * nrow(ss_list$sample),
				 ps_area = ss_list$sample$ps_area,
                 area_loc_id = first.changes(ss_list$sample$ps_area), 
                 area_loc_size = rle(ss_list$sample$ps_area)$lengths,
				 bot_SR = bot_SR, 
				 prior_scale_res = prior_scale_res, 
				 FPC = ( (Ni - ni) / (ni * Ni * (ni - 1)) ),
				 w_ss2 = ss_list$sample$w_ss2)

rm(x_u, ni, Ni)

## Fit the stan Model ## -------------------------------------------------------

m_s <- Sys.time()
s1LN2_fit <- sampling(s1LN2_comp,
                    data = s1LN2_data,
                    chains = n.ch, 
                    iter = n.iter,
                    warmup = n.warm,
					refresh = 0)
s1LN2_loo <- loo(s1LN2_fit)
# s1LN2_loo$estimates[1,1]
s1LN2_d <- round(as.numeric(Sys.time() - m_s, units = "mins"), 2) 

## Check HMC errors/performance ## ---------------------------------------------

# check HMC diagnostics
#check_hmc_diagnostics(s1LN2_fit)

# Check Rhat's
s1LN2_btess <- suppressMessages(get_bulktailESS(s1LN2_fit))
s1LN2_c <- cbind(as.data.frame(summary(s1LN2_fit)$summary),
                s1LN2_btess$df)
s1LN2_div <- length(which(get_divergent_iterations(s1LN2_fit)))
s1LN2_its <- extract(s1LN2_fit)

# check convergence for parameters of interest
s1LN2_c_ip <- s1LN2_c[str_detect(s1LN2_c$parameter, "mu_pd"),]
# stop the current repetition if convergence has not been achieved
if(mean(as.numeric(ifelse(s1LN2_c_ip$Rhat<1.01, 1, 0)), na.rm = TRUE) != 1){
  message("Convergence was not achieved for model: ", MC, ". Proportion of important parameters with Rhat<1.01 is ", mean(as.numeric(ifelse(s1LN2_c_ip$Rhat<1.01, 1, 0)), na.rm = TRUE), 
          ". Max Rhat is ", max(s1LN2_c_ip$Rhat, na.rm = T), ".")
		  s1LN2_convergence <- c(s1LN2_convergence, 0)
  next
}
if(s1LN2_div > 0){
  message("Convergence was not achieved for model: ", MC, ". There were ", s1LN2_div, " divergent transition(s).")
		  s1LN2_convergence <- c(s1LN2_convergence, 0)
  next
}
s1LN2_convergence <- c(s1LN2_convergence, 1)

## Create "measures" object ## -------------------------------------------------

true_prop <- arrange(ss_list$true_prop, ps_area)

assign(paste0(MC, "_measures"), 
       s1LN2_fit %>% 
         spread_draws(mu_pd[ps_area]) %>%     
         # median and highest density continuous interval
         median_hdci(.simple_names = F) %>% 
         # rename variables
         rename_with(~gsub("mu_pd", "s1LN2", .x)) %>% 
         rename(s1LN2.median = s1LN2) %>% 
         # remove useless columns
         dplyr::select(-.width, -.point, -.interval) %>% 
         right_join(.,dplyr::select(true_prop, ps_area), by = "ps_area") %>% 
         mutate(
           # Confidence interval size
           s1LN2_ci_size = s1LN2.upper - s1LN2.lower,
           # relative root mean square error
           s1LN2_RRMSE = getRRMSE(s1LN2_its$mu_pd, true_prop$prop),
           # variance of posterior
           s1LN2_VAR = getVAR(s1LN2_its$mu_pd, true_prop$prop),
           # relative bias
           s1LN2_ARB = getARB(s1LN2_its$mu_pd, true_prop$prop),
         )
)

# add to model parameters list - mpl
mpl$s1LN2[[JaS]] <- s1LN2_fit %>%
					gather_draws(beta_u[i], sigma_e, sigma_res) %>%
					median_hdci() %>% 
					ungroup() %>% 
					mutate(model = MC)

# add to draws_u list
ar$s1LN2[[JaS]] <- s1LN2_fit %>% 
  spread_draws(mu_pd[ps_area]) %>% 
  rename(s1LN2 = mu_pd) %>% 
  left_join(.,dplyr::select(true_prop, ps_area, prop, area), by = "ps_area") %>%
  ungroup() %>% 
  dplyr::select(-ps_area) %>%
  group_by(area) %>% 
  mutate(JaS = JaS)
