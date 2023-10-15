##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                        SCRIPT CONTROL: GET RUNNING                       ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start_time_full <- Sys.time()
library(tidyverse)
library(rstan)
library(tidybayes)
library(MASS)
library(rjags)
library(MCMCvis)
library(loo)

wd <- "C:/r_proj/two_stage_saemodel_proportions/simulation_study/r_src/"
source(paste0(wd, "functions.R"))

## GRAND PARAMETERS ##
n.iter <- 4000
n.warm <- 2000 
n.ch <- 4
prop_draws <- 0.5
prior_scale_res <- 0.5

## ---- Define the simulation parameters ---- ##
sim_grid <- data.frame(scenario = c("1I","2I", "3I", "4I", "5I", "6I"),
						seed = c(45, 45, 47, 47, 48, 48),
						inform = rep(1, 6),
						alpha_survey = c(0.8, 1.5, 0.8, 1.5, 0.8, 1.5), 	# was normally 1.5
						alpha_nonsurvey = rep(1, 6), 						# was normally 1.5
						gamma = c(0.05, 0.05, 0.01, 0.01, 0.01, 0.01),
						m = rep(60, 6),
						TP_lower = c(0.35, 0.35, 0.1, 0.1, 0.6, 0.6),
						TP_upper = c(0.65, 0.65, 0.4, 0.4, 0.9, 0.9))
QaS <- which(sim_grid$scenario == scen)

## ---- Compile the stan models ---- ## 

NORM_comp <- stan_model(file=paste0(wd, "modelcode/NORM_stan.stan"))
LOG_comp <- stan_model(file=paste0(wd, "modelcode/LOG_stan.stan"))
BETA_comp <- stan_model(file=paste0(wd, "modelcode/BETA_stan.stan"))
ELN_comp <- stan_model(file=paste0(wd, "modelcode/ELN_stan.stan"))
s1LN_comp <- stan_model(file=paste0(wd, "modelcode/s1LN_stan.stan"))
s1LN2_comp <- stan_model(file=paste0(wd, "modelcode/s1LN2_stan.stan"))

# start the timer
start_time_full <- Sys.time()

# define grand list object - overall performance
# each element is a performance metric
# element 1 of each element is the first scenario performance metrics 
sim_list <- list(scenario = list(),
                 # FP: per area
                 fp_pa = list(),# FP: frequentist point estimates
                 # FP: global
                 fp_global = list(), 
                 # SPM: global measures
                 spm_global = list(), # SPM: summarize_point-estimates_median
                 # SPM: per area
                 spm_pa = list(),
                 # CD: global measures
                 cd_global = list(), # CD: combine draws
                 # CD: per area 
                 cd_pa = list(),
				 # median model coefficients
				 mpl = list(),
                 # median number of sampled areas
                 mbar = list(),
                 # median of median area sampled sizes
                 nbarbar = list(),
                 # Convergence checks
                 Rhat_okay = list(),
                 Rhat_max = list(),
                 num_divs_max = list(),
				 prop_conv = list(),
				 total_reps = list(),
				 # how many do the models take
                 median_duration = list())

# Load workhorse functions to create census and take sample
#source(paste0(wd, "generateCENSUS.R"))
#source(paste0(wd, "takeSAMPLE.R"))

# initial update
sim_list$scenario[[QaS]] <- sim_grid[QaS,]

## STEP 1: Setup census ## 
#GC <- generateCENSUS(sim_grid, QaS)
all_data_surveys <- readRDS(paste0(wd, "data/", scen, ".rds")) 

# Define list to store results from each rep
  # pr: per rep
  # 14 elements where each element will have `no_reps` number of elements
  pr <-  list(m = list(),               # nu. of sampled areas
              nbar = list(),            # median sampled size of areas
              SPL = list(),             # simulation parameters
			  LOG_c = list(),          	# ALL parameter summary matrix from LOG model
              BETA_c = list(),          # ALL parameter summary matrix from BETA model
              BIN_c = list(),           # ALL parameter summary matrix from BIN model
			  ELN_c = list(),           # ALL parameter summary matrix from ELN model
			  s1LN_c = list(),        	# ALL parameter summary matrix from s1LN model
			  #s1LN2_c = list(),        # ALL parameter summary matrix from s1LN2 model
			  s2LN_c = list(),       	# ALL parameter summary matrix from s2LN model
			  #s2LN2_c = list(),       	# ALL parameter summary matrix from s2LN2 model
              spm_pa = list(),          # Summarize Point-estimate Median (spm) per area
              spm_global = list(),      # global measures for rep
			  nom_ci = list(),          # nominal confidence interval coverage
			  btess = list(),           # proportion of parameters with good bulk and tail ESS
              Rhat_okay = list(),       # proportion of parameters which Rhat<1.01
              Rhat_max = list(),        # max Rhat across ALL parameters
              Rhat_okay_ip = list(),    # proportion of IMPORTANT parameters which Rhat<1.01
              num_divs = list(),        # the number of divergence iterations in Stan
              true_prop = list(),       # current true_prop dataset
              model_duration = list(),  # how long for the model to fit  
			  loo = list(),             # LOOCV for models
			  loo_prop_bad_k = list(),  # LOOCV diagnostic
			  SR = list(),				# Smoothing Ratio (SR)
			  #SR2 = list(),				# Smoothing Ratio (SR)
			  WOLSB = list(),			# Regression coefficient from WOLS of pd vs d
			  WCOR = list(),			# Weighted pearson correlation between pd and d
			  samp_var_sm = list(),		# dataset comparing the reduction in sampling variance under TSLN
			  bias_red = list(),		# comparing the absolute bias between the HT_Direct and P_o from TSLN
			  ratio_topten = list() )	# Median ratio of reduction in sampling variance for the largest HT_el_VAR values
	pr_all <- pr
  
# Define list to collect posterior draws across all reps
  # ar: all reps
  ar <- list(LOG = list(),
             BETA = list(),
             BIN = list(),
			 ELN = list(),
			 #s1LN2 = list(),
			 #s2LN2 = list(),
			 s1LN = list(),
			 s2LN = list())
			 
# Define list to collect posterior medians of variance terms 
# and fixed coefficients across all reps
  # mpl: model parameter list
  mpl <- list(LOG = list(),
              BETA = list(),
              BIN = list(),
			  ELN = list(),
			  #s1LN2 = list(),
			  #s2LN2 = list(),
			  s1LN = list(),
			  s2LN = list())

# Define convergence vectors
  # each vector will be varying lengths based on how many times we believe they
  # have converged during simulation
            BETA_convergence = as.numeric()
			LOG_convergence = as.numeric()
            BIN_convergence = as.numeric()
			ELN_convergence = as.numeric()
			#s1LN2_convergence = as.numeric()
			#s2LN2_convergence = as.numeric()
			s1LN_convergence = as.numeric()
			s2LN_convergence = as.numeric()
  
## ---- Start taking repeated sampled from census ---- ##
JaS <- 1
repeat{
  
message("Now starting rep ", JaS, " of ", no_reps)
  
## STEP 2: Take stratified sample ## 
#ss_list <- takeSAMPLE(GC, pert = 0.01)
ss_list <- all_data_surveys[[JaS]]
pr$true_prop[[JaS]] <- ss_list$true_prop
  
## STEP 3: Fit the models ##
source(paste0(wd, "BETA_modelrun.R"))
message("Finished model: ", MC)

#source(paste0(wd, "s1LN2_modelrun.R"))
#message("Finished model: ", MC)

#source(paste0(wd, "s2LN2_modelrun.R"))
#message("Finished model: ", MC)

source(paste0(wd, "ELN_modelrun.R"))
message("Finished model: ", MC)

source(paste0(wd, "BIN_modelrun.R"))
message("Finished model: ", MC)

source(paste0(wd, "LOG_modelrun.R"))
message("Finished model: ", MC)

source(paste0(wd, "s1LN_modelrun.R"))
message("Finished model: ", MC)

source(paste0(wd, "s2LN_modelrun.R"))
message("Finished model: ", MC)

## STEP 4: Add results to PER REP LIST ##
source(paste0(wd, "step4.R"))

## Stopping rule
if(JaS == no_reps){
  break
}

# Next JaS value
JaS <- JaS + 1

}

## STEP 6: Collapse repetitions ## 
source(paste0(wd, "step6.R"))

total_run_time <- as.numeric(Sys.time() - start_time_full, units = "mins")
message("Total run time was ", 
        round(total_run_time, 2), 
        " mins. Each repetition took ", round(total_run_time/no_reps, 2), 
        " mins on average.")

# Save output
saveRDS(sim_list, file = paste0("TSLN2/outputs/", cur_date, "/r/scen", scen, "_sim_list.rds"))
saveRDS(ar, file = paste0("TSLN2/outputs/", cur_date, "/r/scen", scen, "_ar.rds"))
saveRDS(pr, file = paste0("TSLN2/outputs/", cur_date, "/r/scen", scen, "_pr.rds"))
saveRDS(pr_all, file = paste0("TSLN2/outputs/", cur_date, "/r/scen", scen, "_pr_all.rds"))








