
# packages
library(tidyverse)

wd <- "C:/r_proj/two_stage_saemodel_proportions/simulation_study/r_src/"
source(paste0(wd, "functions.R"))
source(paste0(wd, "generateCENSUS.R"))
source(paste0(wd, "takeSAMPLE.R"))
nu_reps <- 10

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

## Setup list
main_list <- list()

## LOOP - over scenarios #### --------------------------------------------------
for(scen in c("1I","2I", "3I", "4I", "5I", "6I")){
  
message("Starting ", scen)

# select QaS
QaS <- which(sim_grid$scenario == scen)

# generate census
main_list[[scen]]$GC <- generateCENSUS(sim_grid, QaS)

## LOOP - over repetitions #### ------------------------------------------------

  pb <- txtProgressBar(min = 0, max = nu_reps, style = 3)
  for(i in 1:nu_reps){
    
    main_list[[scen]]$ss_list[[i]] <- takeSAMPLE(main_list[[scen]]$GC, pert = 0.001)
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  message(paste0("---- Saving survey's"))
  saveRDS(main_list[[scen]]$ss_list, paste0(wd, "data/", scen, ".rds"))
}

## END SCRIPT ## ---------------------------------------------------------------













