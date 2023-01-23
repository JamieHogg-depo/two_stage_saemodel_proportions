##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                           STEP 1: GENERATE CENSUS                        ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @param grid data.frame of simulations parameters
#' @param QaS index from simulation for loop
generateCENSUS <- function(grid, QaS){
  
## PARAMS ## -------------------------------------------------------------------
seed <- grid$seed[QaS]
M <- 100                          			# Total number of areas
m <- grid$m[QaS]                  			# Number of areas to sample
N_lower <- 500 #450               			# population size per area
N_upper <- 3000 #1400
TP_lower <- grid$TP_lower[QaS]              # TP: True Proportion
TP_upper <- grid$TP_upper[QaS]            	# TP: True Proportion
alpha_survey <- grid$alpha_survey[QaS]    	 # magnitude of affect for unit-level survey covariate
alpha_nonsurvey <- grid$alpha_nonsurvey[QaS] # magnitude of affect for unit-level nonsurvey covariate
gamma <- grid$gamma[QaS]          			 # magnitude of affect for area-level covariates
SF <- 0.004                       			 # sampling fraction
inform <- grid$inform[QaS]        			 # whether sample is informative 

## Generated PARAMS ## ---------------------------------------------------------
set.seed(seed)
N_i <- round(runif(M, N_lower, N_upper))
TP <- seq(from = TP_lower, to = TP_upper, length.out = M)
N <- sum(N_i)

SPL <- list(M = M,
            m = m,
            alpha = alpha, 
            gamma = gamma,
            inform = inform,
            SF = SF,
            N = N, 
            seed = seed)

## STAGE 1: Outcome ## ---------------------------------------------------------

# --- Area level counts of outcome ---
ll <- list()
for(j in 1:M){
  temp <- rbinom(1, N_i[j], TP[j])
  ll[[j]] <- data.frame(y = c(1,0), 
                        N = c(temp,N_i[j]-temp),
                        area = j)
}

# create census
census_raw <- bind_rows(ll) %>% 
  uncount(N)

# add sample sizes in each area
fixed_area_ss <- data.frame(area = 1:M,
                            N_i = N_i,
                            # increase the area_ss based on the number of areas
                            # that will be sampled
                            area_ss = round((M/m)*SF * N_i))
SPL$n <- sum(fixed_area_ss$area_ss)
SPL$nbar <- median(fixed_area_ss$area_ss)
SPL$nmax <- max(fixed_area_ss$area_ss)
SPL$nmin <- min(fixed_area_ss$area_ss)

# return population size
N <- nrow(census_raw)

## STAGE 2: Unit-level covariates ## -------------------------------------------

## --- Three unit level covariates ---
y <- census_raw$y
# get residuals
e1 <- rnorm(N, 0, 1)
e2 <- rnorm(N, 0, 1)

# generate continuous covariates
x1 <- y + alpha_survey*e1
census_raw$x2 <- scale(y + alpha_nonsurvey*e2)[,1]

# convert to categorical
census_raw$x1 <- quant_groups(x1, 3)

#remove objects
rm(e1, e2, x1, y)

## STAGE 3: Area-Level Covariates ## -------------------------------------------

## --- Area-level aggregate dataset ---
Ys <- aggregate(y~area, data = census_raw, mean)$y
# get residuals
v1 <- rnorm(M, 0, 1)

# generate continuous covariates
k1 <- scale(jlogit(Ys) + gamma * v1)[,1]

# create census_agg
census_raw_agg <- census_raw %>% 
  group_by(area) %>% 
  summarise(Ys = mean(y),
            N_i = n()) %>% 
  mutate(k1 = k1, 
         pi_area = N_i/SPL$N)

# add to census
census_raw <- census_raw_agg %>% 
  dplyr::select(area, k1) %>% 
  right_join(.,census_raw, by = "area")

# remove objects
rm(Ys, k1, v1, ll)

## STAGE 4: Add z metric and calculate sample probabilities and weights --------

census_raw <- census_raw %>% 
  #mutate(z = inform * y + 0.5 * rexp(N)) %>% 
  mutate(z = inform * ifelse(y == 1, 0, 1) + 0.8 * rexp(N)) %>% 
  left_join(.,fixed_area_ss, by = "area") %>% 
  group_by(area) %>% 
  mutate(area_ss = mean(area_ss),
         sum_z = sum(z)) %>% 
  ungroup() %>% 
  mutate(pi_ij = z/sum_z,
         w = (1/area_ss)*(1/pi_ij))

## STAGE 4.1: Output list ## ---------------------------------------------------
rm(j, N_upper, N_lower, TP, TP_lower, TP_upper, temp,
   seed, N_i, N)

out_list <- list()
out_list$census <- census_raw
out_list$census_agg <- census_raw_agg
out_list$SPL <- SPL
out_list$fixed_area_ss <- fixed_area_ss
return(out_list)

}

## ---- END ---- ##