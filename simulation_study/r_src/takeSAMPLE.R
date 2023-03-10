##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                       STEP 2: TAKE STRATIFIED SAMPLE                     ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @param input_list list output from generateCENSUS()
takeSAMPLE <- function(input_list, pert = 0.001){

# Extract objects from input list
census <- input_list$census
census_agg <- input_list$census_agg
fixed_area_ss <- input_list$fixed_area_ss
SPL <- input_list$SPL

## STAGE 5: Take sample ## -----------------------------------------------------

ll <- list()
selected_areas <- sample(1:SPL$M, SPL$m, prob = census_agg$pi_area)
for(i in 1:SPL$m){
  temp_p <- census[census$area == selected_areas[i],]$pi_ij
  area_ss <- fixed_area_ss[fixed_area_ss$area == selected_areas[i],]$area_ss
  sample_ids <- sample(1:length(temp_p), area_ss, prob = temp_p)
  ll[[i]] <- census[census$area == selected_areas[i],][sample_ids,]
}
sample <- bind_rows(ll) %>% 
  arrange(area)
  
SPL$n <- nrow(sample)

## STAGE 5.1: Add psueo area indicators ## -------------------------------------

ps_area_id <- data.frame(area = unique(sample$area),
                         ps_area = 1:length(unique(sample$area)))
temp2 <- which(!unique(census$area) %in% unique(sample$area))

if(nrow(ps_area_id) != input_list$SPL$M){
  temp2 <- which(!unique(census$area) %in% unique(sample$area))
  ps_area_id <- rbind(ps_area_id,
                      data.frame(area = temp2,
                                 ps_area = nrow(ps_area_id)+1:length(temp2)))
}

# update census
census <- census %>% 
  inner_join(.,ps_area_id, by = "area") %>% 
  arrange(ps_area)

# update census_agg
census_agg <- census_agg %>% 
  inner_join(.,ps_area_id, by = "area") %>% 
  arrange(ps_area)

# update sample
sample <- sample %>% 
  inner_join(.,ps_area_id, by = "area") %>% 
  arrange(ps_area)

## STAGE 6: reweight the sample weights to sum to area population
sample <- sample %>% 
  group_by(ps_area) %>% 
  mutate(area_sum_w = sum(w),
		 f_i = n()/sum(w)) %>% 
  ungroup() %>% 
  mutate(w_ss2 = (f_i * w)^2,
		 w = N_i * (w/area_sum_w)) %>% 
  dplyr::select(-sum_z, -pi_ij, -area_sum_w, -f_i)

## STAGE 6.1: Get the HT estimates ## ------------------------------------------

# population counts for sampled areas
pop_sampled <- sample %>% 
  group_by(ps_area) %>% 
  summarise(nn = n()) %>% 
  left_join(., census_agg, by = "ps_area") %>% 
  dplyr::select(ps_area, N_i)

# new dataset
direct <- jdirect(y = as.numeric(sample$y), 
                      dom = as.numeric(sample$ps_area),
                      domsize = as.matrix(pop_sampled[,1:2]), 
                      sweight = sample$w) %>% 
    mutate(ps_area = 1:length(unique(sample$ps_area)),
           el_VAR = VAR / ( (Direct^2) * (1 - Direct )^2 ),
           el_Direct = jlogit( Direct ),
           stable = (Direct != 0 & Direct != 1))

# ensure that all direct estimates are non-zero
# here we add some jitter
direct$Direct <- ifelse(direct$Direct<pert, pert, direct$Direct)
direct$Direct <- ifelse(direct$Direct>(1-pert), (1-pert), direct$Direct)
direct$el_Direct <- ifelse(is.infinite(direct$el_Direct), NA, direct$el_Direct)
direct$el_Direct <- ifelse(direct$Direct == (1-pert), NA, direct$el_Direct)

# ensure that all VAR's are less than 0.25 and not exactly zero
direct$VAR <- ifelse(direct$VAR>0.249, NA, direct$VAR)
direct$VAR <- ifelse(direct$VAR == 0, NA, direct$VAR)
direct$el_VAR <- ifelse(is.infinite(direct$el_VAR), NA, direct$el_VAR)

# set CV nan's to zero
direct$CV[is.nan(direct$CV)]<-0

# add el_SD
direct$el_SD <- sqrt(direct$el_VAR)
direct$el_Direct = jlogit( direct$Direct )

# define simplier direct dataframe
df2 <- direct %>% 
  dplyr::select(ps_area, Direct, SD, CV, VAR, el_Direct, el_VAR, el_SD, stable) %>% 
  setNames(c("ps_area", paste0("HT_", names(.)[-1])))

## STAGE 6.2: Create other objects ## ------------------------------------------

# Create sample_agg
sample_agg <- sample %>% 
  group_by(ps_area) %>% 
  summarise(area_ss = n(),
            y_count = sum(y),
            prop_samp = y_count/area_ss,
            prop_samp_w = sum(y*w)/sum(w),
            sum_w = sum(w))

# Create true_prop
true_prop <- census %>% 
  group_by(ps_area, area) %>% 
  summarise(area_pop = n(), 
            prop = sum(y)/n(),
            .groups = "drop") %>% 
  left_join(.,sample_agg, by = "ps_area") %>% 
  mutate(missing = is.na(prop_samp)) %>% 
  left_join(.,df2, by = "ps_area")

# replace NAs for y_count and area_ss with zeros
true_prop[is.na(true_prop$area_ss),]$area_ss <- 0
true_prop[is.na(true_prop$y_count),]$y_count <- 0

## STAGE 6.3: Output list ## ---------------------------------------------------

SPL$nbar = median(sample_agg$area_ss)
SPL$nmax = median(sample_agg$area_ss)
SPL$nmin = median(sample_agg$area_ss)

out_list <- list()
out_list$census <- census
out_list$census_agg <- census_agg
out_list$sample <- sample
out_list$sample_agg <- sample_agg
out_list$true_prop <- true_prop
out_list$SPL <- SPL
out_list$direct <- direct
return(out_list)

}

## ---- END ---- ##