##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                                  FUNCTIONS                               ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @param param posterior draws (its x areas)
#' @param truth vector of true proportions
getVAR <- function(param, truth){
  out <- rep(NA, length(truth))
  if(ncol(param) != length(truth)){
    out[1:ncol(param)] <- apply(param, 2, var)
    return(out)
  }
  return(apply(param, 2, var))
}

#' @param param posterior draws (its x areas)
#' @param truth vector of true proportions
getRRMSE <- function(param, truth){
  out <- rep(NA, length(truth))
  for(i in 1:ncol(param)){
    out[i] <- (sqrt(mean((truth[i]-param[,i])^2)))/truth[i]
  }
  return(out)
}

#' @param param posterior draws (its x areas)
#' @param truth vector of true proportions
getARB <- function(param, truth){
  out <- rep(NA, length(truth))
  for(i in 1:ncol(param)){
    out[i] <- abs((mean(truth[i]-param[,i]))/truth[i])
  }
  return(out)
}

# fp: frequentist point estimates
#' @param listofmeasures current pr$spm_pa object
#' @returns df with area indicator and each model fARB and fRRMSE
getfp <- function(listofmeasures){
  list(
    getfBP_s(HT_Direct, listofmeasures),
    getfBP_s(BETA.median, listofmeasures),
	getfBP_s(LOG.median, listofmeasures),
    getfBP_s(BIN.median, listofmeasures),
	getfBP_s(ELN.median, listofmeasures),
	#getfBP_s(s1LN2.median, listofmeasures),
	#getfBP_s(s2LN2.median, listofmeasures),
	getfBP_s(s1LN.median, listofmeasures),
	getfBP_s(s2LN.median, listofmeasures)) %>% 
    reduce(left_join, by = c("area", "missing", "n_reps") ) %>% 
    rename_with(~gsub(".median", "", .x))
}

#' @param .data data.frame with all iterations and true proportion
#' @param var1-var4 non-character vectors for the variables
tidyPERF <- function(.data, var1, var2 = NULL, var3 = NULL, var4 = NULL){
  .data %>% 
    group_by(area) %>% 
    summarise(vv1_r = sqrt(mean((prop-{{var1}})^2, na.rm = T)),
              vv2_r = sqrt(mean((prop-{{var2}})^2, na.rm = T)),
              vv3_r = sqrt(mean((prop-{{var3}})^2, na.rm = T)),
              vv4_r = sqrt(mean((prop-{{var4}})^2, na.rm = T)),
              vv1_a = mean(prop-{{var1}}, na.rm = T),
              vv2_a = mean(prop-{{var2}}, na.rm = T),
              vv3_a = mean(prop-{{var3}}, na.rm = T),
              vv4_a = mean(prop-{{var4}}, na.rm = T),
              prop = mean(prop),
              "{{var1}}_VAR" := var_j({{var1}}, na.rm = T),
              "{{var2}}_VAR" := var_j({{var2}}, na.rm = T),
              "{{var3}}_VAR" := var_j({{var3}}, na.rm = T),
              "{{var4}}_VAR" := var_j({{var4}}, na.rm = T)) %>% 
    mutate("{{var1}}_RRMSE" := vv1_r/prop,
           "{{var2}}_RRMSE" := vv2_r/prop,
           "{{var3}}_RRMSE" := vv3_r/prop,
           "{{var4}}_RRMSE" := vv4_r/prop,
           "{{var1}}_ARB" := abs(vv1_a/prop),
           "{{var2}}_ARB" := abs(vv2_a/prop),
           "{{var3}}_ARB" := abs(vv3_a/prop),
           "{{var4}}_ARB" := abs(vv4_a/prop)) %>% 
    dplyr::select(-prop, -contains("vv")) %>% 
    dplyr::select(-contains("NULL"))
}

#' @param x numeric vector
var_j <- function(x, na.rm = F){
  if(is.null(x)){
    return(NULL)
  } else{
    var(x, na.rm = na.rm)
  }
}

#' @param var non-character vector defining the parameter of interest
#' @param listofmeasures current pr$spm_pa object
#' @returns df with area indicator and various frequentist metrics
getfBP_s <- function(var, listofmeasures){
  bind_rows(listofmeasures) %>% 
    group_by(area) %>% 
    mutate(br = prop - {{var}} ) %>% 
    ungroup() %>% 
    dplyr::select(area, br, missing, prop, {{var}}) %>% 
    group_by(area, missing) %>% 
    summarise(n_reps = n(),
				prop = mean(prop, na.rm = T),
			  "{{var}}_mean" := mean( {{var}}, na.rm = TRUE ),
			  "{{var}}_sd" := sd( {{var}}, na.rm = TRUE ),
			  "{{var}}_upper" := mean( {{var}}, na.rm = TRUE ) + 1.96 * sd( {{var}}, na.rm = TRUE ),
			  "{{var}}_lower" := mean( {{var}}, na.rm = TRUE ) - 1.96 * sd( {{var}}, na.rm = TRUE ),
              "{{var}}_fRRMSE" := (sqrt(mean(br^2, na.rm = TRUE)))/prop,
              "{{var}}_fARB" := abs((mean(br, na.rm = TRUE))/prop),
			  "{{var}}_biasminus" := prop - mean( {{var}}, na.rm = TRUE ),
			  "{{var}}_varInFrac" := ( sd( {{var}}, na.rm = TRUE ) )^2, 
				.groups = "drop" ) %>% 
    dplyr::select( -prop )
}

#' @param d vector 
first.changes <- function(d) {
  p <- cumsum(rle(d)$lengths) + 1
  c(1, p[-length(p)])
}

#' @param p vector of probabilities
jlogit <- function(p){
  log(p/(1-p))
}

#' @param eta vector of unconstrained values
jinvlogit <- function(eta){
  exp(eta)/(1+exp(eta))
}

#' @param x numeric vector
#' @param groups integer value specifying the number of categories required
quant_groups <- function(x, groups, ...) {
  
  # Calculate quantiles
  quantiles <- quantile(x, probs = seq(0, 1, 1 / groups), na.rm = TRUE, ...)
  
  # Create quantile groups
  out <- cut(x, breaks = quantiles, include.lowest = TRUE, ...)
  
  # relabel outcome
  levels(out) <- as.character(1:groups)
  
  # return groups
  return(out)
  
}

#' @param fit fitted Stan or JAGS model
get_bulktailESS <- function(fit, chains = c(1,2,3,4), Stan = TRUE){
  
  # array of iterations
  param_array <- as.array(fit)
  
  # number of chains
  n__chains <- length(chains)
  
  # complete the for loop
  if(Stan == TRUE){
    df <- data.frame(parameter = names(param_array[1,1,]),
                     ess_bulk = as.numeric(NA),
                     ess_tail = as.numeric(NA),
                     Rhat_j = as.numeric(NA))
    for(i in 1:nrow(df)){
      df$ess_bulk[i] <- ess_bulk(param_array[,chains,i])
      df$ess_tail[i] <- ess_tail(param_array[,chains,i])
      df$Rhat_j[i] <- Rhat(param_array[,chains,i])
    }
  }else{
    df <- data.frame(parameter = attr(param_array, "dimnames")$var,
                     ess_bulk = as.numeric(NA),
                     ess_tail = as.numeric(NA),
                     Rhat_j = as.numeric(NA))
    for(i in 1:nrow(df)){
      df$ess_bulk[i] <- ess_bulk(param_array[,i,chains])
      df$ess_tail[i] <- ess_tail(param_array[,i,chains])
      df$Rhat_j[i] <- Rhat(param_array[,i,chains])
    }
  }
  
  # props with reasonable ESS
  prop_ESS <- data.frame(bulk_ESS_below = mean(df$ess_bulk < 100 * n__chains, na.rm = T),
                         tail_ESS_below = mean(df$ess_tail < 100 * n__chains, na.rm = T))
  
  
  # helpful message
  message(paste0("All ESS should be > ", 100 * n__chains, 
                 ". Bulk_ESS below threshold: ", 
                 round(100*mean(df$ess_bulk < 100 * n__chains, na.rm = T), 2), 
                 "%. Tail_ESS below threshold: ", 
                 round(100*mean(df$ess_tail < 100 * n__chains, na.rm = T), 2), 
                 "%. Rhat's above 1.01: ",
                 round(100*mean(df$Rhat_j > 1.01, na.rm = T), 2),
                 "%."))
  
  # return the dataset
  return(list(df = df, 
              prop_ESS = prop_ESS))
}

# Stolen from sae::direct
jdirect <- function(y,dom,sweight,domsize, data, replace=FALSE){
    
    result <- data.frame(Domain=0,SampSize=0,Direct=0,
						VAR=0,SD=0,CV=0)
    
    did     <- unique(dom)    # unique identifiers of domains
    Dsample <- length(did)    # number of domains in sample
    
    # Calculate HT direct estimator for sampled domains   
    nds      <-rep(0,Dsample)   # domain sample sizes
    dirds    <-rep(0,Dsample)   # domain direct estimators 
    vardirds <-rep(0,Dsample)   # variances of direct estimators
    
    for (d in 1:Dsample){
      yd       <- y[dom==did[d]]
      nds[d]   <- length(yd)

      sweightd <- sweight[dom==did[d]]
	  w <- nds[d] * ( sweightd / sum(sweightd) )
      domsized <- domsize[(domsize[,1]==did[d]),2]
      
      dirds[d] <- sum(yd*sweightd)/domsized  
        
      # Approximated unbiased estimator of variance of HT direct estimator
      vardirds[d]<-(1/nds[d]) * (1-(nds[d]/domsized)) * (1/(nds[d]-1)) * sum( w^2 * (yd - dirds[d])^2 )
      
    }
    
    dird    <- dirds
    vardird <- vardirds
    nd      <- nds
    
    # Percent coeficient of variation of direct estimator
    cvdird<-100*sqrt(vardird)/dird
    
    result   <- data.frame(Domain=did,SampSize=nd,Direct=dird,
							VAR=vardird, SD=sqrt(vardird), CV=cvdird)
    roworder <- order(result[,1])
    
    result   <- result[roworder,]
    return(result)
  }
