########################################################################################
##### Name: gev_bootstrap.R                                                              
##### Description: Code for bootstrapping return levels
##### Notes: The function here requires the use of parallel computation             
########################################################################################
library(quantreg)
library(foreach)
library(doParallel)


# Calculate bias correction constant zBC using and return bootstrap samples
# Returns bootstrap samples and zBC matrix

zBCcalc <- function(snowf_data,year_return = c(25,50,75,100),
                           B = 1e4,stationary = T,model = 1,original_model){
  
  ############################
  # Set up data
  ###########################
  snow_list_set <- snowf_data
  
  series_lengths <- unlist(lapply(snow_list_set,length))
  
  n_data_midway <- series_lengths['Midway'] 
  
  n_data_ohare <- series_lengths['Ohare'] 
  
  n_data_parkforest <- series_lengths['ParkForest'] 
  
  num_stations <- length(snow_list_set)
  
  num_return_years <- length(year_return)
  
  original_returns <- original_model$returns
  
  zBC <- matrix(0,num_stations,num_return_years)
  
  colnames(zBC) <- colnames(original_returns)
  
  rownames(zBC) <- rownames(original_returns)
  
  # Detect cores for parallel computation
  cores <- detectCores()
  
  cl <- makeCluster(cores[1]-1)
  
  registerDoParallel(cl)
  
  # Comute B bootstrap estimates for return levels
  
  boots <- foreach(i = 1:B) %dopar% {
    
    source("./code/gev_mle_functions.R") 
    
    set.seed(i)
    
    index_midway <- sample(1:n_data_midway,replace=T)
    
    index_ohare <- sample(1:n_data_ohare,replace=T)
    
    index_parkforest <- sample(1:n_data_parkforest,replace=T)
    
    # resampling snow data
    
    resample_midway <- snow_list_set$Midway[index_midway]
    
    resample_ohare <- snow_list_set$Ohare[index_ohare]
    
    resample_parkforest <- snow_list_set$ParkForest[index_parkforest]
    
    resample_list <- list(resample_midway,resample_ohare,resample_parkforest)
    
    names(resample_list) <- names(snow_list_set)
    
    # This resamples years for non-stationary GEV case
    
    boot_year_midway <- index_midway/10
    
    boot_year_ohare <- index_ohare/10
    
    boot_year_parkforest <- index_parkforest/10
    
    boot_list <- list(boot_year_midway,boot_year_ohare,boot_year_parkforest)
    
    names(boot_list) <- names(snow_list_set)
    
    
    # Calculating bootstrap return levels based on settings
    
    return_boot <- switch (model,
         model1_gev_nlme(resample_list,boot_list,stationary = stationary),
         model2_gev_nlme(resample_list,boot_list,stationary = stationary),
         model3_gev_nlme(resample_list,boot_list,stationary = stationary))
    
    
    #if(model == 1){
      
    #  return_boot <- model1_gev_nlme(resample_list,boot,stationary = stationary)
      
    #}else if(model == 2){
      
    #  return_boot <- model2_gev_nlme(resample_list,boot,stationary = stationary)
      
    #}else{
      
    #  return_boot <- model3_gev_nlme(resample_list,boot,stationary = stationary)
      
    #}
    
    

  }
  
  stopCluster(cl)

  # Matrix of bootstraps for each station
  
  return_boot_midway <- matrix(0,B,num_return_years)
  
  return_boot_ohare <- matrix(0,B,num_return_years)
  
  return_boot_parkforest <- matrix(0,B,num_return_years)
  
  trends_boot <- matrix(0,B,1)
  
  for(i in 1:B) {
    
    return_boot_midway[i,] <- boots[[i]]$returns[1,]
    
    return_boot_ohare[i,] <- boots[[i]]$returns[2,]
    
    return_boot_parkforest[i,] <- boots[[i]]$returns[3,]
    
    trends_boot[i] <- boots[[i]]$trend
  }
  
  ### Get zBC - bias correction constant
  
  
  for(i in 1:num_return_years){
    
    zBC[1,i] <- qnorm(sum(return_boot_midway[,i] < original_returns[1,i])/B,mean=0,sd=1)
    
    zBC[2,i] <- qnorm(sum(return_boot_ohare[,i] < original_returns[2,i])/B,mean=0,sd=1)
    
    zBC[3,i] <- qnorm(sum(return_boot_parkforest[,i] < original_returns[3,i])/B,mean=0,sd=1)
    
  }
  
  results <- list('Sample'=list('midway' = return_boot_midway,
                                'ohare' = return_boot_ohare,
                                'parkforest' = return_boot_parkforest),
                  'zBC'=zBC,
                  'trend_boot' = trends_boot)
  
  return(results)
  
} 


# Calculate acceleration constant cA

cAcalc <- function(snowf_data,year_return = c(25,50,75,100),
                   num_years = 61,stationary = T,model = 1){
  
  ##########################################
  # Set up constants and data
  ##########################################
  
  snow_list_set <- snowf_data
  
  series_lengths <- unlist(lapply(snow_list_set,length))
  
  n_data_midway <- series_lengths['Midway'] 
  
  n_data_ohare <- series_lengths['Ohare'] 
  
  n_data_parkforest <- series_lengths['ParkForest'] 
  
  num_stations <- length(snow_list_set)
  
  all_years <- 1:num_years/10
  
  num_return_years <- length(year_return)
  
  # Compute delete-1 jacknife return levels
  
  # Set up cores 
  
  cores <- detectCores()
  
  cl <- makeCluster(cores[1]-1)
  
  registerDoParallel(cl)
  
  
  jacks <- foreach(i = 1:num_years) %dopar% {
    
    source("./code/gev_mle_functions.R") 
    
    jack_sample_midway <- snow_list_set$Midway[-i]
    
    jack_sample_ohare <- snow_list_set$Ohare[-i]
    
    jack_sample_parkforest <- snow_list_set$ParkForest[-i]
    
    jack_sample_list <- list(jack_sample_midway,jack_sample_ohare,jack_sample_parkforest)
    
    if(!stationary){
      
      jack_year <- list(all_years[-i],all_years[-i],all_years[-i])
      
    }else{
      
      jack_year <- list(all_years,all_years,all_years)
      
    }
    
    names(jack_year) <- names(jack_sample_list) <- names(snow_list_set)
    
    # Calculating bootstrap return levels based on settings
    
    return_jack <- switch(model,
                          model1_gev_nlme(jack_sample_list,jack_year,stationary = stationary),
                          model2_gev_nlme(jack_sample_list,jack_year,stationary = stationary),
                          model3_gev_nlme(jack_sample_list,jack_year,stationary = stationary))
    
 
  }
  

  stopCluster(cl)
  
  return_jack_midway <- matrix(0,num_years,num_return_years)
  
  return_jack_ohare <- matrix(0,num_years,num_return_years)
  
  return_jack_parkforest <- matrix(0,num_years,num_return_years)
  
  for(i in 1:num_years){
    
    return_jack_midway[i,] <- jacks[[i]]$returns[1,]
    
    return_jack_ohare[i,] <- jacks[[i]]$returns[2,]
    
    return_jack_parkforest[i,] <- jacks[[i]]$returns[3,]
    
    ### Calculating acceleration constant cA
    cA <- matrix(0,num_stations,num_return_years)
    
    for(i in 1:num_return_years) {
      
      summation_midway <- return_jack_midway[,i] - mean(return_jack_midway[,i])
      
      summation_ohare <- return_jack_ohare[,i] - mean(return_jack_ohare[,i])
      
      summation_parkforest <- return_jack_parkforest[,i] - mean(return_jack_parkforest[,i])
      
      cA[1,i] <- (1/6)*(sum(summation_midway^3)/(sum(summation_midway^2)^(3/2)))
      
      cA[2,i] <- (1/6)*(sum(summation_ohare^3)/(sum(summation_ohare^2)^(3/2)))
      
      cA[3,i] <- (1/6)*(sum(summation_parkforest^3)/(sum(summation_parkforest^2)^(3/2)))
      
    }
  }

  return(cA)
}

# Calculate bootstrap CI for return levels

ReturnLvlBootCI <- function(samples,zBC,cA,significance = 0.95,
                            year_return = c(25,50,75,100)){
  
  Z <- qnorm(1-(1-significance)/2)
  
  lower <- zBC + (zBC-Z)/(1-(cA*(zBC-Z)))
  
  upper <- zBC + (zBC+Z)/(1-(cA*(zBC+Z)))
  
  quantile_BCa_midway <- cbind(pnorm(lower[1,]),pnorm(upper[1,]))
  
  quantile_BCa_ohare <- cbind(pnorm(lower[2,]),pnorm(upper[2,]))
  
  quantile_BCa_parkforest <- cbind(pnorm(lower[3,]),pnorm(upper[3,]))
  
  midway_ci_return <- matrix(NA,nrow = length(year_return),3)
  
  ohare_ci_return <- matrix(NA,nrow = length(year_return),3)
  
  parkforest_ci_return <- matrix(NA,nrow = length(year_return),3)
  
  for(K in 1:length(year_return)){
    
      # Returns for midway
      
      midway_ci_return[K,1] <- quantile(samples$midway[,K], quantile_BCa_midway[K,1])
      
      midway_ci_return[K,2] <- quantile(samples$midway[,K], 0.5)
      
      midway_ci_return[K,3] <- quantile(samples$midway[,K], quantile_BCa_midway[K,2])
      
      # Returns for ohare
      
      ohare_ci_return[K,1] <- quantile(samples$ohare[,K], quantile_BCa_ohare[K,1])
      
      ohare_ci_return[K,2] <- quantile(samples$ohare[,K], 0.5)
      
      ohare_ci_return[K,3] <- quantile(samples$ohare[,K], quantile_BCa_ohare[K,2])
      
      # Returns for parkforest
      
      parkforest_ci_return[K,1] <- quantile(samples$parkforest[,K], quantile_BCa_parkforest[K,1])
      
      parkforest_ci_return[K,2] <- quantile(samples$parkforest[,K], 0.5)
      
      parkforest_ci_return[K,3] <- quantile(samples$parkforest[,K], quantile_BCa_parkforest[K,2]) 
  }
  
  CI_list <- list('midway' = midway_ci_return,'ohare' = ohare_ci_return,'parkforest' = parkforest_ci_return)
    
  return(CI_list)
  
}
