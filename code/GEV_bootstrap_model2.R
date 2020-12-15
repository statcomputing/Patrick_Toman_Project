########################################################################################
##### Name: gev_bootstrap.R                                                              
##### Description: 
##### Purpose: boostrap of return levels for gev return levels    
########################################################################################

########################################################################################
#### Bootsrap for Model 2
########################################################################################

library(quantreg)
library(foreach)
library(doParallel)

## setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)


### Paramter set up
gev <- T

stationary <- F

# Significance level for confidence interval
significance <- 0.95

iteration <- 10000

series_lengths <- unlist(lapply(snow_list_set,length))

n_data_midway <- series_lengths['Midway'] 

n_data_ohare <- series_lengths['Ohare'] 

n_data_parkforest <- series_lengths['ParkForest'] 

num_stations <- length(snow_list_set)

num_return_years <- length(year_return)

num_years <- length(snow_list_set$Midway)

### Moving Block Bootstrap resampling

moving_block <- function(size,num_total,seed) {
  
  set.seed(seed)
  
  num_blocks <- num_total - size + 1
  
  blocks_out <- floor(num_total/size)
  
  block_index <- sample(1:num_blocks,size=blocks_out,replace=T)
  
  index <- numeric()
  
  for(i in 1:blocks_out) {
    
    index[(1+size*(i-1)):(size*i)] <- seq(from=block_index[i],to=block_index[i]+size-1)
  }
  
  index
}

########################################################################################
##### Computing bootstrap return levels                                            #####
########################################################################################

### Computing Bootstrap Return Level Samples for B iterations for GEV
if(gev){
  
  boots <- foreach(i = 1:iteration,.packages="quantreg") %dopar% {
    
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
    if(stationary){
      
      return_boot <- model2_gev_nlme(resample_list,boot_list,stationary = T)
      
    }else{
      
      return_boot <- model2_gev_nlme(resample_list,boot_list,stationary = F) 
      
    }
    
  }
  
} 

# Distributing bootstrap output into each station

return_boot_midway <- matrix(0,iteration,num_return_years)

return_boot_ohare <- matrix(0,iteration,num_return_years)

return_boot_parkforest <- matrix(0,iteration,num_return_years)

trends_boot <- matrix(0,iteration,1)

for(i in 1:iteration) {
  
  return_boot_midway[i,] <- boots[[i]]$returns[1,]
  
  return_boot_ohare[i,] <- boots[[i]]$returns[2,]
  
  return_boot_parkforest[i,] <- boots[[i]]$returns[3,]
  
  trends_boot[i] <- boots[[i]]$trend
}

### Calculating bias-correction constant, zBC

zBC <- matrix(0,num_stations,num_return_years)

for(i in 1:num_return_years){
  
  zBC[1,i] <- qnorm(sum(return_boot_midway[,i] < stationary_mod2$returns[1,i])/iteration,mean=0,sd=1)
  
  zBC[2,i] <- qnorm(sum(return_boot_ohare[,i] < stationary_mod2$returns[2,i])/iteration,mean=0,sd=1)
  
  zBC[3,i] <- qnorm(sum(return_boot_parkforest[,i] < stationary_mod2$returns[3,i])/iteration,mean=0,sd=1)
  
}


all_years <- 1:num_years/10 


if(gev){
  
  jacks <- foreach(i = 1:num_years,.packages="quantreg") %dopar% {
    
    jack_sample_midway <- snow_list_set$Midway[-i]
    
    jack_sample_ohare <- snow_list_set$Ohare[-i]
    
    jack_sample_parkforest <- snow_list_set$ParkForest[-i]
    
    jack_sample_list <- list(jack_sample_midway,jack_sample_ohare,jack_sample_parkforest)
    
    jack_year <- list(all_years,all_years,all_years)
    
    names(jack_year) <- names(jack_sample_list) <- names(snow_list_set)
    
    if(stationary){
      
      return_jack <- model2_gev_nlme(jack_sample_list,jack_year,stationary = T)
      
    }else{
      
      return.jack <- model2_gev_nlme(jack_sample_list,jack_year,stationary = F)
      
    }
  }
} 



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


### Calculating Quantiles for BCa
Z <- qnorm(1-(1-significance)/2)

lower <- zBC + (zBC-Z)/(1-(cA*(zBC-Z)))

upper <- zBC + (zBC+Z)/(1-(cA*(zBC+Z)))

quantile_BCa_midway <- cbind(pnorm(lower[1,]),pnorm(upper[1,]))

quantile_BCa_ohare <- cbind(pnorm(lower[2,]),pnorm(upper[2,]))

quantile_BCa_parkforest <- cbind(pnorm(lower[3,]),pnorm(upper[3,]))


### Computing Confidence Intervals
CI <- matrix(0,15,7)

CI[,1] <- c(1,year_return,2,year_return,3,year_return)


# For each CI, left: (alpha/2)%, middle: median, right: (1-alpha/2)%

for(i in 1:num_return_years){
  
  # Results for Midway 
  
  CI[1+i,2] <- quantile(return_boot_midway[,i], (1-significance)/2)
  
  CI[1+i,3] <- quantile(return_boot_midway[,i], 0.5)
  
  CI[1+i,4] <- quantile(return_boot_midway[,i], 1-(1-significance)/2)
  
  CI[1+i,5] <- quantile(return_boot_midway[,i], quantile_BCa_midway[i,1])
  
  CI[1+i,6] <- quantile(return_boot_midway[,i], 0.5)
  
  CI[1+i,7] <- quantile(return_boot_midway[,i], quantile_BCa_midway[i,2])
  
  # Results for Ohare 
  
  CI[6+i,2] <- quantile(return_boot_ohare[,i], (1-significance)/2)
  
  CI[6+i,3] <- quantile(return_boot_ohare[,i], 0.5)
  
  CI[6+i,4] <- quantile(return_boot_ohare[,i], 1-(1-significance)/2)
  
  CI[6+i,5] <- quantile(return_boot_ohare[,i], quantile_BCa_ohare[i,1])
  
  CI[6+i,6] <- quantile(return_boot_ohare[,i], 0.5)
  
  CI[6+i,7] <- quantile(return_boot_ohare[,i], quantile_BCa_ohare[i,2])
  
  # Results for Park Forest 
  
  CI[11+i,2] <- quantile(return_boot_parkforest[,i], (1-significance)/2)
  
  CI[11+i,3] <- quantile(return_boot_parkforest[,i], 0.5)
  
  CI[11+i,4] <- quantile(return_boot_parkforest[,i], 1-(1-significance)/2)
  
  CI[11+i,5] <- quantile(return_boot_parkforest[,i], quantile_BCa_parkforest[i,1])
  
  CI[11+i,6] <- quantile(return_boot_parkforest[,i], 0.5)
  
  CI[11+i,7] <- quantile(return_boot_parkforest[,i], quantile_BCa_parkforest[i,2])
  
}

### Rounding calculated values to third decimals
CI <- round(CI,3)

### Stop the cluster
stopCluster(cl)


### Displaying bootstrapped linear trends
hist(trends_boot)
c(quantile(trends_boot, c((1-significance)/2,1-(1-significance)/2)))

CI_list <- list('Midway' = CI[2:5,], 'Ohare' = CI[7:10,], 'ParkForest' = CI[12:15,])

