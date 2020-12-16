##### Name: GEV_full_model.R                                                              
##### Description: GEV Model 1 from paper
##### Purpose: Function to obtain MLEs of model 1 GEV parameters via nlm/optim  
########################################################################################

data <- split(annual_max_snowf_df,f=annual_max_snowf_df$name)

n.years <- nrow(annmax_midway)

year_return <- c(25,50,75,100)

year_set <- 1:n.years/10

snow_list_set <- list(data$Midway$ann_max_snowf,
                      data$Ohare$ann_max_snowf,
                      data$ParkForest$ann_max_snowf)

names(snow_list_set) <- names(data)

time_list_set <- list(year_set,year_set,year_set)

names(time_list_set) <- names(data)

model_full_gev_nlme <- function(snow_list,time_list,stationary = T){
  
  num_stations <- length(snow_list)
  
  if(!stationary){
    
    s1 <- snow_list$Midway
    
    s2 <- snow_list$Ohare
    
    s3 <- snow_list$ParkForest
    
    t1 <- time_list$Midway
    
    t2 <- time_list$Ohare
    
    t3 <- time_list$ParkForest 
    
  }else{
    
    s1 <- snow_list$Midway
    
    s2 <- snow_list$Ohare
    
    s3 <- snow_list$ParkForest
    
    
  }
  
  
  num_obs <- length(s1)
  
  gev_ll <- function(theta){
    
    mu1 <- theta[1]
    
    mu2 <- theta[2]
    
    mu3 <- theta[3]
    
    sc1 <- theta[4]
    
    sc2 <- theta[5]
    
    sc3 <- theta[6]
    
    xi1 <- theta[7]
    
    xi2 <- theta[8]
    
    xi3 <- theta[9]
    
    if(!stationary){
      
      mu_t <- theta[10]
      
    }
    
    if(stationary){
      
      y1 <- 1 + xi1 * (s1 - mu1)/sc1
      
      y2 <- 1 + xi2 * (s2 - mu2)/sc2
      
      y3 <- 1 + xi3 * (s3 - mu3)/sc3
      
    }else{
      
      y1 <- 1 + xi1 * (s1 - (mu1 + mu_t*t1))/sc1
      
      y2 <- 1 + xi2 * (s2 - (mu2 + mu_t*t2))/sc2
      
      y3 <- 1 + xi3 * (s3 - (mu3 + mu_t*t3))/sc3
      
    }
    
    
    if(any(y1) <= 0 || any(y2 <= 0) || any(y3 <= 0) || any(sc1 <= 0) || any(sc2 <= 0) || any(sc3 <= 0)){
      
      nllh <- 1e6
      
      
    }else{
      
      nllh <- num_obs*log(sc1) + num_obs*log(sc2) + num_obs*log(sc3) + sum(y1^(-1/xi1)) + sum(y2^(-1/xi2)) + sum(y3^(-1/xi3)) +
        sum(log(y1) * (1/xi1 + 1)) + sum(log(y2) * (1/xi2 + 1)) + sum(log(y3)*(1/xi3 + 1))
      
    }
    
    return(nllh)
    
  }  
  
  sig1_init <- sqrt(6 * var(s1))/pi
  
  sig2_init <- sqrt(6 * var(s2))/pi
  
  sig3_init <- sqrt(6 * var(s3))/pi
  
  mu1_init <- mean(s1) - 0.57722 * sig1_init
  
  mu2_init <- mean(s2) - 0.57722 * sig2_init
  
  mu3_init <- mean(s3) - 0.57722 * sig3_init 
  
  xi1_init <- 0.1
  
  xi2_init <- 0.1
  
  xi3_init <- 0.1
  
  mu_t_init <- 0
  
  if(stationary){
    
    init_par <- c(mu1_init,mu2_init,mu3_init,
                  sig1_init,sig2_init,sig3_init,
                  xi1_init,xi2_init,xi3_init) 
    
  }else{
    
    init_par <- c(mu1_init,mu2_init,mu3_init,
                  sig1_init,sig2_init,sig3_init,
                  xi1_init,xi2_init,xi3_init,mu_t_init) 
    
  }
  
  fit <- nlm(gev_ll,init_par,hessian = T)
  
  if(stationary){
    
    mu_t <- 0
    
  }else{
    
    mu_t <- fit$estimate[6]
    
  }
  
  return_mat <- matrix(0,num_stations,length(year_return))
  
  return_mat[1,] <- (fit$estimate[1] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[4]/fit$estimate[7])*(1-(-log(1-(1/year_return)))^(-fit$estimate[7])))
  
  return_mat[2,] <- (fit$estimate[2] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[5]/fit$estimate[8])*(1-(-log(1-(1/year_return)))^(-fit$estimate[8])))
  
  return_mat[3,] <- (fit$estimate[3] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[6]/fit$estimate[9])*(1-(-log(1-(1/year_return)))^(-fit$estimate[9])))
  
  return(list(returns=return_mat, trend=mu_t, output=fit))
  
}

stationary_mod_full <- model_full_gev_nlme(snow_list=snow_list_set)

nonstationary_mod_full <- model_full_gev_nlme(snow_list_set,time_list_set,stationary = F)

stationary_mod_full_AICc <- AICc(stationary_mod_full,snow_list_set)
nonstationary_mod_full_AICc <- AICc(nonstationary_mod_full,snow_list_set)

stationary_mod_full$output$estimate
nonstationary_mod_full$output$estimate

stationary_mod_full_AICc
nonstationary_mod_full_AICc
