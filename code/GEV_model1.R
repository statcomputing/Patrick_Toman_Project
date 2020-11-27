########################################################################################
##### Name: GEV_model1.R                                                              
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

model1_gev_nlme <- function(snow_list,time_list,stationary = T){
  
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
  
  
  num_obs <- length(c(s1,s2,s3))
  
  gev_ll <- function(theta){
    
    mu1 <- theta[1]
    
    mu2 <- theta[2]
    
    mu3 <- theta[3]
    
    sc <- theta[4]
    
    xi <- theta[5]
    
    if(!stationary){
      
      mu_t <- theta[6]
      
    }
    
    if(stationary){
      
      y1 <- 1 + xi * (s1 - mu1)/sc
      
      y2 <- 1 + xi * (s2 - mu2)/sc
      
      y3 <- 1 + xi * (s3 - mu3)/sc
      
    }else{
      
      y1 <- 1 + xi * (s1 - (mu1 + mu_t*t1))/sc
      
      y2 <- 1 + xi * (s2 - (mu2 + mu_t*t2))/sc
      
      y3 <- 1 + xi * (s3 - (mu3 + mu_t*t3))/sc
      
    }
    
    
    if(any(y1) <= 0 || any(y2 <= 0) || any(y3 <= 0) || any(sc <= 0)){
      
      nllh <- 1e6
      
      
    }else{
      
      nllh <- num_obs*log(sc) + sum(y1^(-1/xi)) + sum(y2^(-1/xi)) + sum(y3^(-1/xi)) +
        sum(log(y1) * (1/xi + 1)) + sum(log(y2) * (1/xi + 1)) + sum(log(y3)*(1/xi + 1))
      
    }
    
    return(nllh)
    
  }  
  
  siginit <- sqrt(6 * var(c(s1,s2,s3)))/pi
  
  mu1_init <- mean(s1) - 0.57722 * siginit
  
  mu2_init <- mean(s2) - 0.57722 * siginit
  
  mu3_init <- mean(s3) - 0.57722 * siginit 
  
  xi_init <- 0.1
  
  mu_t_init <- 0
  
  if(stationary){
   
    init_par <- c(mu1_init,mu2_init,mu3_init,siginit,xi_init) 
    
  }else{
    
    init_par <- c(mu1_init,mu2_init,mu3_init,siginit,xi_init,mu_t_init)
    
  }
  
  fit <- nlm(gev_ll,init_par,hessian = T)
  
  if(stationary){
    
    mu_t <- 0
    
  }else{
    
   mu_t <- fit$estimate[6]
    
  }
  
  return_mat <- matrix(0,num_stations,length(year_return))
  
  return_mat[1,] <- (fit$estimate[1] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[4]/fit$estimate[5])*(1-(-log(1-(1/year_return)))^(-fit$estimate[5])))
  
  return_mat[2,] <- (fit$estimate[2] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[4]/fit$estimate[5])*(1-(-log(1-(1/year_return)))^(-fit$estimate[5])))
  
  return_mat[3,] <- (fit$estimate[3] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[4]/fit$estimate[5])*(1-(-log(1-(1/year_return)))^(-fit$estimate[5])))
  
  return(list(returns=return_mat, trend=mu_t, output=fit))

}

stationary_mod1 <- model1_gev_nlme(snow_list=snow_list_set)

nonstationary_mod1 <- model1_gev_nlme(snow_list_set,time_list_set,stationary = F)
