########################################################################################
##### Name: GEV_model3.R                                                              
##### Description: GEV Model 3 from paper
##### Purpose: Function to obtain MLEs of model 2 GEV parameters via nlm/optim. In this  
#####          particular function we fit the model to all of the same variables
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

model3_gev_nlme <- function(snow_list,time_list,stationary = T){
  
  num_stations <- length(snow_list)
  
  #g1_name <- paste0(names(snow_list)[groups$g1[1]],'-',names(snow_list)[groups$g1[2]])
  
  #g2_name <- names(snow_list)[groups$g2]
  
  
  if(!stationary){
    
    s1 <- unlist(snow_list)
    
    t1 <- unlist(time_list)
    
    
  }else{
    
    
    s1 <- unlist(snow_list)
    
  }
  
  
  num_obs <- length(c(s1))
  
  gev_ll <- function(theta){
    
    mu1 <- theta[1]
    
    sc <- theta[2]
    
    xi <- theta[3]
    
    if(!stationary){
      
      mu_t <- theta[4]
      
    }
    
    if(stationary){
      
      y1 <- 1 + xi * (s1 - mu1)/sc
      
    }else{
      
      y1 <- 1 + xi * (s1 - (mu1 + mu_t*t1))/sc
      
    }
    
    
    if(any(y1) <= 0 || any(sc <= 0)){
      
      nllh <- 1e6
      
      
    }else{
      
      nllh <- num_obs*log(sc) + sum(y1^(-1/xi)) + sum(log(y1) * (1/xi + 1))  
      
    }
    
    return(nllh)
    
  }  
  
  siginit <- sqrt(6 * var(c(s1)))/pi
  
  mu1_init <- mean(s1) - 0.57722 * siginit
  
  xi_init <- 0.1
  
  mu_t_init <- 0
  
  if(stationary){
    
    init_par <- c(mu1_init,siginit,xi_init) 
    
  }else{
    
    init_par <- c(mu1_init,siginit,xi_init,mu_t_init)
    
  }
  
  fit <- nlm(gev_ll,init_par,hessian = T)
  
  if(stationary){
    
    mu_t <- 0
    
  }else{
    
    mu_t <- fit$estimate[4]
    
  }
  
  return_mat <- matrix(0,length(snow_list),length(year_return))
  
  return_mat[1,] <- (fit$estimate[1] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[2]/fit$estimate[3])*(1-(-log(1-(1/year_return)))^(-fit$estimate[3])))
  
  return_mat[2,] <- (fit$estimate[1] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[2]/fit$estimate[3])*(1-(-log(1-(1/year_return)))^(-fit$estimate[3])))
  
  return_mat[3,] <- (fit$estimate[1] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[2]/fit$estimate[3])*(1-(-log(1-(1/year_return)))^(-fit$estimate[3])))
  
  
  return(list(returns=return_mat, trend=mu_t, output=fit))
  
}

stationary_mod3 <- model3_gev_nlme(snow_list=snow_list_set)

nonstationary_mod3 <- model3_gev_nlme(snow_list_set,time_list_set,stationary = F)


stationary_mod3_AICc <- AICc(stationary_mod3,snow_list_set)
nonstationary_mod3_AICc <- AICc(nonstationary_mod3,snow_list_set)

stationary_mod3_AICc
nonstationary_mod3_AICc
