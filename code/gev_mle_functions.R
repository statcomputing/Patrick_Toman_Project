########################################################################################
##### Name: gev_mle_functions.R                                                              
##### Description: GEV MLE functions (non-stationary and stationary) plus AICc functions
####  Purpose - Each function returs parameter estimates and return levels
########################################################################################

########################################################################################
### AICc for model comparison
########################################################################################


AICc <- function(gev_mod,snow_data){
  
  ll <- gev_mod$output$minimum
  
  n <- length(unlist(snow_data))
  
  p <- length(gev_mod$output$estimate)
  
  aicc <- 2*(p - ll + (p*(p+1)) / (n-p-1))  
  
  return(aicc)
  
}


########################################################################################
### Full Model described in section 4.1 - All parameters and return levels are estimated
###                                       at each station
########################################################################################


full_model_gev_nlme <- function(snow_list,time_list,n.years = 61,year_return = c(25,50,75,100),stationary = T){
  
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
  
  colnames(return_mat) <- as.character(year_return)
  
  rownames(return_mat) <- names(snow_list)
  
  
  return(list(returns=return_mat, trend=mu_t, output=fit))
  
}

########################################################################################
### Reduced Model 1 Section 4.2 - All parameters and return levels are estimated
###                               at each station
########################################################################################


model1_gev_nlme <- function(snow_list,time_list,year_return = c(25,50,75,100),n.years = 61,stationary = T){
  
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
  
  colnames(return_mat) <- as.character(year_return)
  
  rownames(return_mat) <- names(snow_list)
  
  return(list(returns=return_mat, trend=mu_t, output=fit))
  
}


########################################################################################
### Reduced Model 2 Section 4.2 - All parameters and return levels are estimated
###                               at each station
########################################################################################

model2_gev_nlme <- function(snow_list,time_list,n.years = 61,year_return = c(25,50,75,100)
                            ,stationary = T,groups = list('g1'=c(1),'g2'=c(2,3))){
  
  num_stations <- length(snow_list)
  
  g1_name <- paste0(names(snow_list)[groups$g1[1]])
  
  g2_name <- names(snow_list)[groups$g2]
  
  
  if(!stationary){
    
    s1 <- unlist(snow_list[groups$g1])
    
    s2 <- unlist(snow_list[groups$g2])
    
    t1 <- c(unlist(time_list[groups$g1[1]]),unlist(time_list[groups$g1[2]]))
    
    t2 <- unlist(time_list[groups$g2])
    
    
  }else{
    
    
    s1 <- unlist(snow_list[groups$g1])
    
    s2 <- unlist(snow_list[groups$g2])
    
  }
  
  
  num_obs <- length(c(s1,s2))
  
  gev_ll <- function(theta){
    
    mu1 <- theta[1]
    
    mu2 <- theta[2]
    
    sc <- theta[3]
    
    xi <- theta[4]
    
    if(!stationary){
      
      mu_t <- theta[5]
      
    }
    
    if(stationary){
      
      y1 <- 1 + xi * (s1 - mu1)/sc
      
      y2 <- 1 + xi * (s2 - mu2)/sc
      
    }else{
      
      y1 <- 1 + xi * (s1 - (mu1 + mu_t*t1))/sc
      
      y2 <- 1 + xi * (s2 - (mu2 + mu_t*t2))/sc
      
    }
    
    
    if(any(y1) <= 0 || any(y2 <= 0) ||  any(sc <= 0)){
      
      nllh <- 1e6
      
      
    }else{
      
      nllh <- num_obs*log(sc) + sum(y1^(-1/xi)) + sum(y2^(-1/xi))  +
        sum(log(y1) * (1/xi + 1)) + sum(log(y2) * (1/xi + 1)) 
      
    }
    
    return(nllh)
    
  }  
  
  siginit <- sqrt(6 * var(c(s1,s2)))/pi
  
  mu1_init <- mean(s1) - 0.57722 * siginit
  
  mu2_init <- mean(s2) - 0.57722 * siginit
  
  xi_init <- 0.1
  
  mu_t_init <- 0
  
  if(stationary){
    
    init_par <- c(mu1_init,mu2_init,siginit,xi_init) 
    
  }else{
    
    init_par <- c(mu1_init,mu2_init,siginit,xi_init,mu_t_init)
    
  }
  
  fit <- nlm(gev_ll,init_par,hessian = T)
  
  if(stationary){
    
    mu_t <- 0
    
  }else{
    
    mu_t <- fit$estimate[5]
    
  }
  
  return_mat <- matrix(0,length(snow_list),length(year_return))
  
  return_mat[1,] <- (fit$estimate[1] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[3]/fit$estimate[4])*(1-(-log(1-(1/year_return)))^(-fit$estimate[4])))
  
  
  return_mat[2,] <- (fit$estimate[2] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[3]/fit$estimate[4])*(1-(-log(1-(1/year_return)))^(-fit$estimate[4])))
  
  
  return_mat[3,] <- (fit$estimate[2] + mu_t*(year_return + n.years)/10 - 
                       (fit$estimate[3]/fit$estimate[4])*(1-(-log(1-(1/year_return)))^(-fit$estimate[4])))
  
  colnames(return_mat) <- as.character(year_return)
  
  rownames(return_mat) <- names(snow_list)
  
  
  return(list(returns=return_mat, trend=mu_t, output=fit,'g1'=g1_name,'g2'=g2_name))
  
}

########################################################################################
### Reduced Model 3 Section 4.2 - All parameters and return levels are estimated
###                               at each station
########################################################################################

model3_gev_nlme <- function(snow_list,time_list,n.years = 61,year_return = c(25,50,75,100),stationary = T){
  
  num_stations <- length(snow_list)
  
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

