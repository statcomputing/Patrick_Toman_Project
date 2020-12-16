########################################################################################
##### Name: smith_method_gev.R                                                              
##### Description: Smith's method adjust the standard errors of MLE fit to account for 
#####              spatial correlations between stations   
########################################################################################

########################################################################################
### Smith Correction for Reduced GEV Model 1
########################################################################################


model1_gev_smith <- function(model_fit,snow_list,stationary=T){

  model_fit <- model_fit$output
  
  H <- model_fit$hessian
  
  num_years <- length(snow_list[[1]])
  
  num_stations <- length(snow_list)
  
  mu_s1_mle <- model_fit$estimate[1]
  
  mu_s2_mle <- model_fit$estimate[2]
  
  mu_s3_mle <- model_fit$estimate[3]
  
  sigma_mle <- model_fit$estimate[4]
  
  xi_mle <- model_fit$estimate[5]
  
  trend_mle <- 0 
  
  if(!stationary){
    
    trend_mle <- model_fit$estimate[6]
    
  }
  
  se_original <- sqrt(diag(solve(H)))
  
  # Initialize vectors for new MLE/SE estimation via Smith's Method
  
  mu_s1_year <- numeric()
  
  mu_s2_year <- numeric()
  
  mu_s3_year <- numeric()
  
  sigma_year <- numeric()
  
  xi_year <- numeric()
  
  trend_year <- numeric()
  
  for(i in 1:num_years){
    
    mu_s1_mle_tmp <- mu_s1_mle
    
    mu_s2_mle_tmp <- mu_s2_mle
    
    mu_s3_mle_tmp <- mu_s3_mle
    
    sigma_mle_tmp <- sigma_mle
    
    xi_mle_tmp <- xi_mle
    
    diff_s1 <- snow_list[[1]][i]- mu_s1_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff_s2 <- snow_list[[2]][i]- mu_s2_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff_s3 <- snow_list[[3]][i]- mu_s3_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff <- c(diff_s1, diff_s2, diff_s3)
    
    diff <- diff[!is.na(diff)]
    
    diff_xi_s1_sc <- 1+xi_mle*diff_s1/sigma_mle
    
    diff_xi_s2_sc <- 1+xi_mle*diff_s2/sigma_mle
    
    diff_xi_s3_sc <- 1+xi_mle*diff_s3/sigma_mle
    
    diff_xi_sc <- c(diff_xi_s1_sc, diff_xi_s2_sc,diff_xi_s3_sc)
    
    diff_xi_sc <- diff_xi_sc[!is.na(diff_xi_sc)]
    
    # partial derivative of log likelihood function w.r.t. mu.t parameter
    if(!stationary){
    
      trend_year[i]<- -sum(diff_xi_sc^(-1-1/xi_mle))*c(1:num_years/10)[i]/sigma_mle +
        c(1:num_years/10)[i]*(1+xi_mle)*sum(1/(xi_mle*diff+sigma_mle))
    
      }else{
      
      trend_year[i]<-0
    
    }
    
    # partial derivatives of log likelihood function w.r.t. mu parameters
    mu_s1 <- -sum(diff_xi_s1_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s1+sigma_mle))
    
    mu_s2 <- -sum(diff_xi_s2_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s2+sigma_mle))
    
    mu_s3 <- -sum(diff_xi_s3_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s3+sigma_mle))
    
    # Combining partial derivatives as needed
    
    mu_s1_year[i] <- mu_s1
    
    mu_s2_year[i] <- mu_s2
    
    mu_s3_year[i] <- mu_s3
    
    # partial derivative of log likelihood function w.r.t. scale parameter
    sigma_year[i]<- -length(diff)/sigma_mle-sum(diff*(diff_xi_sc^(-1-1/xi_mle)))/sigma_mle^2+
      (1+xi_mle)*sum(diff/(xi_mle*sigma_mle*diff+sigma_mle^2))
    
    # partial derivative of log likelihood function w.r.t. shape parameter
    xi_year[i]<-sum(log(diff_xi_sc)/xi_mle^2-(1+1/xi_mle)*diff/(sigma_mle*diff_xi_sc))-
      sum(diff_xi_sc^(-1/xi_mle)*(log(diff_xi_sc)/xi_mle^2-diff/(xi_mle*sigma_mle*diff_xi_sc)))
    
  }
  
  # Setting up the gradient vector
  if(sum(trend_year)==0){
    
    new_par_matrix <- cbind(mu_s1_year,mu_s2_year,mu_s3_year,sigma_year,xi_year)
    
    colnames(new_par_matrix) <- c('mu_s1','mu_s2','mu_s3',
                                  'sigma','xi')
    
  
    }else{
      
      new_par_matrix <- cbind(mu_s1_year,mu_s2_year,mu_s3_year,sigma_year,xi_year,trend_year)
      
      colnames(new_par_matrix) <- c('mu_s1','mu_s2','mu_s3',
                                    'sigma','xi','trend')
    
  }
  
  # Calculating the V matrix
  
  if(any(is.na(new_par_matrix))){
    
    V <- cov(new_par_matrix,use='pairwise.complete.obs')*num_years
    
  }else{
   
     V <- cov(new_par_matrix)*num_years
     
  }
  
  
  # Calculating the 'corrected' SE
  new_hessian <- (solve(H) %*% V) %*% solve(H)
  
  se_new <- sqrt(diag(new_hessian))
  
  return(list(old=se_original, new=se_new))
  
  
}

########################################################################################
### Smith Correction for Reduced GEV Model 2 
########################################################################################


model2_gev_smith <- function(model_fit,snow_list,stationary=T){
  
  model_fit <- model_fit$output
  
  H <- model_fit$hessian
  
  num_years <- length(snow_list[[1]])
  
  num_stations <- length(snow_list)
  
  mu_s1_mle <- model_fit$estimate[1]
  
  mu_s2_mle <- model_fit$estimate[2]
  
  sigma_mle <- model_fit$estimate[3]
  
  xi_mle <- model_fit$estimate[4]
  
  trend_mle <- 0 
  
  if(!stationary){
    
    trend_mle <- model_fit$estimate[5]
    
  }
  
  se_original <- sqrt(diag(solve(H)))
  
  # Initialize vectors for new MLE/SE estimation via Smith's Method
  
  mu_s1_year <- numeric()
  
  mu_s2_year <- numeric()
  
  sigma_year <- numeric()
  
  xi_year <- numeric()
  
  trend_year <- numeric()
  
  for(i in 1:num_years){
    
    mu_s1_mle_tmp <- mu_s1_mle
    
    mu_s2_mle_tmp <- mu_s2_mle
    
    sigma_mle_tmp <- sigma_mle
    
    xi_mle_tmp <- xi_mle
    
    diff_s1 <- snow_list[[1]][i]- mu_s1_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff_s2 <- snow_list[[2]][i]- mu_s1_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff_s3 <- snow_list[[3]][i]- mu_s2_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff<-c(diff_s1, diff_s2,diff_s3)
    
    diff <- diff[!is.na(diff)]
    
    diff_xi_s1_sc <- 1+xi_mle*diff_s1/sigma_mle
    
    diff_xi_s2_sc <- 1+xi_mle*diff_s2/sigma_mle
    
    diff_xi_s3_sc <- 1+xi_mle*diff_s3/sigma_mle
    
    diff_xi_sc <- c(diff_xi_s1_sc, diff_xi_s2_sc,diff_xi_s3_sc)
    
    diff_xi_sc <- diff_xi_sc[!is.na(diff_xi_sc)]
    
    # partial derivative of log likelihood function w.r.t. mu.t parameter
    if(!stationary){
      
      trend_year[i] <- (-sum(diff_xi_sc^(-1-1/xi_mle))*c(1:num_years/10)[i]/sigma_mle +
        c(1:num_years/10)[i]*(1+xi_mle)*sum(1/(xi_mle*diff+sigma_mle)))
      
    }else{
      
      trend_year[i] <- 0
      
    }
    
    # partial derivatives of log likelihood function w.r.t. mu parameters
    mu_s1 <- -sum(diff_xi_s1_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s1+sigma_mle))
    
    mu_s2 <- -sum(diff_xi_s2_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s2+sigma_mle))
    
    mu_s3 <- -sum(diff_xi_s3_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s3+sigma_mle))
    
    # Combining partial derivatives as needed
    
    mu_s1_year[i] <- mu_s1 + mu_s2
    
    mu_s2_year[i] <- mu_s3
    
    # partial derivative of log likelihood function w.r.t. scale parameter
    sigma_year[i]<- -length(diff)/sigma_mle-sum(diff*(diff_xi_sc^(-1-1/xi_mle)))/sigma_mle^2+
      (1+xi_mle)*sum(diff/(xi_mle*sigma_mle*diff+sigma_mle^2))
    
    # partial derivative of log likelihood function w.r.t. shape parameter
    xi_year[i]<-sum(log(diff_xi_sc)/xi_mle^2-(1+1/xi_mle)*diff/(sigma_mle*diff_xi_sc))-
      sum(diff_xi_sc^(-1/xi_mle)*(log(diff_xi_sc)/xi_mle^2-diff/(xi_mle*sigma_mle*diff_xi_sc)))
    
  }
  
  # Setting up the gradient vector
  if(sum(trend_year)==0){
    
    new_par_matrix <- cbind(mu_s1_year,mu_s2_year,sigma_year,xi_year)
    
    colnames(new_par_matrix) <- c('mu_s1','mu_s2','sigma','xi')
    
    
  }else{
    
    new_par_matrix <- cbind(mu_s1_year,mu_s2_year,sigma_year,xi_year,trend_year)
    
    colnames(new_par_matrix) <- c('mu_s1','mu_s2','sigma','xi','trend')
    
  }
  
  # Calculating the V matrix
  
  if(any(is.na(new_par_matrix))){
    
    V <- cov(new_par_matrix,use='pairwise.complete.obs')*num_years
    
  }else{
    
    V <- cov(new_par_matrix)*num_years
    
  }
  
  
  # Calculating the 'corrected' SE
  new_hessian <- (solve(H) %*% V) %*% solve(H)
  
  se_new <- sqrt(diag(new_hessian))
  
  return(list(old=se_original, new=se_new))
  
  
}


########################################################################################
### Smith Correction for Reduced GEV Model 3 
########################################################################################

model3_gev_smith <- function(model_fit,snow_list,stationary=T){
  
  model_fit <- model_fit$output
  
  H <- model_fit$hessian
  
  num_years <- length(snow_list[[1]])
  
  num_stations <- length(snow_list)
  
  mu_s1_mle <- model_fit$estimate[1]
  
  sigma_mle <- model_fit$estimate[2]
  
  xi_mle <- model_fit$estimate[3]
  
  trend_mle <- 0 
  
  if(!stationary){
    
    trend_mle <- model_fit$estimate[4]
    
  }
  
  se_original <- sqrt(diag(solve(H)))
  
  # Initialize vectors for new MLE/SE estimation via Smith's Method
  
  mu_s1_year <- numeric()
  
  sigma_year <- numeric()
  
  xi_year <- numeric()
  
  trend_year <- numeric()
  
  for(i in 1:num_years){
    
    mu_s1_mle_tmp <- mu_s1_mle
    
    sigma_mle_tmp <- sigma_mle
    
    xi_mle_tmp <- xi_mle
    
    diff_s1 <- snow_list[[1]][i]- mu_s1_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff_s2 <- snow_list[[2]][i]- mu_s1_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff_s3 <- snow_list[[3]][i]- mu_s1_mle_tmp - trend_mle*c(1:num_years/10)[i]
    
    diff<-c(diff_s1, diff_s2,diff_s3)
    
    diff <- diff[!is.na(diff)]
    
    diff_xi_s1_sc <- 1+xi_mle*diff_s1/sigma_mle
    
    diff_xi_s2_sc <- 1+xi_mle*diff_s2/sigma_mle
    
    diff_xi_s3_sc <- 1+xi_mle*diff_s3/sigma_mle
    
    diff_xi_sc <- c(diff_xi_s1_sc, diff_xi_s2_sc,diff_xi_s3_sc)
    
    diff_xi_sc <- diff_xi_sc[!is.na(diff_xi_sc)]
    
    # partial derivative of log likelihood function w.r.t. mu.t parameter
    if(!stationary){
      
      trend_year[i] <- (-sum(diff_xi_sc^(-1-1/xi_mle))*c(1:num_years/10)[i]/sigma_mle +
                          c(1:num_years/10)[i]*(1+xi_mle)*sum(1/(xi_mle*diff+sigma_mle)))
      
    }else{
      
      trend_year[i] <- 0
      
    }
    
    # partial derivatives of log likelihood function w.r.t. mu parameters
    mu_s1 <- -sum(diff_xi_s1_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s1+sigma_mle))
    
    mu_s2 <- -sum(diff_xi_s2_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s2+sigma_mle))
    
    mu_s3 <- -sum(diff_xi_s3_sc^(-1-1/xi_mle))/sigma_mle+sum((1+xi_mle)/(xi_mle*diff_s3+sigma_mle))
    
    # Combining partial derivatives as needed
    
    mu_s1_year[i] <- mu_s1 + mu_s2 + mu_s3
    
    # partial derivative of log likelihood function w.r.t. scale parameter
    sigma_year[i] <- -length(diff)/sigma_mle-sum(diff*(diff_xi_sc^(-1-1/xi_mle)))/sigma_mle^2+
      (1+xi_mle)*sum(diff/(xi_mle*sigma_mle*diff+sigma_mle^2))
    
    # partial derivative of log likelihood function w.r.t. shape parameter
    xi_year[i] <- sum(log(diff_xi_sc)/xi_mle^2-(1+1/xi_mle)*diff/(sigma_mle*diff_xi_sc))-
      sum(diff_xi_sc^(-1/xi_mle)*(log(diff_xi_sc)/xi_mle^2-diff/(xi_mle*sigma_mle*diff_xi_sc)))
    
  }
  
  # Setting up the gradient vector
  if(sum(trend_year)==0){
    
    new_par_matrix <- cbind(mu_s1_year,sigma_year,xi_year)
    
    colnames(new_par_matrix) <- c('mu_s1','sigma','xi')
    
    
  }else{
    
    new_par_matrix <- cbind(mu_s1_year,sigma_year,xi_year,trend_year)
    
    colnames(new_par_matrix) <- c('mu_s1','sigma','xi','trend')
    
  }
  
  # Calculating the V matrix
  
  if(any(is.na(new_par_matrix))){
    
    V <- cov(new_par_matrix,use='pairwise.complete.obs')*num_years
    
  }else{
    
    V <- cov(new_par_matrix)*num_years
    
  }
  
  
  # Calculating the 'corrected' SE
  new_hessian <- (solve(H) %*% V) %*% solve(H)
  
  se_new <- sqrt(diag(new_hessian))
  
  return(list(old=se_original, new=se_new))

}

