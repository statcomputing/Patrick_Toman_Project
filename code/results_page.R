########################################################################################
##### Name: run_models.R                                                              
##### Description: Generates all output from report including plots 
########################################################################################

rm(list=ls())
library(tidyverse)
library(lubridate)
source("./code/gev_mle_functions.R") 
source("./code/smith_method_gev.R")
source("./code/gev_bootstrap.R")

# Read in all data and set up a list of dataframes

annual_max_snowf_df <- read_csv("data/final_snowfall_data.csv")

ann_max_snowf_split <- split(annual_max_snowf_df,f=annual_max_snowf_df$name)

# Number of years

n.years <- nrow(ann_max_snowf_split$Midway)

# Return level periods of interest

return_levels <- c(25,50,75,100)

# make trend covariate

year_set <- 1:n.years/10

# Set up a list of containing vectors for the annual maxes of each dataframe 

snowf_list <- list(ann_max_snowf_split$Midway$ann_max_snowf,
                      ann_max_snowf_split$Ohare$ann_max_snowf,
                      ann_max_snowf_split$ParkForest$ann_max_snowf)


names(snowf_list) <- names(ann_max_snowf_split)

# Set up a list containing a vector of trend covariates one -  for each station

trend_list <- list(year_set,year_set,year_set)

names(trend_list) <- names(ann_max_snowf_split)

#####################################################################################
### SECTION 4 - MODEL FITTING ###
#####################################################################################

#########################################################
# Full Model (section 4.1)
#########################################################

full_model <- full_model_gev_nlme(snowf_list,trend_list)

#########################################################
# stationary reduced models (section 4.2)
#########################################################

# Model Fits

reduced_model1_stationary <- model1_gev_nlme(snowf_list,trend_list)

reduced_model2_stationary <- model2_gev_nlme(snowf_list,trend_list)

reduced_model3_stationary <- model3_gev_nlme(snowf_list,trend_list)

# Smith's Correction for Spatial Dependence

reduced_model1_stationary_sc <- model1_gev_smith(reduced_model1_stationary,snowf_list)

reduced_model2_stationary_sc <- model2_gev_smith(reduced_model2_stationary,snowf_list)

reduced_model3_stationary_sc <- model3_gev_smith(reduced_model3_stationary,snowf_list)


#########################################################
# non-stationary reduced models (section 4.2)
#########################################################


# Model Fits
stationary_set <- F

reduced_model1_non_stationary <- model1_gev_nlme(snowf_list,trend_list,
                                                 stationary = stationary_set)

reduced_model2_non_stationary <- model2_gev_nlme(snowf_list,trend_list,
                                                 stationary = stationary_set)

reduced_model3_non_stationary <- model3_gev_nlme(snowf_list,trend_list,
                                                 stationary = stationary_set)

# Smith's Correction for Spatial Dependence

reduced_model1_non_stationary_sc <- model1_gev_smith(reduced_model1_non_stationary,
                                                     snowf_list,stationary = stationary_set)

reduced_model2_non_stationary_sc <- model2_gev_smith(reduced_model2_non_stationary,
                                                     snowf_list,stationary = stationary_set)

reduced_model3_non_stationary_sc <- model3_gev_smith(reduced_model3_non_stationary,
                                                     snowf_list,stationary = stationary_set)

###################################################################################################
### SECTION 5 - BOOTSTRAP RETURN LEVELS ###
###################################################################################################

##################################
### stationary models ###
##################################

# Model 1

zBC_stationary_mod1 <- zBCcalc(snowf_list,year_return = return_levels,model = 1,
                               original_model = reduced_model1_stationary)

cA_stationary_mod1 <- cAcalc(snowf_list)

stationary_mod1_ci <- ReturnLvlBootCI(samples = zBC_stationary_mod1$Sample,zBC = zBC_stationary_mod1$zBC,cA = cA_stationary_mod1)

# Model 2

zBC_stationary_mod2 <- zBCcalc(snowf_list,year_return = return_levels,model = 2,
                               original_model = reduced_model2_stationary)

cA_stationary_mod2 <- cAcalc(snowf_list,model = 2)

stationary_mod2_ci <- ReturnLvlBootCI(samples = zBC_stationary_mod2$Sample,
                                      zBC = zBC_stationary_mod2$zBC,cA = cA_stationary_mod2)

# Model 3

zBC_stationary_mod3 <- zBCcalc(snowf_list,year_return = return_levels,model = 3,
                               original_model = reduced_model3_stationary)

cA_stationary_mod3 <- cAcalc(snowf_list,model = 3)

stationary_mod3_ci <- ReturnLvlBootCI(samples = zBC_stationary_mod3$Sample,
                                      zBC = zBC_stationary_mod3$zBC,cA = cA_stationary_mod3)

##################################
### non-stationary models ###
##################################

# Model 1

zBC_nonstationary_mod1 <- zBCcalc(snowf_list,year_return = return_levels,model = 1,
                               original_model = reduced_model1_non_stationary,stationary = F)

cA_nonstationary_mod1 <- cAcalc(snowf_list)

nonstationary_mod1_ci <- ReturnLvlBootCI(samples = zBC_nonstationary_mod1$Sample,
                                         zBC = zBC_nonstationary_mod1$zBC,cA = cA_nonstationary_mod1)

# Model 2

zBC_nonstationary_mod2 <- zBCcalc(snowf_list,year_return = return_levels,model = 2,
                               original_model = reduced_model2_non_stationary,stationary = F)

cA_nonstationary_mod2 <- cAcalc(snowf_list,model = 2)

nonstationary_mod2_ci <- ReturnLvlBootCI(samples = zBC_nonstationary_mod2$Sample,
                                         zBC = zBC_nonstationary_mod2$zBC,cA = cA_stationary_mod2)

# Model 3

zBC_nonstationary_mod3 <- zBCcalc(snowf_list,year_return = return_levels,model = 3,
                               original_model = reduced_model3_non_stationary,stationary = F)

cA_nonstationary_mod3 <- cAcalc(snowf_list,model = 3)

nonstationary_mod3_ci <- ReturnLvlBootCI(samples = zBC_nonstationary_mod3$Sample,
                                         zBC = zBC_nonstationary_mod3$zBC,cA = cA_nonstationary_mod3)

###################################################################################################
### FINAL RESULTS ###
###################################################################################################

# Function to make table of parameter estimates and standard errors for a GEV model (except for full model)

build_results_table <- function(gev_model,gev_model_se){
  
 result <- rbind(gev_model$output$estimate,
                 gev_model_se$old,
                 gev_model_se$new)
 
  rownames(result) <- c('Estimate','se_naive','se_smith')
 
  return(result) 
}


##################################
### parameter estimates & AICc ###
##################################

### full model 

full_model_params <- rbind(full_model$output$estimate,sqrt(diag(solve(full_model$output$hessian))))

full_model_aicc <- AICc(full_model,snowf_list)

full_model_ll <- full_model$output$minimum

### stationary models 

# model 1

stationary_mod1_params <- build_results_table(reduced_model1_stationary,reduced_model1_stationary_sc)

stationary_mod1_aicc <- AICc(reduced_model1_stationary,snowf_list)

stationary_mod1_ll <- reduced_model1_stationary$output$minimum

# model 2

stationary_mod2_params <- build_results_table(reduced_model2_stationary,reduced_model2_stationary_sc)

stationary_mod2_aicc <- AICc(reduced_model2_stationary,snowf_list)

stationary_mod2_ll <- reduced_model2_stationary$output$minimum

# model 3

stationary_mod3_params <- build_results_table(reduced_model3_stationary,reduced_model3_stationary_sc)

stationary_mod3_aicc <- AICc(reduced_model3_stationary,snowf_list)

stationary_mod3_ll <- reduced_model3_stationary$output$minimum

### Non-stationary models

# model 1

nonstationary_mod1_params <- build_results_table(reduced_model1_non_stationary,reduced_model1_non_stationary_sc)

nonstationary_mod1_aicc <- AICc(reduced_model1_non_stationary,snowf_list)

nonstationary_mod1_ll <- reduced_model1_non_stationary$output$minimum

# model 2

nonstationary_mod2_params <- build_results_table(reduced_model2_non_stationary,reduced_model2_non_stationary_sc)

nonstationary_mod2_aicc <- AICc(reduced_model2_non_stationary,snowf_list)

nonstationary_mod2_ll <- reduced_model2_non_stationary$output$minimum

# model 3

nonstationary_mod3_params <- build_results_table(reduced_model3_non_stationary,reduced_model3_non_stationary_sc)

nonstationary_mod3_aicc <- AICc(reduced_model3_non_stationary,snowf_list)

nonstationary_mod3_ll <- reduced_model3_non_stationary$output$minimum

##################################
### Bootstrap BCa Intervals ###
##################################
# Stationary Models

k25_stationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = stationary_mod1_ci$midway[1,],
                                        "O'hare" = stationary_mod1_ci$ohare[1,],
                                        "Park" = stationary_mod1_ci$parkforest[1,])

k50_stationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = stationary_mod1_ci$midway[2,],
                                        "O'hare" = stationary_mod1_ci$ohare[2,],
                                        "Park" = stationary_mod1_ci$parkforest[2,])

k75_stationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = stationary_mod1_ci$midway[3,],
                                        "O'hare" = stationary_mod1_ci$ohare[3,],
                                        "Park" = stationary_mod1_ci$parkforest[3,])

k100_stationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                         "Midway" = stationary_mod1_ci$midway[4,],
                                         "O'hare" = stationary_mod1_ci$ohare[4,],
                                         "Park" = stationary_mod1_ci$parkforest[4,])

stationary_plt_list <- list("K25" = k25_stationary_mod1,"K50" = k50_stationary_mod1,
                            "K75" = k75_stationary_mod1,"K100" = k100_stationary_mod1)

stationary_plt_df <- bind_rows(stationary_plt_list,.id = "K") 

stationary_plt_df %>% 
  pivot_longer(cols = c(3:5)) -> stationary_plt_df

stationary_plt_df$K <- factor(as.numeric(gsub("K","",stationary_plt_df$K)),labels = paste0("K",c(25,50,75,100)))


colnames(stationary_plt_df) <- c("K","%","Station","Max Snowfall")

stationary_plt_df %>% 
  group_by(K,Station) %>% 
  ggplot(aes(x = Station,y = `Max Snowfall`,group = Station,color=Station)) + 
  geom_boxplot() + theme_minimal() +
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 45)) + 
  facet_wrap(~K, nrow = 1, ncol = 4,scales = 'free_x') + 
  ggtitle("95% Bootstrap CI for Return Levels \n (Stationary)")

# Non-stationary Models

k25_nonstationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = nonstationary_mod1_ci$midway[1,],
                                        "O'hare" = nonstationary_mod1_ci$ohare[1,],
                                        "Park" = nonstationary_mod1_ci$parkforest[1,])

k50_nonstationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = nonstationary_mod1_ci$midway[2,],
                                        "O'hare" = nonstationary_mod1_ci$ohare[2,],
                                        "Park" = nonstationary_mod1_ci$parkforest[2,])

k75_nonstationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = nonstationary_mod1_ci$midway[3,],
                                        "O'hare" = nonstationary_mod1_ci$ohare[3,],
                                        "Park" = nonstationary_mod1_ci$parkforest[3,])

k100_nonstationary_mod1 <- cbind.data.frame("%" = c("2.5%","Median","97.5%"),
                                        "Midway" = nonstationary_mod1_ci$midway[4,],
                                        "O'hare" = nonstationary_mod1_ci$ohare[4,],
                                        "Park" = nonstationary_mod1_ci$parkforest[4,])

nonstationary_plt_list <- list("K25" = k25_nonstationary_mod1,"K50" = k50_nonstationary_mod1,
                            "K75" = k75_nonstationary_mod1,"K100" = k100_nonstationary_mod1)

nonstationary_plt_df <- bind_rows(nonstationary_plt_list,.id = "K") 

nonstationary_plt_df %>% 
  pivot_longer(cols = c(3:5)) -> nonstationary_plt_df

nonstationary_plt_df$K <- factor(as.numeric(gsub("K","",nonstationary_plt_df$K)),labels = paste0("K",c(25,50,75,100)))

colnames(nonstationary_plt_df) <- c("K","%","Station","Max Snowfall")

nonstationary_plt_df %>% 
  group_by(K,Station) %>% 
  ggplot(aes(x = Station,y = `Max Snowfall`,group = Station,color=Station)) + 
  geom_boxplot() + theme_minimal() +
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 45)) + 
  facet_wrap(~K, nrow = 1, ncol = 4,scales = 'free_x') + 
  ggtitle("95% Bootstrap CI for Return Levels \n (Non-stationary)")


