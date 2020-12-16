########################################################################################
##### Name: run_models.R                                                              
##### Description: Generates all output from report except for plots 
########################################################################################
rm(list=ls())
library(tidyverse)
library(lubridate)
source("./code/gev_mle_functions.R")
source("./code/smith_method_gev.R")

# Read in all data and set up a list of dataframes

annual_max_snowf_df <- read_csv("data/final_snowfall_data.csv")

ann_max_snowf_split <- split(annual_max_snowf_df,f=annual_max_snowf_df$name)

# Number of years

n.years <- nrow(ann_max_snowf_split$Midway)

# Return level periods of interest

year_return <- c(25,50,75,100)

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

reduced_model1_non_stationary_sc <- model1_gev_smith(reduced_model1_stationary,
                                                     snowf_list,stationary = stationary_set)

reduced_model2_non_stationary_sc <- model2_gev_smith(reduced_model2_stationary,
                                                     snowf_list,stationary = stationary_set)

reduced_model3_non_stationary_sc <- model3_gev_smith(reduced_model3_stationary,
                                                     snowf_list,stationary = stationary_set)

#####################################################################################
### SECTION 5 - BOOTSTRAP RETURN LEVELS ###
#####################################################################################

