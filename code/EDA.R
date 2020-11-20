########################################################################################
##### Name: EDA.R                                                                  #####
#####  Description: Exploratory Data Analysis for Chicago Snowfall Data            #####
########################################################################################
rm(list=ls())
library(tidyverse)
library(lubridate)

snowfall_df <- read_csv("data/SnowfallCHI_v5.csv")

snowfall_df %>% 
  group_by(station,name) %>% 
  summarise(min_date = min(date),max_date = max(date),
            num_obs = n())


snowfall_df %>% 
  summarise(ct = n(),snowf_missing = sum(is.na(snowf)),snowf_pct_missing = sum(is.na(snowf))/n(),
            snowdep_missing = sum(is.na(snowdep)),snowdep_pct_missing = sum(is.na(snowdep))/n()) -> missingness_df 

snowfall_plt <- snowfall_df %>% 
  ggplot(aes(x=date,y=snowf,color=name)) +
  geom_line() + theme_minimal() + theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_x_date(date_breaks = '5 years') +facet_wrap(~name,nrow = 3,scales = 'free_y')

snowfall_plt


snowfall_df %>% 
  select(-snowdep) %>% 
  filter(snowf > 0) %>% 
  nest(data = c(date,snowf)) -> snowdays_df 


