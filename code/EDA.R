########################################################################################
##### Name: EDA.R                                                                  #####
#####  Description: Exploratory Data Analysis for Chicago Snowfall Data            #####
########################################################################################
rm(list=ls())
library(tidyverse)
library(lubridate)

snowfall_df <- read_csv("data/SnowfallCHI_v5.csv")

# Summary Statistics

snowfall_df %>% 
  group_by(station,name) %>% 
  summarise(min_date = min(date),max_date = max(date),
            num_obs = n())


snowfall_df %>% 
  group_by(station,name,month(date)) %>% 
  summarise(ct = n(),snowf_missing = sum(is.na(snowf)),snowf_pct_missing = sum(is.na(snowf))/n(),
            snowdep_missing = sum(is.na(snowdep)),snowdep_pct_missing = sum(is.na(snowdep))/n()) -> missingness_df 
snowfall_df %>% 
  group_by(station,name,'year'=year(date)) %>% 
  summarise(ann_max_snowf = max(snowf,na.rm = T)) -> annual_max_snowf_df 

snowfall_df$missing_ts <- 0

snowfall_df$missing_ts[which(is.na(snowfall_df$snowf)==T)] <- 1



# Plots
snowfall_plt <- snowfall_df %>% 
  ggplot(aes(x=date,y=snowf,color=name)) +
  geom_line() + theme_minimal() + theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_x_date(date_breaks = '5 years') +xlab('Date')+ylab('Snowfall (inches)') + 
  ggtitle('Daily Snowfall')+ facet_wrap(~name,nrow = 3,scales = 'free_y')

missing_plt <- snowfall_df %>% 
  ggplot(aes(x=date,y=missing_ts,color=name)) +
  geom_line() + theme_minimal() + theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_x_date(date_breaks = '5 years') +xlab('Date')+ylab('Missing Indicator') + 
  ggtitle('Missingness Plots')+ facet_wrap(~name,nrow = 3,scales = 'free_y')


ann_max_snowfall_plt <- annual_max_snowf_df %>% 
  ggplot(aes(x=as.numeric(year),y=ann_max_snowf,color=name)) +
  geom_line() + theme_minimal() + theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab('Date')+ylab('Maximum Snowfall (inches)')+ ggtitle('Annual Maxes')+
  facet_wrap(~name,nrow = 3,scales = 'free_y')

ann_max_snowfall_plt
missing_plt
snowfall_plt


snowfall_df %>% 
  select(-snowdep) %>% 
  filter(snowf > 0) %>% 
  nest(data = c(date,snowf)) -> snowdays_df 

# Splits into groups of consecutive and non-consecutive days

date_split <- function(df){
  
  split(df$date, cumsum(c(TRUE, diff(df$date) != 1)))
  
}

snowdays_df %>% 
  mutate('date_splits'=lapply(snowdays_df$data,date_split)) -> snowdays_df

