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
#missing_plt
#snowfall_plt

snowfall_df %>% 
  select(-snowdep) %>% 
  filter(snowf > 0) %>% 
  nest(data = c(date,snowf)) -> snowdays_df 

# Splits into groups of consecutive and non-consecutive days

date_split <- function(df){
  
  split(df$date, cumsum(c(TRUE, diff(df$date) != 1)))
  
}



snow_events_list <- list()
date_splits_list <- list()

station_names <- snowdays_df$name

for(i in 1:length(station_names)){
  
  df_tmp <- snowdays_df %>% 
    filter(name == station_names[i]) 
  
  events_tmp <- date_split(df_tmp$data[[1]])
  
  df_tmp %>% 
    unnest(cols=c(data)) -> df_temp 
  
  snowfall_totals_tmp <- list()
  
  for(j in 1:length(events_tmp)){
    
    df_temp %>% 
      filter(date %in% events_tmp[[j]]) %>% 
      summarise(snowfall_total = sum(snowf)) -> snowf_tot_tmp
    
    snowfall_totals_tmp[[j]] <- snowf_tot_tmp$snowfall_total
    
    names(snowfall_totals_tmp)[j] <- as.character(events_tmp[[j]][1])
    
  }
  
  snow_events_list[[i]] <- snowfall_totals_tmp
  
  date_splits_list[[i]] <- events_tmp
  
  names(snow_events_list)[i] <- names(date_splits_list)[i] <- station_names[i]
  
}

########################################################################################
# Required Data based on authors data file layout 

# Annual maxes 
annmax_midway <- annual_max_snowf_df %>% 
  filter(name=='Midway') %>% 
  select(ann_max_snowf)

annmax_ohare <- annual_max_snowf_df %>% 
  filter(name=='Ohare') %>% 
  select(ann_max_snowf)

annmax_parkforest <- annual_max_snowf_df %>% 
  filter(name=='ParkForest') %>% 
  select(ann_max_snowf)

# Indices of days where snowfall events occur 

days_midway <- snowfall_df %>%
  filter(name=='Midway') 

days_midway <- which(days_midway$date %in% ymd(names(snow_events_list$Midway)))

days_ohare <- snowfall_df %>%
  filter(name=='Ohare') 

days_ohare <- which(days_ohare$date %in% ymd(names(snow_events_list$Ohare)))

days_parkforest <- snowfall_df %>%
  filter(name=='ParkForest') 

days_parkforest <- which(days_parkforest$date %in% ymd(names(snow_events_list$ParkForest)))

# Durations of snowevents 

duration_midway <- unlist(lapply(date_splits_list$Midway, length))

duration_ohare <- unlist(lapply(date_splits_list$Ohare, length))

duration_parkforest <- unlist(lapply(date_splits_list$ParkForest, length))

# Snow event values 
nonzero_snow_midway <- unlist(snow_events_list$Midway)

nonzero_snow_ohare <- unlist(snow_events_list$Ohare)

nonzero_snow_parkforest <- unlist(snow_events_list$ParkForest)

# Seasons 

# Function below gives season variables as authors defined variable in their code

make_season <- function(X){

  season_dates <- as.Date(names(unlist(X)))
  
  date_seq <- seq.Date(from = as.Date('1960-07-01'),to =as.Date('2020-06-30'),by='year')
    
  nseasons <- length(date_seq)
  
  season_tmp <- c()
  
  for(i in 1:(nseasons-1)){
    
    date_interval_tmp <- c(date_seq[i],date_seq[i+1])
    
    season_inds <- which((season_dates >= date_interval_tmp[1]) & (season_dates < date_interval_tmp[2]))  
    
    year_tmp <- rep(as.numeric(year(date_interval_tmp[1])),length(season_inds))
    
    season_tmp <- append(season_tmp,year_tmp)
    
  }
  
  return(season_tmp)
  
}

season_midway <- make_season(snow_events_list$Midway)

season_ohare <- make_season(snow_events_list$Ohare)

season_parkforest <- make_season(snow_events_list$ParkForest)

