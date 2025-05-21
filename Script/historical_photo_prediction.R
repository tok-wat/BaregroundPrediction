library(pls)
library(caret)
library(tidyverse)
library(ggplot2)

# ------------------------------------------------------------------------------
# load necessary files

# Load CSV file 
df_historical <- read.csv("./CSV/historical_photo_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  

# head(df_historical)

# Load a script# LoYearad a script
source("./Script/functions.R")

# Load the model
best_non_growing <- readRDS("./Script/models/pls_non_growing.rds")
best_growing <- readRDS("./Script/models/pls_growing.rds")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio","tc_bright","tc_green","tc_wet")

# ------------------------------------------------------------------------------
# non-growing season
df_historical_non_growing <- df_historical %>% 
  filter(Month > 4) # 29 rows

hist_non_growing <- predict(best_non_growing, df_historical_non_growing)
plot_obs_pred(df_historical_non_growing$BareGround, hist_non_growing)

# residuals or prediction error
data.frame(hist_non_growing, df_historical_non_growing$BareGround) %>% 
  mutate(pred_error = hist_non_growing-df_historical_non_growing.BareGround) %>% 
  mutate(Year = df_historical_non_growing$Year) %>%
  mutate(Month = df_historical_non_growing$Month) %>% 
  select(pred_error, Year, Month) %>% 
  plot()

# ------------------------------------------------------------------------------
# growing season
df_historical_growing <- df_historical %>% 
  filter(Month < 5) # 19 rows

# # log transformation
# trans_vars <- func_log_transform(sf_growing, vi)[[2]] # "bsi","evi","msavi", "ndvi", "osavi", "satvi","savi","swirRatio","tc_bright","tc_green" 
# add_log <- func_log_transform(sf_growing, vi)[[4]] # 0.367
# neg_var_names <- func_log_transform(sf_growing, vi)[[5]] # "bsi", "exg" , "ndmi" , "satvi" , "tc_wet"
# 
# df_historical_growing_log <- df_historical_growing %>% 
#   mutate(across(all_of(neg_var_names), ~ .x + add_log)) %>%
#   mutate(across(all_of(trans_vars), log))

hist_growing <- predict(best_growing, df_historical_growing)
plot_obs_pred(df_historical_growing$BareGround, hist_growing)

# residuals or prediction error
data.frame(hist_growing, df_historical_growing$BareGround) %>% 
  mutate(pred_error = hist_growing-df_historical_growing.BareGround) %>% 
  mutate(Year = df_historical_growing$Year) %>%
  mutate(Month = df_historical_growing$Month) %>% 
  select(pred_error, Year, Month) %>% 
  plot()

#-------------------------------------------------------------------------------
# Export Graphs

# export as png files
png("./Script/figures/historical_non_growing.png", width = 1500, height = 1600, res=300)  
plot_obs_pred(df_historical_non_growing$BareGround, hist_non_growing)
dev.off() 

png("./Script/figures/historical_growing.png", width = 1500, height = 1600, res=300)  
plot_obs_pred(df_historical_growing$BareGround, hist_growing)
dev.off() 

png("./Script/figures/historical_photos_non_growing.png", width = 1500, height = 1000, res=300)  
hist(df_historical_non_growing$Year,breaks=1989:2013,
     xlab="Year",
     ylab="Numbers of Photographs",
     xlim=c(1989,2015),
     ylim=c(0,15),
     main=NULL)
dev.off()

png("./Script/figures/historical_photos_growing.png", width = 1500, height = 1000, res=300)  
hist(df_historical_growing$Year,breaks=1989:2013,
     xlab="Year",
     ylab="Numbers of Photographs",
     xlim=c(1989,2015),
     ylim=c(0,15),
     main=NULL)
dev.off()
