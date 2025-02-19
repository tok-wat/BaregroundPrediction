library(tidyverse)
library(caret)
library(randomForest)
library(ggplot2)

############ Load CSV file ##############
df_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 1,006 rows

head(df_ground_photo)

# Load a script
source("./Script/functions.R")

bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
precipitation <- c("X0_2013_precipitation", "X10_2023_precipitation", "X11_2024_precipitation", "X1_2014_precipitation",  
                   "X2_2015_precipitation", "X3_2016_precipitation", "X4_2017_precipitation", "X5_2018_precipitation",  
                   "X6_2019_precipitation", "X7_2020_precipitation", "X8_2021_precipitation", "X9_2022_precipitation" )
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
topography <- c("elevation", "slope", "aspect")

predictors_1 <- c(vi)
predictors_2 <- c(vi, bands)
predictors_3 <- c(vi, bands, minerals)
predictors_4 <- c(vi, bands, minerals, topography)

############ RandomForest Analysis ##############
# random seed
set.seed(123)

df_non_growing <- df_ground_photo %>% 
  filter(Month >= 5)

# split in to 70% training and 30% test data
df_non_growing_train <- df_non_growing %>% 
  sample_frac(0.7)
df_non_growing_test <- df_non_growing %>% 
  anti_join(df_non_growing_train, by = "ID") 

####### Vegetation Indicies Only #################
# Recursive Feature Elimination (RFE) for variable selection
rfe_1 <- rfe_tune(df_non_growing_train, predictors_1) 
# running RandomForest regression with the selected variables with hyper-parameter tuning
rf_1 <- run_rf(df_non_growing_train, predictors_1, rfe_1)
# predictions with test data
rf_prediction_1 <- predict(rf_1, df_non_growing_test[,predictors_1])

## some visulisations of the process
# rfe_1$results %>% 
#   ggplot(aes(x=Variables, y=RMSE))+
#   geom_point()+
#   geom_errorbar(aes(ymin=RMSE-RMSESD, ymax=RMSE+RMSESD))+
#   theme_bw()
# 
# varImpPlot(rf_1)

# plot_obs_pred(df_non_growing_test$BareGround, rf_prediction_1)
# get_mae_rmse(rf_prediction_1, df_non_growing_test$BareGround)

####### Vegetation Indicies + bands #################
rfe_2 <- rfe_tune(df_non_growing_train, predictors_2) 
rf_2 <- run_rf(df_non_growing_train, predictors_2, rfe_2)
rf_prediction_2 <- predict(rf_2, df_non_growing_test[,predictors_2])

####### Vegetation Indicies + bands + minerals #################
rfe_3 <- rfe_tune(df_non_growing_train, predictors_3) 
rf_3 <- run_rf(df_non_growing_train, predictors_3, rfe_3)
rf_prediction_3 <- predict(rf_3, df_non_growing_test[,predictors_3])

####### Vegetation Indicies + bands + minerals + topography #################
rfe_4 <- rfe_tune(df_non_growing_train, predictors_4) 
rf_4 <- run_rf(df_non_growing_train, predictors_4, rfe_4)
rf_prediction_4 <- predict(rf_4, df_non_growing_test[,predictors_4])


######### comparison #########
# obs vs pred plot
rf_pred_list <- list(rf_prediction_1,rf_prediction_2,rf_prediction_3,rf_prediction_4) 

par(mfrow = c(2,2))
for(pred in rf_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}

# a table for RMSE and MAE
df_rf_rmse <- get_rmse_df(rf_pred_list)

######### Export ############
models <- list(rf_1, rf_2, rf_3, rf_4)
for(i in 1:length(models)){
  temp_filename <- paste0("rf_",i)
  saveRDS(models[[i]], file = paste0("./Script/models/", temp_filename, ".rds"))
}

# rmse df
write.csv(df_rf_rmse, "./Script/csv/rf_rmse.csv")

# plot
png("./Script/figures/rf.png", width = 600, height = 700)
par(mfrow = c(2,2))
for(pred in rf_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}
dev.off()  

