library(tidyverse)
library(pls)
library(plsdepot)
library(plsVarSel)
library(plsRglm)

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

############ Pls Analysis ##############
# random seed
set.seed(123)

df_non_growing <- df_ground_photo %>% 
  filter(Month >= 5)

# split in to 70% training and 30% test data
df_non_growing_train <- df_non_growing %>% 
  sample_frac(0.7)
df_non_growing_test <- df_non_growing %>% 
  anti_join(df_non_growing_train, by = "ID") 

######### Run the PLS regression models ####

######### vegetation indicies only #########
pls_1 <- df_non_growing_train %>% 
  pls_variable_selection(predictors_1) %>% 
  pls_best_model(df_non_growing_train, c(predictors_1))
# get_q2(df_non_growing_train, pls_1)

# validation against test data
ncomp <- selectNcomp(pls_1, method = "onesigma")
pls_prediction_1 <- predict(pls_1, df_non_growing_test[,predictors_1], ncomp=ncomp)[,,1]

######### vegetation indicies + bands #########
pls_2 <- df_non_growing_train %>% 
  pls_variable_selection(predictors_2) %>% 
  pls_best_model(df_non_growing_train, c(predictors_2))
# get_q2(df_non_growing_train, pls_2)

# validation against test data
ncomp <- selectNcomp(pls_2, method = "onesigma")
pls_prediction_2 <- predict(pls_2, df_non_growing_test[,predictors_2], ncomp=ncomp)[,,1]

######### vegetation indicies + bands + minerals #########
pls_3 <- df_non_growing_train %>% 
  pls_variable_selection(predictors_3) %>% 
  pls_best_model(df_non_growing_train, c(predictors_3))
# get_q2(df_non_growing_train, pls_3)

# validation against test data
ncomp <- selectNcomp(pls_3, method = "onesigma")
pls_prediction_3 <- predict(pls_3, df_non_growing_test[,predictors_3], ncomp=ncomp)[,,1]

######### vegetation indicies + bands + minerals + topography #########
pls_4 <- df_non_growing_train %>% 
  pls_variable_selection(predictors_4) %>% 
  pls_best_model(df_non_growing_train, c(predictors_4))
# get_q2(df_non_growing_train, pls_4)

# validation against test data
ncomp <- selectNcomp(pls_4, method = "onesigma")
pls_prediction_4 <- predict(pls_4, df_non_growing_test[,predictors_4], ncomp=ncomp)[,,1]

######### comparison #########
# obs vs pred plot
pls_pred_list <- list(pls_prediction_1,pls_prediction_2,pls_prediction_3,pls_prediction_4) 

par(mfrow = c(2,2))
for(pred in pls_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}

# a table for RMSE and MAE
df_pls_rmse <- get_rmse_df(pls_pred_list)

######### Export ############
models <- list(pls_1, pls_2, pls_3, pls_4)
for(i in 1:length(models)){
  temp_filename <- paste0("pls_",i)
  saveRDS(models[[i]], file = paste0("./Script/models/", temp_filename, ".rds"))
}

# rmse df
write.csv(df_pls_rmse, "./Script/csv/pls_rmse.csv")

# plot
png("./Script/figures/pls.png", width = 600, height = 700)
  par(mfrow = c(2,2))
  for(pred in pls_pred_list){
    plot_obs_pred(df_non_growing_test$BareGround,pred)
  }
dev.off()  
