library(tidyverse)
library(pls)
library(plsdepot)
library(plsVarSel)
library(randomForest)

############ Load CSV file ##############
df_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit() %>%  # 1,006 rows
  mutate(cos_aspect = cos(aspect))
head(df_ground_photo)

bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
precipitation <- c("X0_2013_precipitation", "X10_2023_precipitation", "X11_2024_precipitation", "X1_2014_precipitation",  
                   "X2_2015_precipitation", "X3_2016_precipitation", "X4_2017_precipitation", "X5_2018_precipitation",  
                   "X6_2019_precipitation", "X7_2020_precipitation", "X8_2021_precipitation", "X9_2022_precipitation" )
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
topography <- c("elevation", "slope", "cos_aspect")

############ Data Cleaning ##############

############ Functions ##############
# variable selection based on vip value
pls_variable_selection <- function(dataframe, predictors){
  df_temp <- dataframe[,c('BareGround', predictors)]
  pls_initial <- plsr(BareGround ~., data = df_temp, validation="CV", scale = TRUE)
  low_vip_names <-VIP(pls_initial,2) %>% 
    sort() %>% 
    names()
  
  out_df <- data.frame()
  less_than_median <- VIP(pls_initial,2) < median(VIP(pls_initial,2))
  n_loop <- less_than_median %>% 
    unname() %>% 
    as.integer() %>% 
    sum()

  all_variables_model <- plsreg1(df_temp[,predictors],df_temp[,'BareGround'], crosval = TRUE)
  out_df <- c(all_variables_model$Q2[2,5], sum(all_variables_model$R2), "all_variables") %>% 
    t() %>% 
    as.data.frame()
  colnames(out_df) <- c("Q2cum", "R2cum", "excluded")
  
  for(i in 1:n_loop){
    excluded_name <- low_vip_names[1]
    low_vip_names <- low_vip_names[-1]
    plsBandsVi <- plsreg1(df_temp[,c(low_vip_names)],df_temp[,'BareGround'], crosval = TRUE)
    temp_Q2cum <- plsBandsVi$Q2[2,5]
    temp_R2sum <- sum(plsBandsVi$R2)
    output <- c(temp_Q2cum, temp_R2sum, excluded_name)
    out_df <- rbind(out_df, output)
  }
  return(out_df)
}
pls_best_model <- function(pls_variable_selection_output, dataframe, predictor){
  max_index <- which.max(pls_variable_selection_output$Q2cum)
  removed_predictor <- pls_variable_selection_output %>% 
    slice(2:max_index) %>% 
    pull(excluded)
  selected_predictor <- predictor[!predictor %in% removed_predictor]
  best_model <- plsr(BareGround ~., data = dataframe[,c('BareGround', selected_predictor)], validation="CV", scale = TRUE)
  return(best_model)
}
# extract q2 using plsreg1
get_q2 <- function(dataframe, pls_model){
  temp_model <- plsreg1(dataframe[,c(colnames(pls_model$model)[-1])],dataframe[,'BareGround'], crosval = TRUE)
  return(temp_model$Q2)
}
get_coefficients <- function(pls_model){
  ncomp <- selectNcomp(pls_model, method = "onesigma")
  std_coef <- coef(pls_model, type="original", ncomp=ncomp, intercept=TRUE)
  unstd_coef <- c(std_coef[1], std_coef[2:nrow(std_coef)]/pls_model$scale)
  return(unstd_coef)
}
# evaluate prediction accuracy
get_mae_rsme <- function(vec_pred, vec_test){
  rsme <- sqrt(mean((vec_pred - vec_test)^2))
  mae <- mean(abs(vec_pred - vec_test))
  out_vec <- c(mae, rsme)
  names(out_vec) <- c("MAE", "RMSE")
  return(out_vec)
}
############ Analysis ##############
# random seed
set.seed(123)

df_non_growing <- df_ground_photo %>% 
  filter(Month >= 5)

# split in to 70% training and 30% test data
df_non_growing_train <- df_non_growing %>% 
  sample_frac(0.7)
df_non_growing_test <- df_non_growing %>% 
  anti_join(df_non_growing_train, by = "ID") 

# Run the PLS regression models

######### vegetation indicies only #########
predictors <- c(vi)
model_1 <- df_non_growing_train %>% 
  pls_variable_selection(predictors) %>% 
  pls_best_model(df_non_growing_train, c(predictors))
get_q2(df_non_growing_train, model_1)

# validation against test data
ncomp <- selectNcomp(model_1, method = "onesigma")
prediction_1 <- predict(model_1, df_non_growing_test[,predictors], ncomp=ncomp)[,,1]

plot(df_non_growing_test$BareGround,prediction_1)
get_mae_rsme(prediction_1, df_non_growing_test$BareGround)

######### vegetation indicies + bands #########
predictors <- c(vi, bands)
model_2 <- df_non_growing_train %>% 
  pls_variable_selection(predictors) %>% 
  pls_best_model(df_non_growing_train, c(predictors))
get_q2(df_non_growing_train, model_2)

# validation against test data
ncomp <- selectNcomp(model_2, method = "onesigma")
prediction_2 <- predict(model_2, df_non_growing_test[,predictors], ncomp=ncomp)[,,1]

plot(df_non_growing_test$BareGround,prediction_2)
get_mae_rsme(prediction_2, df_non_growing_test$BareGround)

######### vegetation indicies + bands + minerals #########
predictors <- c(vi, bands, minerals)
model_3 <- df_non_growing_train %>% 
  pls_variable_selection(predictors) %>% 
  pls_best_model(df_non_growing_train, c(predictors))
get_q2(df_non_growing_train, model_3)

# validation against test data
ncomp <- selectNcomp(model_3, method = "onesigma")
prediction_3 <- predict(model_3, df_non_growing_test[,predictors], ncomp=ncomp)[,,1]

plot(df_non_growing_test$BareGround,prediction_3)
get_mae_rsme(prediction_3, df_non_growing_test$BareGround)

######### vegetation indicies + bands + minerals + topography #########
predictors <- c(vi, bands, minerals, topography)
model_4 <- df_non_growing_train %>% 
  pls_variable_selection(predictors) %>% 
  pls_best_model(df_non_growing_train, c(predictors))
get_q2(df_non_growing_train, model_4)

# validation against test data
ncomp <- selectNcomp(model_4, method = "onesigma")
prediction_4 <- predict(model_4, df_non_growing_test[,predictors], ncomp=ncomp)[,,1]

plot(df_non_growing_test$BareGround,prediction_4)
get_mae_rsme(prediction_4, df_non_growing_test$BareGround)



get_coefficients(model_4)


# Run the RandomForest regression models