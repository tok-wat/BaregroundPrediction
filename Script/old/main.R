library(tidyverse)
library(pls)
library(plsdepot)
library(plsVarSel)
library(plsRglm)
library(randomForest)
library(ggplot2)

############ Load CSV file ##############
df_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 1,006 rows

head(df_ground_photo)

bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
precipitation <- c("X0_2013_precipitation", "X10_2023_precipitation", "X11_2024_precipitation", "X1_2014_precipitation",  
                   "X2_2015_precipitation", "X3_2016_precipitation", "X4_2017_precipitation", "X5_2018_precipitation",  
                   "X6_2019_precipitation", "X7_2020_precipitation", "X8_2021_precipitation", "X9_2022_precipitation" )
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
topography <- c("elevation", "slope", "aspect")

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

  all_variables_model <- plsreg1(df_temp[,predictors],df_temp[,'BareGround'], comps = NULL,  crosval = TRUE)
  out_df <- c(all_variables_model$Q2[2,5], sum(all_variables_model$R2), "all_variables") %>% 
    t() %>% 
    as.data.frame()
  colnames(out_df) <- c("Q2cum", "R2cum", "excluded")
  
  for(i in 1:n_loop){
    excluded_name <- low_vip_names[1]
    low_vip_names <- low_vip_names[-1]
    plsBandsVi <- plsreg1(df_temp[,c(low_vip_names)],df_temp[,'BareGround'], comps = NULL, crosval = TRUE)
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
get_mae_rmse <- function(vec_pred, vec_test){
  rmse <- sqrt(mean((vec_pred - vec_test)^2))
  mae <- mean(abs(vec_pred - vec_test))
  out_vec <- c(mae, rmse)
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
get_mae_rmse(prediction_1, df_non_growing_test$BareGround)

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
get_mae_rmse(prediction_2, df_non_growing_test$BareGround)

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
get_mae_rmse(prediction_3, df_non_growing_test$BareGround)

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
get_mae_rmse(prediction_4, df_non_growing_test$BareGround)
get_coefficients(model_4)

#オブジェクトを保存しよう！
#プロットを書き出そう！
#RMSEをCSVで書き出そう！

######### plsRbeta ########################
library(plsRbeta)
# same as pls regression
# random seed
set.seed(123)
df_non_growing <- df_ground_photo %>% 
  filter(Month >= 5)

# split in to 70% training and 30% test data
df_non_growing_train <- df_non_growing %>% 
  sample_frac(0.7)
df_non_growing_test <- df_non_growing %>% 
  anti_join(df_non_growing_train, by = "ID") 

#### Function
get_plsbeta_df <- function(dataframe, predictors){
  df_temp <- dataframe %>% 
    select(BareGround, all_of(predictors)) %>% 
    mutate(BareGround = BareGround/100) %>% 
    mutate(BareGround = if_else(BareGround == 0, 1e-8, BareGround)) %>% 
    mutate(BareGround = if_else(BareGround == 1, 1-(1e-8), BareGround))
  return(df_temp) 
}
get_selected_var <- function(ci){
  selected_predictors <- (ci[,7]<0&ci[,8]<0)|(ci[,7]>0&ci[,8]>0)
  temp_pred <- selected_predictors[selected_predictors == TRUE] %>% 
    names()
  return(temp_pred)
}
plot_plsbeta <- function(test_vec, pred_vec){
  plot(test_vec, pred_vec,
       xlim =c(0,100),ylim=c(-20,120),
       xlab = "Photo interpretation (%)", ylab = "Model Prediction (%)" , col=NULL, bg=rgb(0, 0, 1, alpha=0.2), pch=21, cex=1.5)
  temp_lm <-lm(pred_vec ~ test_vec)
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  mtext(text=paste0(" R^2 = ",round(summary(temp_lm)$r.squared,3)),adj=0, line=-1) #side=1に対して相対位置を指定（左端が0、右端が1）
}

####### vegetation indicies only ###############
# Select predictor variables
predictors_1 <- c(vi)

# Preprocess training data for beta regression
# Scale BareGround to the range [0,1]
# Adjust 0 and 1 values with 1e-8 correction
df_plsbeta_cv_1 <- df_non_growing_train %>% get_plsbeta_df(predictors_1)

# Perform 10-fold cross-validation for PLS-beta regression
plsbeta_cv_1 <- PLS_beta_kfoldcv_formula(
  BareGround ~ ., data = df_plsbeta_cv_1, nt = 10, 
  modele = "pls-beta", K = 10, verbose = FALSE
)

# Determine the optimal number of components (Nb_comp) from cross-validation results
cv_results_1 <- kfolds2CVinfos_beta(plsbeta_cv_1)
print(cv_results_1) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_1 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_1, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
plsRglm::boxplots.bootpls(plsbeta_boot_1)

# Select significant variables based on confidence intervals
ci_1 <- confints.bootpls(plsbeta_boot_1, indices = 1:length(predictors_1))
selected_vars_1 <- get_selected_var(ci_1) # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
df_pls_train_1 <- df_non_growing_train %>% get_plsbeta_df(selected_vars_1)
plsbeta_1 <- plsRbeta(
  BareGround ~ ., data = df_pls_train_1, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
)

# Preprocess test data for PLS-beta, using only selected variables
df_pls_test_1 <- df_non_growing_test %>% 
  mutate(intercept = 1) %>% 
  select(intercept, all_of(selected_vars_1))

# Compute predictions using the PLS-beta model
logit_values_1 <- as.matrix(df_pls_test_1) %*% as.vector(plsbeta_1$Coeffs)
plsbeta_prediction_1 <- exp(logit_values_1) / (1 + exp(logit_values_1)) * 100

# Evaluate model performance (MAE, RMSE)
get_mae_rmse(plsbeta_prediction_1, df_non_growing_test$BareGround)

####### vegetation indicies + bands ###############
# Select predictor variables
predictors_2 <- c(vi, bands)

# Preprocess training data for beta regression
# Scale BareGround to the range [0,1]
# Adjust 0 and 1 values with 1e-8 correction
df_plsbeta_cv_2 <- df_non_growing_train %>% get_plsbeta_df(predictors_2)

# Perform 10-fold cross-validation for PLS-beta regression
plsbeta_cv_2 <- PLS_beta_kfoldcv_formula(
  BareGround ~ ., data = df_plsbeta_cv_2, nt = 10, 
  modele = "pls-beta", K = 10, verbose = FALSE
)

# Determine the optimal number of components (Nb_comp) from cross-validation results
cv_results_2 <- kfolds2CVinfos_beta(plsbeta_cv_2)
print(cv_results_2) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_2 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_2, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
plsRglm::boxplots.bootpls(plsbeta_boot_2)

# Select significant variables based on confidence intervals
ci_2 <- confints.bootpls(plsbeta_boot_2, indices = 1:length(predictors_2))
selected_vars_2 <- get_selected_var(ci_2) # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
df_pls_train_2 <- df_non_growing_train %>% get_plsbeta_df(selected_vars_2)
plsbeta_2 <- plsRbeta(
  BareGround ~ ., data = df_pls_train_2, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
)

# Preprocess test data for PLS-beta, using only selected variables
df_pls_test_2 <- df_non_growing_test %>% 
  mutate(intercept = 1) %>% 
  select(intercept, all_of(selected_vars_2))

# Compute predictions using the PLS-beta model
logit_values_2 <- as.matrix(df_pls_test_2) %*% as.vector(plsbeta_2$Coeffs)
plsbeta_prediction_2 <- exp(logit_values_2) / (1 + exp(logit_values_2)) * 100

# plot
plot_plsbeta(df_non_growing_test$BareGround,plsbeta_prediction_2)
# Evaluate model performance (MAE, RMSE)
get_mae_rmse(plsbeta_prediction_2, df_non_growing_test$BareGround)

####### vegetation indicies + bands + minerals ###############
# Select predictor variables
predictors_3 <- c(vi, bands, minerals)

# Preprocess training data for beta regression
# Scale BareGround to the range [0,1]
# Adjust 0 and 1 values with 1e-8 correction
df_plsbeta_cv_3 <- df_non_growing_train %>% get_plsbeta_df(predictors_3)

# Perform 10-fold cross-validation for PLS-beta regression
plsbeta_cv_3 <- PLS_beta_kfoldcv_formula(
  BareGround ~ ., data = df_plsbeta_cv_3, nt = 10, 
  modele = "pls-beta", K = 10, verbose = FALSE
)

# Determine the optimal number of components (Nb_comp) from cross-validation results
cv_results_3 <- kfolds2CVinfos_beta(plsbeta_cv_3)
print(cv_results_3) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_3 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_3, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
plsRglm::boxplots.bootpls(plsbeta_boot_3)

# Select significant variables based on confidence intervals
ci_3 <- confints.bootpls(plsbeta_boot_3, indices = 1:length(predictors_3))
selected_vars_3 <- get_selected_var(ci_3) # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
df_pls_train_3 <- df_non_growing_train %>% get_plsbeta_df(selected_vars_3)
plsbeta_3 <- plsRbeta(
  BareGround ~ ., data = df_pls_train_3, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
)

# Preprocess test data for PLS-beta, using only selected variables
df_pls_test_3 <- df_non_growing_test %>% 
  mutate(intercept = 1) %>% 
  select(intercept, all_of(selected_vars_3))

# Compute predictions using the PLS-beta model
logit_values_3 <- as.matrix(df_pls_test_3) %*% as.vector(plsbeta_3$Coeffs)
plsbeta_prediction_3 <- exp(logit_values_3) / (1 + exp(logit_values_3)) * 100

# Evaluate model performance (MAE, RMSE)
get_mae_rmse(plsbeta_prediction_3, df_non_growing_test$BareGround)
# plot
# plot_plsbeta(df_non_growing_test$BareGround, plsbeta_prediction_3)

####### vegetation indicies + bands + minerals + topography ###############
# Select predictor variables
predictors_4 <- c(vi, bands, minerals, topography)

# Preprocess training data for beta regression
# Scale BareGround to the range [0,1]
# Adjust 0 and 1 values with 1e-8 correction
df_plsbeta_cv_4 <- df_non_growing_train %>% get_plsbeta_df(predictors_4)

# Perform 10-fold cross-validation for PLS-beta regression
plsbeta_cv_4 <- PLS_beta_kfoldcv_formula(
  BareGround ~ ., data = df_plsbeta_cv_4, nt = 10, 
  modele = "pls-beta", K = 10, verbose = FALSE
)

# Determine the optimal number of components (Nb_comp) from cross-validation results
cv_results_4 <- kfolds2CVinfos_beta(plsbeta_cv_4)
print(cv_results_4) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_4 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_4, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
plsRglm::boxplots.bootpls(plsbeta_boot_4)

# Select significant variables based on confidence intervals
ci_4 <- confints.bootpls(plsbeta_boot_4, indices = 1:length(predictors_4))
selected_vars_4 <- get_selected_var(ci_4) # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
df_pls_train_4 <- df_non_growing_train %>% get_plsbeta_df(selected_vars_4)
plsbeta_4 <- plsRbeta(
  BareGround ~ ., data = df_pls_train_4, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
)

# Preprocess test data for PLS-beta, using only selected variables
df_pls_test_4 <- df_non_growing_test %>% 
  mutate(intercept = 1) %>% 
  select(intercept, all_of(selected_vars_4))

# Compute predictions using the PLS-beta model
logit_values_4 <- as.matrix(df_pls_test_4) %*% as.vector(plsbeta_4$Coeffs)
plsbeta_prediction_4 <- exp(logit_values_4) / (1 + exp(logit_values_4)) * 100

# Evaluate model performance (MAE, RMSE)
get_mae_rmse(plsbeta_prediction_4, df_non_growing_test$BareGround)
# plot
# plot_plsbeta(df_non_growing_test$BareGround, plsbeta_prediction_4)




######### RandomForest regression models ############
library(caret)

rfe_tune <- function(dataframe, predictors){
  temp_df <- dataframe %>% 
    select(BareGround, all_of(predictors))  
  temp_rfe_tune <- caret::rfe(x = temp_df[,names(temp_df) %in% predictors],
                              y = temp_df$BareGround,
                              sizes = 1:length(predictors),
                              metric = "RMSE",
                              rfeControl = rfeControl(functions =rfFuncs,
                                                      method = "cv",
                                                      repeats = 10, 
                                                      seeds = set.seed(123))
  )
  return (temp_rfe_tune)
}
run_rf <- function(dataframe, predictors, tuned_rfe){
  temp_df <- dataframe %>% 
    select(BareGround, all_of(predictors(tuned_rfe)))
  best <- tuneRF(temp_df[,-1],temp_df$BareGround,doBest=T)
  temp_rf_model <- randomForest(BareGround~., data = temp_df, mtry=best$mtry, ntree=500)
  return(temp_rf_model)
}

####### Vegetation Indicies Only #################
predictors <- c(vi)

rfe_1 <- rfe_tune(df_non_growing_train, predictors) 
rf_1 <- run_rf(df_non_growing_train, predictors, rfe_1)
rf_prediction_1 <- predict(rf_1, df_non_growing_test[,predictors])

# rfe_1$results %>% 
#   ggplot(aes(x=Variables, y=RMSE))+
#   geom_point()+
#   geom_errorbar(aes(ymin=RMSE-RMSESD, ymax=RMSE+RMSESD))+
#   theme_bw()
# 
# varImpPlot(rf_1)

# plot_plsbeta(df_non_growing_test$BareGround, rf_prediction_1)
# get_mae_rmse(rf_prediction_1, df_non_growing_test$BareGround)

####### Vegetation Indicies + bands #################
predictors <- c(vi,bands)

rfe_2 <- rfe_tune(df_non_growing_train, predictors) 
rf_2 <- run_rf(df_non_growing_train, predictors, rfe_2)
rf_prediction_2 <- predict(rf_2, df_non_growing_test[,predictors])

####### Vegetation Indicies + bands + minerals #################
predictors <- c(vi,bands, minerals)

rfe_3 <- rfe_tune(df_non_growing_train, predictors) 
rf_3 <- run_rf(df_non_growing_train, predictors, rfe_3)
rf_prediction_3 <- predict(rf_3, df_non_growing_test[,predictors])

####### Vegetation Indicies + bands + minerals + topography #################
predictors <- c(vi,bands, minerals, topography)

rfe_4 <- rfe_tune(df_non_growing_train, predictors) 
rf_4 <- run_rf(df_non_growing_train, predictors, rfe_4)
rf_prediction_4 <- predict(rf_4, df_non_growing_test[,predictors])




###### comparison #########

pls_pred_list <- list(prediction_1,prediction_2,prediction_3,prediction_4) 
plsbeta_pred_list <- list(plsbeta_prediction_1,plsbeta_prediction_2,plsbeta_prediction_3,plsbeta_prediction_4) 
rf_pred_list <- list(rf_prediction_1,rf_prediction_2,rf_prediction_3,rf_prediction_4) 

get_rmse_df <- function(pred_list){
  rmse_df <- data.frame(MAE=as.numeric(), RMSE=as.numeric())
  for(pred in pred_list){
    temp_rmse <- get_mae_rmse(pred, df_non_growing_test$BareGround) %>% 
      t() %>% 
      as.data.frame()
    rmse_df <- rbind(rmse_df, temp_rmse)
  }
  return(rmse_df)
}

get_rmse_df(pls_pred_list)
get_rmse_df(plsbeta_pred_list)
get_rmse_df(rf_pred_list)

# plot
par(mfrow = c(2,2))

for(pred in pls_pred_list){
  plot_plsbeta(df_non_growing_test$BareGround,pred)
}

for(pred in plsbeta_pred_list){
  plot_plsbeta(df_non_growing_test$BareGround,pred)
}

for(pred in rf_pred_list){
  plot_plsbeta(df_non_growing_test$BareGround,pred)
}



###### test against field data ######
df_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 1,006 rows
