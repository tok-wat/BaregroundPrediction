library(tidyverse)
library(pls)
library(plsdepot)
library(plsVarSel)
library(plsRglm)
library(plsRbeta)
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

######### plsRbeta Analysis ########################
# random seed
set.seed(123)
df_non_growing <- df_ground_photo %>% 
  filter(Month >= 5)

# split in to 70% training and 30% test data
df_non_growing_train <- df_non_growing %>% 
  sample_frac(0.7)
df_non_growing_test <- df_non_growing %>% 
  anti_join(df_non_growing_train, by = "ID") 



####### vegetation indicies only ###############
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
# print(cv_results_2) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_2 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_2, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
# plsRglm::boxplots.bootpls(plsbeta_boot_2)

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

####### vegetation indicies + bands + minerals ###############
# Select predictor variables
predictors_3 <- c(vi, bands, minerals)
df_plsbeta_cv_3 <- df_non_growing_train %>% get_plsbeta_df(predictors_3)

# Perform 10-fold cross-validation for PLS-beta regression
plsbeta_cv_3 <- PLS_beta_kfoldcv_formula(
  BareGround ~ ., data = df_plsbeta_cv_3, nt = 10, 
  modele = "pls-beta", K = 10, verbose = FALSE
)
cv_results_3 <- kfolds2CVinfos_beta(plsbeta_cv_3)
# print(cv_results_3) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_3 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_3, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
# plsRglm::boxplots.bootpls(plsbeta_boot_3)

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

####### vegetation indicies + bands + minerals + topography ###############
# Select predictor variables
predictors_4 <- c(vi, bands, minerals, topography)
df_plsbeta_cv_4 <- df_non_growing_train %>% get_plsbeta_df(predictors_4)

# Perform 10-fold cross-validation for PLS-beta regression
plsbeta_cv_4 <- PLS_beta_kfoldcv_formula(
  BareGround ~ ., data = df_plsbeta_cv_4, nt = 10, 
  modele = "pls-beta", K = 10, verbose = FALSE
)

# Determine the optimal number of components (Nb_comp) from cross-validation results
cv_results_4 <- kfolds2CVinfos_beta(plsbeta_cv_4)
# print(cv_results_4) # Choose Nb_comp=2 based on Q2Chisqcum_Y

# Perform bootstrapping for PLS-beta model with the selected number of components
ncomp <- 2
plsbeta_boot_4 <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_cv_4, nt = ncomp, 
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
# plsRglm::boxplots.bootpls(plsbeta_boot_4)

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


######### comparison #########
# obs vs pred plot
plsbeta_pred_list <- list(plsbeta_prediction_1,plsbeta_prediction_2,plsbeta_prediction_3,plsbeta_prediction_4) 

par(mfrow = c(2,2))
for(pred in plsbeta_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}

# a table for RMSE and MAE
df_pls_rmse <- get_rmse_df(plsbeta_pred_list)

######### Export ############
models <- list(plsbeta_1, plsbeta_2, plsbeta_3, plsbeta_4)
for(i in 1:length(models)){
  temp_filename <- paste0("plsbeta_",i)
  saveRDS(models[[i]], file = paste0("./Script/models/", temp_filename, ".rds"))
}

# rmse df
write.csv(df_pls_rmse, "./Script/csv/plsbeta_rmse.csv")

# plot
png("./Script/figures/plsbeta.png", width = 600, height = 700)
par(mfrow = c(2,2))
for(pred in plsbeta_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}
dev.off()  


