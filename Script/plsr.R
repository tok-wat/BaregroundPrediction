# Load required libraries
library(sf)
library(tidyverse)
library(jsonlite)
library(MASS)  # for box-cox transformation
library(CAST)  # for knndm and CreateSpacetimeFolds
library(caret)
library(ggplot2)
library(ggforce)
library(RColorBrewer)

##################################
#setwd("C:/Users/turat/Dropbox/01_UCT/04_DegradationMapping/07_Analysis/01_Estimate Cover from Ground Level Photo/30_PLSR_Model/BaregroundPrediction")
setwd("~/01_UCT/04_Degradation Mapping/07_Analysis/01_Estimate Cover from Ground Level Photo/30_PLSR_Model/BaregroundPrediction")

# Load a script
source("./Script/functions.R")

# ------------------------------------------------------------------------------
# Load study area shapefile and transform to UTM Zone 34S (EPSG:32734)
study_area <- st_read("./SHP/StudyArea2024.shp") %>% 
  st_transform(crs = 32734)

# ------------------------------------------------------------------------------
# Define predictor variable groups
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio","tc_bright","tc_green","tc_wet")

# Load ground-level photo metadata and convert to spatial data frame
sf_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>%
  mutate(
    # Parse JSON from .geo column to extract coordinates
    geo_parsed = map(.geo, ~ fromJSON(.)$coordinates),
    longitude = map_dbl(geo_parsed, 1),
    latitude = map_dbl(geo_parsed, 2)
  ) %>%
  dplyr::select(ID, longitude, latitude, BioRegion, Biome, Year, Month, 
                BareGround, all_of(vi)) %>%  # Drop unnecessary columns
  na.omit() %>%  # Remove rows with missing values (resulting in 1,017 rows)
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%  # Convert to sf object in WGS84
  st_transform(crs = 32734)  # Transform to UTM for spatial blocking

# ------------------------------------------------------------------------------
# Split data by growing season
sf_growing     <- sf_ground_photo %>% filter(Month < 5)
sf_non_growing <- sf_ground_photo %>% filter(Month >= 5)

# Selected Vegetation Indices to Apply Box-Cox transformation
growing_transform_vi <- c("bsi","evi","exg","ndmi","ndvi","swirRatio","tc_bright","tc_wet")
non_growing_transform_vi <- c("exg", "swirRatio","tc_wet")
# c("exg","ndmi", "swirRatio","tc_bright","tc_wet")

# growing_transform_viの各列に対してパラメータ取得
growing_params <- map(growing_transform_vi, function(col) {
  params <- get_boxcox_params(sf_growing[[col]], sf_growing$BareGround)
  tibble(variable = col, lambda = params$lambda, shift_x = params$shift_x)
}) %>% 
  list_rbind()

# 変換列を追加
sf_growing <- sf_growing %>%
  mutate(across(all_of(growing_transform_vi), 
                .fns = function(col) {
                  p <- growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                }, 
                .names = "{.col}_boxcox" 
                ))

# non-growing_transform_viの各列に対してパラメータ取得
non_growing_params <- map(non_growing_transform_vi, function(col) {
  params <- get_boxcox_params(sf_non_growing[[col]], sf_non_growing$BareGround)
  tibble(variable = col, lambda = params$lambda, shift_x = params$shift_x)
}) %>% 
  list_rbind()

# 変換列を追加
sf_non_growing <- sf_non_growing %>%
  mutate(across(all_of(non_growing_transform_vi), 
                .fns = function(col) {
                  p <- non_growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                }, 
                .names = "{.col}_boxcox" 
  ))
# ------------------------------------------------------------------------------
# Run the pls model pipeline
library(plsVarSel)

#### 1. growing season #########################################################
sf_growing_split <- func_data_split(sf_growing, 0.8, k=5, seed=123) #80:20 split & spatial blocking for CV
growing_tr_vi <- c(setdiff(vi, growing_transform_vi), paste0(growing_transform_vi, "_boxcox"))

#### 1-1 growing season without transformation #################################
full_growing <- caret::train(
  x=sf_growing_split$train_df[,vi], 
  y=sf_growing_split$train_df$BareGround,
  method="pls",
  preProcess = c('center', 'scale'),        # Standardize predictors before training
  importance = TRUE,                        # Enable variable importance calculation
  metric = "Qsquared",
  trControl = trainControl(
    method = "cv",                          # Cross-validation
    summaryFunction = func_q2,
    selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
    index = sf_growing_split$spatial_cv_block$indx_train, 
    indexOut = sf_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
    savePredictions = "final"               # Save final predictions
  ))

growing_results <- func_pls_var_selection(sf_growing_split, 
                                          vi, 
                                          full_growing$bestTune$ncomp)

#### 1-2 growing season with transformation ####################################
full_growing_tr <- caret::train(
  x=sf_growing_split$train_df[,growing_tr_vi], 
  y=sf_growing_split$train_df$BareGround,
  method="pls",
  preProcess = c('center', 'scale'),        # Standardize predictors before training
  importance = TRUE,                        # Enable variable importance calculation
  metric = "Qsquared",
  trControl = trainControl(
    method = "cv",                          # Cross-validation
    summaryFunction = func_q2,
    selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
    index = sf_growing_split$spatial_cv_block$indx_train, 
    indexOut = sf_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
    savePredictions = "final"               # Save final predictions
  ))

growing_tr_results <- func_pls_var_selection(sf_growing_split, 
                                          growing_tr_vi, 
                                          full_growing$bestTune$ncomp)

#### 2. non-growing season #####################################################
sf_non_growing_split <- func_data_split(sf_non_growing, 0.8, k=5, seed=123) #80:20 split & spatial blocking for CV
non_growing_tr_vi <- c(setdiff(vi, non_growing_transform_vi), paste0(non_growing_transform_vi, "_boxcox"))

#### 2-1 non-growing season without transformation #############################
# run a full model
full_non_growing <- caret::train(
                          x=sf_non_growing_split$train_df[,vi], 
                          y=sf_non_growing_split$train_df$BareGround,
                          method="pls",
                          preProcess = c('center', 'scale'),        # Standardize predictors before training
                          importance = TRUE,                        # Enable variable importance calculation
                          metric = "Qsquared",
                          trControl = trainControl(
                            method = "cv",                          # Cross-validation
                            summaryFunction = func_q2,
                            selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
                            index = sf_non_growing_split$spatial_cv_block$indx_train, 
                            indexOut = sf_non_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
                            savePredictions = "final"               # Save final predictions
                          ))
# full_non_growing
non_growing_results <- func_pls_var_selection(sf_non_growing_split, 
                                              vi, 
                                              full_non_growing$bestTune$ncomp)

#### 2-2 non-growing season with Box-Cox transformation ########################
# run a full model
full_non_growing_tr <- caret::train(
  x=sf_non_growing_split$train_df[,non_growing_tr_vi], 
  y=sf_non_growing_split$train_df$BareGround,
  method="pls",
  preProcess = c('center', 'scale'),        # Standardize predictors before training
  importance = TRUE,                        # Enable variable importance calculation
  metric = "Qsquared",
  trControl = trainControl(
    method = "cv",                          # Cross-validation
    summaryFunction = func_q2,
    selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
    index = sf_non_growing_split$spatial_cv_block$indx_train, 
    indexOut = sf_non_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
    savePredictions = "final"               # Save final predictions
  ))

non_growing_tr_results <- func_pls_var_selection(sf_non_growing_split, 
                                                 non_growing_tr_vi, 
                                                 full_non_growing$bestTune$ncomp)

#### 3. Results ################################################################

## Just to have a look
# growing_results
# growing_tr_results
# non_growing_results
# non_growing_tr_results

#### 3-1 variable selection #################################################### 
# backward elimination based on the order of VIP in the full models
# Variables with VIP > 1.0 are retained as they are generally considered useful in prediction
# Best model is selected where Q2 is the highest

growing_summary <- rbind(
  growing_results$metrics %>% 
    slice(1:sum(growing_results$vip_full_model>=1.0)) %>% # number of variables VIP >1.0
    mutate(model = "vi"),
  growing_tr_results$metrics %>% 
     slice(1:sum(growing_tr_results$vip_full_model>=1.0)) %>% 
     mutate(model="transformed")
  )

non_growing_summary <- rbind(
  non_growing_results$metrics %>% 
    slice(1:sum(non_growing_results$vip_full_model>=1.0)) %>% # number of variables VIP >1.0
    mutate(model = "vi"),
  non_growing_tr_results$metrics %>% 
    slice(1:sum(non_growing_tr_results$vip_full_model>=1.0)) %>% 
    mutate(model="transformed")
)

write.csv(non_growing_summary, "./CSV/Non-growing_var_selection.csv")
write.csv(growing_summary, "./CSV/Growing_var_selection.csv")

#### 3-2 Identify the best set of variables in different models ################ 
top_growing <- growing_summary %>% 
  arrange(desc(Qsquared)) %>% 
  head()

top_non_growing <- non_growing_summary %>% 
  arrange(desc(Qsquared)) %>% 
  head()

top_growing
top_non_growing

#### 3-3 Run the best models again and save them ###############################
best_growing <- func_caret_pls(dataframe=sf_growing_split$train_df,
                               predictors = str_split(top_growing[1,]$variables, pattern = ",\\s*")[[1]],
                               blocks=sf_growing_split$spatial_cv_block,
                               n_comp=full_growing$bestTune$ncomp)

best_non_growing <- func_caret_pls(dataframe=sf_non_growing_split$train_df,
                                   predictors = str_split(top_non_growing[1,]$variables, pattern = ",\\s*")[[1]],
                                   blocks=sf_non_growing_split$spatial_cv_block,
                                   n_comp=full_non_growing$bestTune$ncomp)

# save models
saveRDS(best_non_growing, file = "./Script/models/pls_non_growing.rds")
saveRDS(best_growing, file = "./Script/models/pls_growing.rds")

# extract coefficients and parameters for box-cox transformation to predict in GEE
get_caret_coefficients(best_growing)
# growing_params

get_caret_coefficients(best_non_growing)


#### 4. Validation #############################################################
# ------------------------------------------------------------------------------
# Photo-Model Validation
valid_growing <- predict(best_growing, sf_growing_split$valid_df)
plot_obs_pred(sf_growing_split$valid_df$BareGround, valid_growing)

valid_non_growing <- predict(best_non_growing, sf_non_growing_split$valid_df)
plot_obs_pred(sf_non_growing_split$valid_df$BareGround, valid_non_growing)

# ------------------------------------------------------------------------------
# Transect Data Validation
df_transect <- read.csv("./CSV/transect_extracted.csv")

df_transect_growing <- df_transect %>% 
  filter(Month < 5) %>% 
  dplyr::select(BareGround,all_of(vi)) %>%
  na.omit()  # 48 rows

df_transect_non_growing <- df_transect %>% 
  filter(Month >= 5) %>% 
  dplyr::select(BareGround,all_of(vi)) %>%
  na.omit()  # 52 rows

# growing season
# Box-cox transformatioin
df_transect_growing <- df_transect_growing %>%
  mutate(across(all_of(growing_transform_vi), 
                .fns = function(col) {
                  p <- growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                }, 
                .names = "{.col}_boxcox" 
  ))

transect_growing <- predict(best_growing, df_transect_growing)
plot_obs_pred(df_transect_growing$BareGround, transect_growing)

# non-growing season
transect_non_growing <- predict(best_non_growing, df_transect_non_growing)
plot_obs_pred(df_transect_non_growing$BareGround, transect_non_growing)
