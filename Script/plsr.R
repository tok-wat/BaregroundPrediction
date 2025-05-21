# Load required libraries
library(sf)
library(tidyverse)
library(jsonlite)
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
# Load ground-level photo metadata and convert to spatial data frame
sf_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>%
  mutate(
    # Parse JSON from .geo column to extract coordinates
    geo_parsed = map(.geo, ~ fromJSON(.)$coordinates),
    longitude = map_dbl(geo_parsed, 1),
    latitude = map_dbl(geo_parsed, 2)
  ) %>%
  dplyr::select(-geo_parsed, -system.index, -.geo, -vegType) %>%  # Drop unnecessary columns
  na.omit() %>%  # Remove rows with missing values (resulting in 1,017 rows)
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%  # Convert to sf object in WGS84
  st_transform(crs = 32734)  # Transform to UTM for spatial blocking

# ------------------------------------------------------------------------------
# Data wrangling

# Define predictor variable groups
bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio","tc_bright","tc_green","tc_wet")
# minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
# soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
# topography <- c("elevation", "slope", "aspect")


# ------------------------------------------------------------------------------
# Split data by growing season
sf_non_growing <- sf_ground_photo %>% filter(Month >= 5)
sf_growing     <- sf_ground_photo %>% filter(Month < 5)

# # Apply log-transformation
# sf_non_growing_log <- func_log_transform(sf_non_growing, c(vi, bands))[[3]]
# sf_growing_log <- func_log_transform(sf_growing, c(vi, bands))[[3]]

# ------------------------------------------------------------------------------
# Run the pls model pipeline
library(plsVarSel)

# non-growing season
sf_non_growing_split <- func_data_split(sf_non_growing, 0.8, k=5, seed=1234) #80:20 split & spatial blocking for CV

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
full_non_growing

non_growing_results <- func_pls_var_selection(sf_non_growing_split, vi, full_non_growing$bestTune$ncomp)

############## lower accuracy ###################
# sf_non_growing_log_split <- func_data_split(sf_non_growing_log, 0.8, k=5, seed=1234)
# non_growing_log_results <- func_pls_var_selection(sf_non_growing_log_split, vi, 2)
# non_growing_log_results 
 
non_growing_bands_results <- func_pls_var_selection(sf_non_growing_split, c(vi,bands),2) # includes optical bands


# growing season
sf_growing_split <- func_data_split(sf_growing, 0.8, k=5, seed=1234)
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
full_growing
growing_results <- func_pls_var_selection(sf_growing_split, vi, full_growing$bestTune$ncomp)

############### lower accuracy ###################
# sf_growing_split <- func_data_split(sf_growing, 0.8, k=5, seed=1234)
# growing_results <- func_pls_var_selection(sf_growing_split, vi,2)
#  
# growing_log_bands_results <- func_pls_var_selection(sf_growing_log_split, c(vi, bands),2)
 
 
# results
non_growing_results
# non_growing_log_results
# non_growing_bands_results
# 
growing_results
# growing_log_results
# growing_log_bands_results

non_growing_summary <- rbind(
  non_growing_results$metrics %>% 
    slice(1:8) %>% 
    mutate(model = "vi")#,
  # non_growing_log_results$metrics %>% 
  #   slice(1:8) %>% 
  #   mutate(model="vi_log"),
  # non_growing_bands_results$metrics %>% 
  #   slice(1:9) %>% 
  #   mutate(model="bands")
)

growing_summary <- rbind(
  growing_results$metrics %>% 
    slice(1:8) %>% 
    mutate(model = "vi")#,
  # growing_log_results$metrics %>%
  #   slice(1:8) %>%
  #   mutate(model="vi_log"),
  # growing_log_bands_results$metrics %>%
  #   slice(1:9) %>%
  #   mutate(model="bands")
)

write.csv(non_growing_summary, "./CSV/Non-growing_var_selection.csv")
write.csv(growing_summary, "./CSV/Growing_var_selection.csv")

# Select variables based on the Q2 values
top_non_growing <- non_growing_summary %>% 
  arrange(desc(Qsquared)) %>% 
  head()

top_growing <- growing_summary %>% 
  arrange(desc(Qsquared)) %>% 
  head()

top_non_growing
top_growing

best_non_growing <- func_caret_pls(dataframe=sf_non_growing_split$train_df,
                                   predictors = str_split(top_non_growing[1,]$variables, pattern = ",\\s*")[[1]],
                                   blocks=sf_non_growing_split$spatial_cv_block,
                                   n_comp=full_non_growing$bestTune$ncomp)

best_growing <- func_caret_pls(dataframe=sf_growing_split$train_df,
                               predictors = str_split(top_growing[1,]$variables, pattern = ",\\s*")[[1]],
                               blocks=sf_growing_split$spatial_cv_block,
                               n_comp=full_growing$bestTune$ncomp)

# save models
saveRDS(best_non_growing, file = "./Script/models/pls_non_growing.rds")
saveRDS(best_growing, file = "./Script/models/pls_growing.rds")

# extract coefficients to predict in GEE
get_caret_coefficients(best_non_growing)
get_caret_coefficients(best_growing)

# func_log_transform(sf_growing, c(vi, bands))[[2]] # log-transformed variables
# func_log_transform(sf_growing, c(vi, bands))[[4]] # value to add to the negative VI: None of them

################### Validation #################################################
# ------------------------------------------------------------------------------
# Model Validation
valid_non_growing <- predict(best_non_growing, sf_non_growing_split$valid_df)
plot_obs_pred(sf_non_growing_split$valid_df$BareGround, valid_non_growing)

valid_growing <- predict(best_growing, sf_growing_split$valid_df)
plot_obs_pred(sf_growing_split$valid_df$BareGround, valid_growing)

# ------------------------------------------------------------------------------
# Transect Data Validation
df_transect <- read.csv("./CSV/transect_extracted.csv")

df_transect_non_growing <-   df_transect %>% 
  filter(Month >= 5) %>% 
  select(BareGround,all_of(vi)) %>%
  na.omit()  # 52 rows

df_transect_growing <-  df_transect %>% 
  filter(Month < 5) %>% 
  select(BareGround,all_of(vi)) %>%
  na.omit()  # 48 rows

# non-growing season
transect_non_growing <- predict(best_non_growing, df_transect_non_growing)
plot_obs_pred(df_transect_non_growing$BareGround, transect_non_growing)

# growing season
transect_growing <- predict(best_growing, df_transect_growing)
plot_obs_pred(df_transect_growing$BareGround, transect_growing)


# forget about log transformation for now
# log_info <- func_log_transform(sf_growing, vi)
# trans_vars <- log_info[[2]]
# add_values <- log_info[[4]]
# 
# df_transect_growing_log <- df_transect_growing %>%
#   mutate(across(all_of(names(add_values)), ~ .x + add_values[cur_column()])) %>%
#   mutate(across(all_of(trans_vars), log))
# 
