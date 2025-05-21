# Load required libraries
library(sf)
library(tidyverse)
library(jsonlite)
library(CAST)  # for knndm and CreateSpacetimeFolds
library(caret)
library(ggplot2)

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
# Function to create spatio-temporal CV folds for use in caret::train
# Input: 
#   - sf_object: an sf object with spatial coordinates and a 'Year' column
#   - study_area: an sf polygon representing the study area (same CRS as sf_object)
#   - k: number of folds for cross-validation (default is 5)
# Output:
#   - A list of training indices compatible with caret::trainControl(index = ...)

create_spatiotemporal_index <- function(sf_object, study_area, k = 5) {
  # Check required column
  if (!"Year" %in% colnames(sf_object)) {
    stop("The input sf object must contain a 'Year' column.")
  }
  
  # Step 1: Spatial blocking using KNNDM
  indices_knndm <- knndm(sf_object, study_area, k = k)
  
  # Step 2: Add spatial cluster information to sf object
  sf_object$clusters <- indices_knndm$clusters
  p <- plot(geodist(sf_object, study_area, cvfolds = indices_knndm$indx_test)) +
    scale_x_log10(labels = round) +
    ggtitle("Prediction-to-sample and CV distance distribution") +
    theme_minimal()
  print(p)
  
  # Step 3: Create spatio-temporal folds
  st_fold <- CreateSpacetimeFolds(
    sf_object,
    spacevar = "clusters",
#    timevar = "Year",
    k = k
  )

  # Return the caret-compatible list of training indices
  return(st_fold)
}

sf_non_growing <- sf_ground_photo %>% 
  filter(Month >= 5)
sf_growing <- sf_ground_photo %>% 
  filter(Month < 5) 


# split in to 75% training and 25% test data
sf_non_growing_train <- sf_non_growing %>% 
  sample_frac(0.75)
df_non_growing_valid <- sf_non_growing %>% 
  anti_join(st_drop_geometry(sf_non_growing_train), by = "ID") %>% 
  st_drop_geometry()

sf_growing_train <- sf_growing %>% 
  sample_frac(0.75)
df_growing_valid <- sf_growing %>% 
  anti_join(st_drop_geometry(sf_growing_train), by = "ID") %>% 
  st_drop_geometry()


st_block_non_growing <- sf_non_growing_train %>% 
  create_spatiotemporal_index(study_area, k = 5)

st_block_growing <- sf_growing_train %>% 
  create_spatiotemporal_index(study_area, k = 5)

st_block_non_growing

bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect")

vec_pred_1 <- vi
vec_pred_2 <- vi[! vi %in% c("exg", "bsi")] 
vec_pred_3 <- c(vi,bands) 
vec_pred_4 <- c(vi,bands)[! c(vi,bands) %in% c("exg", "bsi","blue")]



######## 
set.seed(123) 
df_non_growing_train <- st_drop_geometry(sf_non_growing_train)
df_growing_train <- st_drop_geometry(sf_growing_train)



caret_pls_1 <- func_caret_pls(df_non_growing_train, vec_pred_1, st_block_non_growing)
caret_pls_2 <- func_caret_pls(df_non_growing_train, vec_pred_2, st_block_non_growing)
caret_pls_3 <- func_caret_pls(df_non_growing_train, vec_pred_3, st_block_non_growing)
caret_pls_4 <- func_caret_pls(df_non_growing_train, vec_pred_4, st_block_non_growing)

caret_pls_1$results
caret_pls_2$results
caret_pls_3$results
caret_pls_4$results
caret_pls_5$results
caret_pls_6$results
caret_pls_7$results
caret_pls_8$results
caret_pls_9$results

caret_pls_2$bestTune


get_caret_coefficients(caret_pls_3)

####### CLEANED CODE WITH ChatGPT
########## Random Forest Modeling Pipeline with Reproducibility ##########
# Load required libraries
library(collinear)
library(caret)
library(dplyr)

# ------------------------------------------------------------------------------
# Function to select predictors based on Variance Inflation Factor (VIF)
# Input: 
#   - dataframe: a dataframe with response and explanatory variables. 
#                The column name for the response variable must be "BareGround"
#   - predictors: a vector of column names for the set of expalantory variables
# Output:
#   - A vector of selected predictors

func_vif <- function(dataframe, predictors) {
  vif_select(
    df = dataframe,
    predictors = predictors,
    preference_order = c(
      "satvi", "swirRatio", "ndvi", "osavi", "ndmi",
      "evi", "savi", "msavi", "exg", "bsi",
      "red", "green", "blue", "nir", "swir1", "swir2"
    ),
    max_vif = 10,
    quiet = TRUE
  ) %>%
    as.vector()
}

# ------------------------------------------------------------------------------
# Function to perform Recursive Feature Elimination (RFE) with Random Forest
# Input: 
#   - dataframe: a dataframe with response and explanatory variables. 
#                The column name for the response variable must be "BareGround"
#   - predictors: a vector of column names for the set of expalantory variables
#   - blocks: a list of spatial blocks developed from CAST::knndm
#   - seed: random seed
# Output:
#   - A list of
#     [[1]]: a dataframe containing metrics of the models (eg RMSE, R^2) 
#     [[2]]: number of predictors of the simplest model within one sigma from the best model
#     [[3]]: a vector of selected predictors
func_rfe <- function(dataframe, predictors, blocks, seed = 123) {
  set.seed(seed)  # Ensure reproducibility
  
  output <- rfe(
    x = dataframe[, predictors],
    y = dataframe$BareGround,
    sizes = 1:length(predictors),
    rfeControl = rfeControl(
      functions = rfFuncs,
      method = "cv",
      index = blocks$index,
      indexOut = blocks$indexOut
    )
  )
  
  oneSE_nparam <- oneSE(output$results, metric = "RMSE", maximize = FALSE, num = 5)
  selected_pred <- output$optVariables[1:oneSE_nparam]
  
  return(list(
    results = output$results,
    n_selected = oneSE_nparam,
    selected_predictors = selected_pred
  ))
}

# ------------------------------------------------------------------------------
# Function to train Random Forest model using caret
# Input: 
#   - dataframe: a dataframe with response and explanatory variables. 
#                The column name for the response variable must be "BareGround"
#   - predictors: a vector of column names for the set of expalantory variables
#   - blocks: a list of spatial blocks developed from CAST::knndm
#   - seed: random seed
func_caret_rf <- function(dataframe, predictors, blocks, seed = 123) {
  set.seed(seed)  # Ensure reproducibility
  
  caret_model <- train(
    x = dataframe[, predictors],
    y = dataframe$BareGround,
    method = "rf",
    preProcess = c("center", "scale"),
    importance = TRUE,
    trControl = trainControl(
      method = "cv",
      index = blocks$index,
      indexOut = blocks$indexOut,
      savePredictions = "final"
    )
  )
  
  return(caret_model)
}

# Full pipeline function: VIF → RFE → Random Forest
run_rf_pipeline <- function(dataframe, initial_predictors, blocks, seed = 123) {
  selected_vif <- func_vif(dataframe, initial_predictors)
  rfe_result <- func_rfe(dataframe, selected_vif, blocks, seed = seed)
  rf_model <- func_caret_rf(dataframe, rfe_result$selected_predictors, blocks, seed = seed)
  
  return(list(
    vif_predictors = selected_vif,
    rfe_result = rfe_result,
    rf_model = rf_model
  ))
}

# -------- Example Usage --------
result_9 <- run_rf_pipeline(df_non_growing_train, vec_pred_9, st_block_non_growing, seed = 123)

# Access components:
result_9
########################################################





library(collinear)
func_vif <- function(dataframe,predictors){
  vif_select(
      df = dataframe,
      predictors = predictors,
      preference_order = c("satvi","swirRatio","ndvi","osavi","ndmi",
                           "evi", "savi","msavi","exg","bsi",
                           "red","green", "blue","nir","swir1","swir2"),
      #preference_order = c("satvi","swirRatio","savi","msavi"),
      max_vif = 10,
      quiet = TRUE
  ) %>% 
    as.vector()
}

# RFE実行（変数選択）
func_rfe <- function(dataframe, predictors, blocks){
  output <- rfe(
    x = dataframe[,predictors],
    y = dataframe$BareGround,
    sizes = 1:length(predictors),
    rfeControl = rfeControl(
      functions = rfFuncs,
      method = "cv",
      index = blocks$index,
      indexOut = blocks$indexOut)
  )
  oneSE_nparam <- oneSE(output$results, "RMSE", maximize=FALSE, num=5)
  selected_pred <- output$optVariables[1:oneSE_nparam]
  return(list(
    output$results,
    oneSE_nparam,
    selected_pred)
    )
}

func_caret_rf <- function(dataframe, predictors, blocks){
  caret_model <- train(x = dataframe[,predictors],
                       y = dataframe$BareGround,
                       method="rf",
                       preProcess = c('center','scale'),
                       importance=TRUE,
                       trControl=trainControl(method="cv",
                                              index = blocks$index,
                                              indexOut = blocks$indexOut,
                                              savePredictions = "final")
                        )
  return(caret_model)
}


v9 <- func_vif(df_non_growing_train, vec_pred_9)
rfe_9 <- func_rfe(df_non_growing_train, v9, st_block_non_growing)    
rf_9 <- func_caret_rf(df_non_growing_train, rfe_9[[3]],st_block_non_growing)
rf_9


rfe(
  x = df_non_growing_train[,v9],
  y = df_non_growing_train$BareGround,
  sizes = 1:length(v9),
  rfeControl = rfeControl(
    functions = rfFuncs,
    method = "cv",
    index = st_block_non_growing$index,
    indexOut = st_block_non_growing$indexOut)
)

v2 <- func_vif(df_non_growing_train, vec_pred_2)
v5 <- func_vif(df_non_growing_train, vec_pred_5)

v9

func_rf_pls <- function(dataframe, predictors, blocks){
  X <- dataframe[,predictors]
  y <- dataframe$BareGround
  caret_model <- train(y~X
    method="rf",
    preProcess = c('center','scale'),
    importance=TRUE,
    trControl=trainControl(method="cv",
                           index = blocks$index,
                           indexOut = blocks$indexOut,
                           savePredictions = "final")
  )
}



rfe_fit
rfe_fit$optVariables[1:5]

rfe_fit$results %>% 
  data.frame() %>% 
  mutate(bar_min = RMSE-RMSESD) %>%
  mutate(bar_max = RMSE+RMSESD) %>% 
  ggplot(aes(x = Variables, y = RMSE)) +
  geom_errorbar(aes(ymin=bar_min, ymax=bar_max),
                position = position_dodge(width = 1), width = 0, linewidth=1) +
  geom_point(position = position_dodge(width = 1), size=3) +
  geom_hline(aes(
    yintercept=rfe_fit$results[best(rfe_fit$results,"RMSE", maximize=FALSE),2] +
      rfe_fit$results[best(rfe_fit$results,"RMSE", maximize=FALSE),5]/sqrt(5)
    ))+
  theme_bw()




caret_rf_1 <- func_rf_pls(df_non_growing_train, vec_pred_1, st_block_non_growing)
caret_rf_2 <- func_rf_pls(df_non_growing_train, vec_pred_2, st_block_non_growing)
caret_rf_3 <- func_rf_pls(df_non_growing_train, vec_pred_3, st_block_non_growing)
caret_rf_4 <- func_rf_pls(df_non_growing_train, vec_pred_4, st_block_non_growing)
caret_rf_5 <- func_rf_pls(df_non_growing_train, vec_pred_5, st_block_non_growing)
caret_rf_6 <- func_rf_pls(df_non_growing_train, vec_pred_6, st_block_non_growing)
caret_rf_7 <- func_rf_pls(df_non_growing_train, vec_pred_7, st_block_non_growing)
caret_rf_8 <- func_rf_pls(df_non_growing_train, vec_pred_8, st_block_non_growing)
caret_rf_9 <- func_rf_pls(df_non_growing_train, func_vif(df_non_growing_train, vec_pred_9), st_block_non_growing)


caret_rf_1$results
caret_rf_2$results
caret_rf_3$results
caret_rf_4$results
caret_rf_5$results
caret_rf_6$results
caret_rf_7$results
caret_rf_8$results
caret_rf_9$results

##############

caret_growing_pls_1 <- func_caret_pls(df_growing_train, vec_pred_1, st_block_growing)

caret_growing_rf_1 <- func_rf_pls(df_growing_train, vec_pred_1, st_block_growing)
caret_growing_rf_2 <- func_rf_pls(df_growing_train, vec_pred_2, st_block_growing)
caret_growing_rf_3 <- func_rf_pls(df_growing_train, vec_pred_3, st_block_growing)
caret_growing_rf_4 <- func_rf_pls(df_growing_train, vec_pred_4, st_block_growing)
caret_growing_rf_5 <- func_rf_pls(df_growing_train, vec_pred_5, st_block_growing)

caret_growing_pls_1$results

caret_growing_rf_1$results
caret_growing_rf_2$results
caret_growing_rf_3$results
caret_growing_rf_4$results
caret_growing_rf_5$results



###########
pls_func <- function(sf_train,vec_predictors,index){
  df <- st_drop_geometry(sf_train)
  train(df[,vec_predictors],
        df$BareGround,
        method="pls",
        preProcess = c('center','scale'),
        importance=TRUE,
        trControl=trainControl(method="cv",
                               index = index,
                               savePredictions = "final")
  )
}

pls_1 <- pls_func(sf_non_growing_train,vec_pred_1,st_block_non_growing$index)
pls_1$results[pls_1$bestTune[[1]],]

pls_2 <- pls_func(sf_non_growing_train,vec_pred_2,st_block_non_growing$index)
pls_3 <- pls_func(sf_non_growing_train,vec_pred_3,st_block_non_growing$index)
pls_4 <- pls_func(sf_non_growing_train,vec_pred_4,st_block_non_growing$index)
pls_5 <- pls_func(sf_non_growing_train,vec_pred_5,st_block_non_growing$index)
pls_6 <- pls_func(sf_non_growing_train,vec_pred_6,st_block_non_growing$index)
pls_7 <- pls_func(sf_non_growing_train,vec_pred_7,st_block_non_growing$index)
pls_8 <- pls_func(sf_non_growing_train,vec_pred_8,st_block_non_growing$index)
pls_9 <- pls_func(sf_non_growing_train,vec_pred_9,st_block_non_growing$index)

pls_1$results[pls_1$bestTune[[1]],]
pls_2$results[pls_2$bestTune[[1]],]
pls_3$results[pls_3$bestTune[[1]],]
pls_4$results[pls_4$bestTune[[1]],]
pls_5$results[pls_5$bestTune[[1]],]
pls_6$results[pls_6$bestTune[[1]],]
pls_7$results[pls_7$bestTune[[1]],]
pls_8$results[pls_8$bestTune[[1]],]
pls_9$results[pls_9$bestTune[[1]],]



plot(varImp(pls_1))
plot(varImp(pls_2))


pls1_pred <- predict(pls_1, df_growing_valid[,vec_pred_5])
pls3_pred <- predict(pls_3, df_growing_valid[,vec_pred_2])

predictors <- c(vi,bands,topography,minerals,soil)
rf_pred <- predict(model, df_non_growing_valid[,vi])

plot(pls1_pred,pls3_pred)
plot(pls1_pred,rf_pred)
mean_pred <- (pls1_pred+rf_pred)/2
plot(pls3_pred,df_growing_valid$BareGround)
plot(rf_pred,df_non_growing_valid$BareGround)
plot(mean_pred,df_non_growing_valid$BareGround)



bestMtry <- tuneRF(x,y, stepFactor = 1.5, improve = 1e-5, ntree = 500)
model <- train(st_drop_geometry(sf_ground_photo)[,vi],
               st_drop_geometry(sf_ground_photo)$BareGround,
               method="rf",
               tuneGrid=data.frame("mtry"=2), 
               importance=TRUE,
               trControl=trainControl(method="cv",
                                      index = st_fold$index,
                                      savePredictions = "final"))

model

model2 <- train(st_drop_geometry(sf_ground_photo)[,predictors],
               st_drop_geometry(sf_ground_photo)$BareGround,
               method="rf",
               tuneGrid=data.frame("mtry"=2), 
               importance=TRUE,
               trControl=trainControl(method="cv",
                                      index = indices_knndm$indx_train,
                                      savePredictions = "final"))

model2
plot(varImp(model))
plot(varImp(model2))

ffsmodel <- ffs(st_drop_geometry(sf_ground_photo)[,predictors],
                st_drop_geometry(sf_ground_photo)$BareGround,
                method="rf", 
                tuneGrid=data.frame("mtry"=2),
                verbose=FALSE,
                ntree=500, 
                trControl=trainControl(method="cv",
                                       index = st_fold$index,
                                       savePredictions = "final"))
ffsmodel
?ffs
ffsmodel2 <- ffs(st_drop_geometry(sf_ground_photo)[,predictors],
                st_drop_geometry(sf_ground_photo)$BareGround,
                method="rf", 
                tuneGrid=data.frame("mtry"=2),
                verbose=FALSE,
                ntree=500, 
                trControl=trainControl(method="cv",
                                       index = indices_knndm$indx_train,
                                       savePredictions = "final"))
ffsmodel2
plot(varImp(ffsmodel))
plot(varImp(ffsmodel2))
ffs_pred <- predict(ffsmodel, st_drop_geometry(sf_ground_photo[,ffsmodel$selectedvars]))
ffs2_pred <- predict(ffsmodel2, st_drop_geometry(sf_ground_photo[,ffsmodel2$selectedvars]))
model_pred <- predict(model, st_drop_geometry(sf_ground_photo[,predictors]))

plot(model_pred,ffs2_pred)
mean(abs(model_pred-ffs2_pred))
 