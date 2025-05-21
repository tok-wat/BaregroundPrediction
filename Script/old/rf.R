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
# Data wrangling

# Define predictor variable groups
bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio","tc_bright","tc_green","tc_wet")
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect")


# ------------------------------------------------------------------------------
# Split data by growing season
sf_non_growing <- sf_ground_photo %>% filter(Month >= 5)
sf_growing     <- sf_ground_photo %>% filter(Month < 5)

# non-growing season
sf_non_growing_split <- func_data_split(sf_non_growing, 0.75, k=5, seed=1234) #70:30 split & spatial blocking for CV
# growing season
sf_growing_split <- func_data_split(sf_growing, 0.75, k=5, seed=1234)

# ------------------------------------------------------------------------------
# Run the random forest model pipeline
library(collinear)

func_rf_pipeline <- function(dataframe, predictors, blocks, seed){
  selected_vars <- func_colinear(dataframe, predictors)
  temp_rfe <- func_rfe(
    dataframe = dataframe,
    predictors = selected_vars, 
    blocks = blocks, 
    seed=seed)
  rf_model <- func_caret_rf(
    dataframe = dataframe,
    predictors = temp_rfe$selected_predictors, 
    blocks = blocks, 
    seed=seed)
  return(list(cor_output = selected_vars,
              RFE_output = temp_rfe, 
              best_model = rf_model))
}

func_rf_pipeline(dataframe = sf_non_growing_split$train_df,
                 predictors = c(vi,bands,soil,minerals),
                 blocks = sf_non_growing_split$spatial_cv_block,
                 seed=1234)

func_rf_pipeline(dataframe = sf_growing_split$train_df,
                 predictors = c(vi,bands,soil,minerals),
                 blocks = sf_growing_split$spatial_cv_block,
                 seed=1234)


# Compared to PLS-regression, random forest performs better for non-growing season but 
# doesn't perform well in growing season - perhaps data is not enough
