library(tidyverse)
library(caret)
library(car)
library(randomForest)
library(blockCV)
library(jsonlite)
library(ggplot2)

############ Load CSV file ##############
df_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>%
  mutate(
    geo_parsed = map(.geo, ~ fromJSON(.)$coordinates),
    longitude = map_dbl(geo_parsed, 1),
    latitude = map_dbl(geo_parsed, 2)
  ) %>%
  select(-geo_parsed) %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 1,017 rows

head(df_ground_photo)

# Load a script
source("./Script/functions.R")

bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi", "ndmi","ndvi","osavi","satvi","savi","swirRatio")
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect")

predictors_1 <- c(vi)
predictors_2 <- c(vi,soil,topography)
predictors_all <- c(vi, bands, minerals, soil, topography)

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

# spatial blocking for cross validation
n_kfold <- 5
spatial_block <- df_non_growing_train %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  cv_spatial(column = "BareGround", # the response column (binary or multi-class)
             k = n_kfold, # number of folds
             size = 100000, # size of the blocks in metres
             selection = "random", # random blocks-to-fold
             iteration = 50 # find evenly dispersed folds
  )  

cv_plot(cv = spatial_block, 
        x = df_non_growing_train %>% 
          sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326))

list_spatial_block <- list(
  "V1"=spatial_block$folds_list[[1]][[2]],
  "V2"=spatial_block$folds_list[[2]][[2]],
  "V3"=spatial_block$folds_list[[3]][[2]],
  "V4"=spatial_block$folds_list[[4]][[2]],
  "V5"=spatial_block$folds_list[[5]][[2]]
)

############ Eliminating highly correlated predictors #######
library(collinear)
selected <- vif_select(
  df = df_non_growing_train,
  predictors = predictors_2,
#  preference_order = c("satvi","swirRatio","savi","msavi","nir","red"),
  preference_order = c("satvi","swirRatio","savi","msavi"),
  max_vif = 10,
  quiet = FALSE
) %>% 
  as.vector()

###### randomforest
ctrl <- trainControl(
  method = "cv",
  number = n_kfold,
  index = list_spatial_block,
  savePredictions = "final",
  allowParallel = TRUE
)

# rfeControl 設定
rfe_ctrl <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  number = n_kfold,
  index = list_spatial_block
)

##############################################################
# 保存用のリスト
n_boot <- 100  # 通常は100回以上が望ましい
results <- vector("list", n_boot)
results_rmse <- vector("list", n_boot)
results_mae <- vector("list", n_boot)
results_r2 <- vector("list", n_boot)

for (i in 1:n_boot) {
  set.seed(100 + i)
  
  # ブートストラップサンプル作成
  boot_index <- sample(1:nrow(df_non_growing_train), replace = TRUE)
  x_boot <- df_non_growing_train[boot_index, selected]
  y_boot <- df_non_growing_train[boot_index, "BareGround"]
  
  # RFEの設定
  rfe_ctrl <- rfeControl(
    functions = rfFuncs,
    method = "cv",
    number = n_kfold,
    index = list_spatial_block
  )
  
  
  # RFE実行（変数選択）
  rfe_fit <- rfe(
    x = x_boot,
    y = y_boot,
    sizes = 1:length(selected),
    rfeControl = rfe_ctrl
  )
  
  # 結果保存
  results[[i]] <- list(
    selected_vars = predictors(rfe_fit))
  results_rmse[[i]] <- list(rmse = rfe_fit$results$RMSE)
  results_mae[[i]] <- list(mae = rfe_fit$results$MAE)
  results_r2[[i]] <- list (r2 = rfe_fit$results$Rsquared)
}

# 結果をデータフレームにまとめる
df_results <- bind_rows(results)
barplot(sort(table(df_results)))

df_rmse <- do.call(cbind, lapply(results_rmse, function(x) x$rmse)) %>% 
  as.data.frame()
df_mae <- do.call(cbind, lapply(results_mae, function(x) x$mae)) %>% 
  as.data.frame()
df_r2 <- do.call(cbind, lapply(results_r2, function(x) x$r2)) %>% 
  as.data.frame()

func_plot_metrics <- function(dataframe){
  df_plot <- cbind(apply(dataframe, 1, mean),
                   apply(dataframe, 1, quantile, probs=0.95),
                   apply(dataframe, 1, quantile, probs=0.05))
  colnames(df_plot)<-c("mean","upper","lower")
  df_plot %>% 
    ggplot(aes(x = 1:nrow(df_plot), y = mean)) +
      geom_errorbar(aes(ymin=lower, ymax=upper),
                    position = position_dodge(width = 1), width = 0, linewidth=1) +
      geom_point(aes(x = 1:nrow(df_plot), y =mean), position = position_dodge(width = 1), size=3) +
      theme_bw()
}

func_plot_metrics(df_rmse)
func_plot_metrics(df_mae)
func_plot_metrics(df_r2)

# 最も重要とされた変数の名前を取得
n_predictors <- 9 #
selected_vars <- table(df_results) %>% 
  data.frame() %>% 
  arrange(desc(Freq)) %>% 
  head(n_predictors) %>% 
  select(selected_vars) %>% 
  unlist() %>% 
  as.character()

final_model <- train(
  x = df_non_growing_train[, selected_vars],
  y = df_non_growing_train$BareGround,
  method = "rf",
  trControl = ctrl,
  importance = TRUE
)

# 結果を表示
print(final_model)

# 変数重要度を確認
varImp(final_model)
plot(varImp(final_model))


############ Load CSV file ##############
df_transect <- read.csv("./CSV/transect_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 52 rows


rf_prediction_holdout <- predict(final_model, df_non_growing_test[,selected_vars])
rf_prediction_transect <- predict(final_model, df_transect[,selected_vars])


plot_obs_pred(df_non_growing_test$BareGround,rf_prediction_holdout)
plot_obs_pred(df_transect$BareGround,rf_prediction_transect)

# a table for RMSE and MAE
get_rmse_df(list(rf_prediction_holdout), df_non_growing_test$BareGround)
get_rmse_df(list(rf_prediction_transect), df_transect$BareGround)
 