library(tidyverse)
library(jsonlite)
library(ggforce)
library(pls)
library(plsdepot)
library(plsVarSel)
library(blockCV)

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
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect")

predictors_1 <- c(vi)
predictors_all <- c(vi, bands, minerals, soil, topography)

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

######### Run the PLS regression models ####
######### vegetation indicies only #########
df_pls_1 <- df_non_growing_train %>% 
  select(BareGround,all_of(predictors_1))  
pls_vi_all <- plsr(BareGround~., data=df_pls_1, validation="CV",scale=TRUE, segments=list_spatial_block) 
summary(pls_vi_all)
plsVarSel::VIP(pls_vi_all,2)

# check VIP
barplot(plsVarSel::VIP(pls_vi_all,2))

func_q2 <- function(pls_model){
  press <- pls_model$validation$PRESS
  tss <- sum((pls_model$model$BareGround-pls_model$Ymeans)^2)
  return(1-press/tss)
}

func_eliminate_predictor <- function(dataframe, elim_predictors){
  temp_df <- dataframe[,!colnames(dataframe) %in% elim_predictors]
  temp_pls <- plsr(BareGround~., data=temp_df, validation="CV",scale=TRUE, segments=list_spatial_block)
  return(temp_pls)
}

pls_vi_exg <- func_eliminate_predictor(df_pls_1,c("exg"))
pls_vi_exg_bsi <- func_eliminate_predictor(df_pls_1,c("exg","bsi"))
pls_vi_exg_bsi_msavi <- func_eliminate_predictor(df_pls_1,c("exg","bsi","msavi"))

func_q2(pls_vi_exg)[2]
func_q2(pls_vi_exg_bsi)[2]
func_q2(pls_vi_exg_bsi_msavi)[2]

#########################################################
df_pls_all <- df_non_growing_train %>% 
  select(BareGround,all_of(predictors_all))  

pls_all <- plsr(BareGround~., data=df_pls_all, validation="CV",scale=TRUE, segments=list_spatial_block) 
summary(pls_all)
pls_all_non_spatial <- plsr(BareGround~., data=df_pls_all, validation="CV",scale=TRUE) 

# check VIP
plsVarSel::VIP(pls_all,2)
plsVarSel::VIP(pls_all_non_spatial,2)

colours <- c(rep("red",length(vi)),
             rep("orange",length(bands)),
             rep("yellow",length(minerals)),
             rep("green",length(soil)),
             rep("blue",length(topography)))
barplot(plsVarSel::VIP(pls_all,2),col=colours)
barplot(plsVarSel::VIP(pls_all_non_spatial,2),col=colours, las=2)


pls_all_topo <- func_eliminate_predictor(df_pls_all,topography)
pls_all_minerals <- func_eliminate_predictor(df_pls_all,minerals)
pls_all_topo_minerals <- func_eliminate_predictor(df_pls_all,c(topography,minerals))
pls_all_soil <- func_eliminate_predictor(df_pls_all,c(soil,minerals))
pls_all_bsi_exg_minerals_topo <- func_eliminate_predictor(df_pls_all,c("exg",
                                                                     "bsi",
                                                                     soil,minerals,topography))
pls_all_less_than_one <- func_eliminate_predictor(df_pls_all,paste(names(plsVarSel::VIP(pls_all,2)<1.0),collapse=","))


func_q2(pls_all)[2]
func_q2(pls_all_non_spatial)[2]

func_q2(pls_all_topo)[2]
func_q2(pls_all_minerals)[2]
func_q2(pls_all_topo_minerals)[2]
func_q2(pls_all_soil)[2]
func_q2(pls_all_bsi_exg_minerals_topo)[2]
func_q2(pls_all_less_than_one)[2]

pls_one <- plsr(BareGround~., data=df_pls_all[,!colnames(df_pls_all) %in% (names(plsVarSel::VIP(pls_all,2))[(plsVarSel::VIP(pls_all,2)<1.0)])],
                validation="CV",scale=TRUE, segments=list_spatial_block)
VIP(pls_one,2)
func_q2(pls_one)[2]

pls_two <- plsr(BareGround~., data=df_pls_all[,!colnames(df_pls_all) %in% c(names(plsVarSel::VIP(pls_all,2))[(plsVarSel::VIP(pls_all,2)<1.0)],
                                                                            names(plsVarSel::VIP(pls_one,2))[(plsVarSel::VIP(pls_one,2)<0.9)])],
                validation="CV",scale=TRUE, segments=list_spatial_block)
VIP(pls_two,2)
func_q2(pls_two)[2]

########################################################
ariable_selection(predictors_1) %>% 
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
pls_5 <- df_non_growing_train %>% 
  pls_variable_selection(predictors_5) %>% 
  pls_best_model(df_non_growing_train, c(predictors_5))
# get_q2(df_non_growing_train, pls_5)

# validation against test data
ncomp <- selectNcomp(pls_5, method = "onesigma")
pls_prediction_5 <- predict(pls_5, df_non_growing_test[,predictors_5], ncomp=ncomp)[,,1]

######### vegetation indicies + bands + minerals + topography #########
pls_8 <- df_non_growing_train %>% 
  pls_variable_selection(predictors_8) %>% 
  pls_best_model(df_non_growing_train, c(predictors_8))
# get_q2(df_non_growing_train, pls_8)

# validation against test data
ncomp <- selectNcomp(pls_8, method = "onesigma")
pls_prediction_8 <- predict(pls_8, df_non_growing_test[,predictors_8], ncomp=ncomp)[,,1]

######### comparison #########
# obs vs pred plot
pls_pred_list <- list(pls_prediction_1,pls_prediction_2,pls_prediction_5,pls_prediction_8) 

par(mfrow = c(2,2))
for(pred in pls_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}

# a table for RMSE and MAE
df_pls_rmse <- get_rmse_df(pls_pred_list, df_non_growing_test$BareGround)


######### Run them all at once #############
pls_function <- function(predictors){
  temp_model <- df_non_growing_train %>% 
    pls_variable_selection(predictors) %>% 
    pls_best_model(df_non_growing_train, c(predictors))
  return(temp_model)
}
pls_prediction <- function(pls_model, predictors){
  temp_ncomp <- selectNcomp(pls_model, method="onesigma")
  temp_pred <- predict(pls_model, df_non_growing_test[,predictors], ncomp=temp_ncomp)[,,1]
  return(temp_pred)
}

predictors_list <- list(predictors_1,predictors_2,predictors_3,predictors_4,
                        predictors_5,predictors_6,predictors_7,predictors_8,
                        predictors_9,predictors_10,predictors_11,predictors_12,
                        predictors_13)
pls_model_list <- list()
pls_pred_list <- list()
for(predictors in predictors_list){
  temp <- pls_function(predictors) %>% 
    list()
  pls_model_list <- c(pls_model_list,temp)
}

for(i in 1:length(predictors_list)){
  temp_pred <- pls_prediction(pls_model_list[[i]],predictors_list[[i]])
  pls_pred_list <- c(pls_pred_list, list(temp_pred))
}


par(mfrow = c(4,4))
for(pred in pls_pred_list){
  plot_obs_pred(df_non_growing_test$BareGround,pred)
}

df_pls_rmse <- get_rmse_df(pls_pred_list, df_non_growing_test$BareGround)
df_pls_rmse
######### Export ############

for(i in 1:length(pls_model_list)){
  temp_filename <- paste0("pls_",i)
  saveRDS(pls_model_list [[i]], file = paste0("./Script/models/", temp_filename, ".rds"))
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

png("./Script/figures/pls_horizontal.png", width = 1000, height = 280)
  par(mfrow = c(1,4))
  for(pred in pls_pred_list){
    plot_obs_pred(df_non_growing_test$BareGround,pred)
  }
dev.off()  


## Details of the Models 
# VIP
vip_list <- list()

for(pls_model in pls_model_list){
  temp <- plsVarSel::VIP(pls_model,2)
  vip_list <- c(vip_list, list(temp))
 }

vip_df <- vip_list %>% 
  purrr::map_df(~as.data.frame(t(.x)), .id = "model") %>% 
  pivot_longer(-model, names_to = "variable", values_to = "VIP") %>% 
  pivot_wider(names_from = "model", values_from = "VIP")

vip_df

# correlation plot
# Function to plot PLS biplot
plot_pls_biplot <- function(pls_model, comps = 1:2, label_scale = 1.2) {
  # Extract scores for the specified components
  S <- scores(pls_model)[, comps, drop = FALSE]
  
  # Calculate correlation between model matrix and scores
  df.cor <- model.matrix(pls_model) %>%
    as.data.frame() %>%
    cor(S) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    select(variable, axis_1 = "Comp 1", axis_2 = "Comp 2")
 
  # Add correlation of the response variable (Y) with the scores
  Y.cor <- cor(pls_model$Yscores[, 1], S)
  df.cor <- df.cor %>%
    bind_rows(data.frame(variable = "Y", axis_1 = Y.cor[1], axis_2 = Y.cor[2])) 
  # Scale the labels for proper placement
  df.cor <- df.cor %>%
    mutate(label_x = axis_1 * label_scale, label_y = axis_2 * label_scale,
           color = ifelse(variable == "Y", "red", "blue"))

  # Create the plot
  ggplot(df.cor, aes(x = axis_1, y = axis_2)) +
    # Add a unit circle
    geom_circle(data = data.frame(x0 = 0, y0 = 0, r = 1), aes(x0 = 0, y0 = 0, r = 1), inherit.aes = FALSE, color = "gray50") +
    # Add arrows representing the components
    geom_segment(aes(x = 0, y = 0, xend = axis_1, yend = axis_2, color = color),
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 1.5, alpha = 0.5, show.legend = FALSE) +
    # Add labels for the variables at the extended positions
    geom_text(aes(x = label_x, y = label_y, label = variable, color = color), size = 4, show.legend = FALSE) +
    # Add horizontal and vertical lines at the origin
    geom_hline(yintercept = 0, color = "gray80") +
    geom_vline(xintercept = 0, color = "gray80") +
    # Set axis limits and aspect ratio
    xlim(-1.2, 1.2) +
    ylim(-1.2, 1.2) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    labs(x = "Axis_1", y = "Axis_2")
}
plot_pls_biplot(pls_1)
