library(pls)
library(tidyverse)
library(ggplot2)
library(randomForest)
library(boot)
library(fwb)

############ Load CSV file ##############
df_transect <- read.csv("./CSV/transect_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 52 rows

head(df_transect)

# Load a script
source("./Script/functions.R")

# Load models
list_models <- list.files("./Script/models", pattern = "\\.rds$", full.names = TRUE)

read_models <- function(file_name = tmp) {
  f <- readRDS(file_name)
  temp_model_name <- basename(file_name)
  temp_model_name <- gsub("\\..+$", "", temp_model_name)
  assign(temp_model_name, f, envir = .GlobalEnv)
  ls() # shows that my_data is not in the current environment
}

for(model in list_models){
  read_models(model)
}


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

############ Prediction ###################
pls_valid_1 <- predict(pls_1, df_transect[,predictors_1], ncomp=2)[,,1]
pls_valid_2 <- predict(pls_2, df_transect[,predictors_2], ncomp=2)[,,1]
pls_valid_3 <- predict(pls_3, df_transect[,predictors_3], ncomp=2)[,,1]
pls_valid_4 <- predict(pls_4, df_transect[,predictors_4], ncomp=2)[,,1]

plsbeta_valid_1 <- predict_plsRbeta(plsbeta_1, df_transect)
plsbeta_valid_2 <- predict_plsRbeta(plsbeta_2, df_transect)
plsbeta_valid_3 <- predict_plsRbeta(plsbeta_3, df_transect)
plsbeta_valid_4 <- predict_plsRbeta(plsbeta_4, df_transect)

rf_valid_1 <- predict(rf_1, df_transect[,predictors_1])
rf_valid_2 <- predict(rf_2, df_transect[,predictors_2])
rf_valid_3 <- predict(rf_3, df_transect[,predictors_3])
rf_valid_4 <- predict(rf_4, df_transect[,predictors_4])


######### comparison #########
# obs vs pred plot
pls_valid_list <- list(pls_valid_1,pls_valid_2, pls_valid_3, pls_valid_4) 
plsbeta_valid_list <- list(plsbeta_valid_1,plsbeta_valid_2, plsbeta_valid_3, plsbeta_valid_4) 
rf_valid_list <- list(rf_valid_1,rf_valid_2,rf_valid_3,rf_valid_4) 

plot_4plots_slope_flat(pls_valid_list, df_transect)
plot_4plots_slope_flat(plsbeta_valid_list, df_transect)
plot_4plots_slope_flat(rf_valid_list, df_transect)

# a table for RMSE and MAE
df_pls_rmse <- get_rmse_df(pls_valid_list, df_transect$BareGround)
df_plsbeta_rmse <- get_rmse_df(plsbeta_valid_list, df_transect$BareGround)
df_rf_rmse <- get_rmse_df(rf_valid_list, df_transect$BareGround)

df_valid_rmse <- rbind(df_pls_rmse, df_plsbeta_rmse, df_rf_rmse) %>% 
  mutate(model=paste0(c(rep("pls_", 4),rep("plsbeta_", 4),rep("rf_",4)),rep(1:4,3))) #"pls_1", "pls_2",...,"rf_4" 
# df_valid_rmse

######### Export ############
# rmse df
write.csv(df_valid_rmse, "./Script/csv/validation_rmse.csv")

# plot
png("./Script/figures/pls_valid.png", width = 600, height = 700)
plot_4plots_slope_flat(pls_valid_list, df_transect)
dev.off()  

png("./Script/figures/plsbeta_valid.png", width = 600, height = 700)
plot_4plots_slope_flat(plsbeta_valid_list, df_transect)
dev.off()  

png("./Script/figures/rf_valid.png", width = 600, height = 700)
plot_4plots_slope_flat(rf_valid_list, df_transect)
dev.off()  


par(mfrow = c(3,4))
plot_list <- c(pls_valid_list, plsbeta_valid_list, rf_valid_list)
for(pred in plot_list){
    plot_slope_flat(df_transect, pred)
  }


######### Bootstrapping ############

pls_boot <- get_dotplot_df(pls_valid_list, df_transect$BareGround)
plsbeta_boot <- get_dotplot_df(plsbeta_valid_list, df_transect$BareGround)
rf_boot <- get_dotplot_df(rf_valid_list, df_transect$BareGround)

df_dotplot <- rbind(pls_boot, plsbeta_boot, rf_boot) %>% 
  mutate(metric=rep(c("MAE","RMSE", "R^2"),4*3)) %>% 
  mutate(model_n=rep(c("VI","+bands","+EMIT","topo"), 3, each=3)) %>% 
  mutate(algo=rep(c("Pls","Pls-beta","Random Forest"),each=3*4)) %>% 
  mutate(model=paste(algo, model_n, sep="_")) %>% 
  mutate(model=as.factor(model))
df_dotplot$model<-factor(df_dotplot$model,
                         levels=c(paste(rep(c("Random Forest","Pls-beta","Pls"),each=4),
                                        rep(c("topo","+EMIT","+bands","VI"), 3),
                                        sep="_")

                                                                    ))
#function for ggplot2
get_dotplot <- function(met){
  if(met=="R^2"){
    temp_hline <- df_dotplot %>% 
      filter(metric==met) %>% 
      select(Original) %>% 
      max()
    ymax<- 1
  }else{
    temp_hline <- df_dotplot %>% 
      filter(metric==met) %>% 
      select(Original) %>% 
      min()
    ymax <- 22
  }

  df_dotplot %>% 
    filter(metric ==met) %>% 
    ggplot(aes(x = model, y = Original, ymin = bca_Lower, ymax = bca_Upper, colour = algo)) +
    geom_vline(xintercept = 1:12, linetype = "dotted")+
    geom_hline(yintercept = temp_hline, linetype = "dotted")+
    geom_errorbar(position = position_dodge(width = 1), width = 0, linewidth=1.5) +
    geom_point(position = position_dodge(width = 1),size=3) +
    geom_point(aes(x = model, y = Mean_boot), position = position_dodge(width = 1),pch=1, size=3) +
    scale_y_continuous(limits=c(0,ymax))+
    coord_flip() +
    ylab(met)+
    theme_bw()
}

get_dotplot("MAE")
get_dotplot("RMSE")
get_dotplot("R^2")

