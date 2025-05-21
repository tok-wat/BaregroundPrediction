library(betareg)
library(pls)
library(tidyverse)
library(ggplot2)

############ Load CSV file ##############
df_transect <- read.csv("./CSV/transect_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # 52 rows

head(df_transect)

# Load a script
source("./Script/functions.R")

# Load the model
pls_1 <- readRDS("./Script/models/pls_1.rds")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
predictors_1 <- c(vi)

# Initial model prediction 
x_pred <- predict(pls_1, df_transect[,predictors_1], ncomp=2)[,,1]/100

# Response variable - field measured fraction of bare ground
y_transect <- df_transect %>% 
  mutate(BareGround = BareGround/100) %>% 
  mutate(BareGround = if_else(BareGround == 0, 1e-8, BareGround)) %>% 
  mutate(BareGround = if_else(BareGround == 1, 1-(1e-8), BareGround)) %>% 
  select(BareGround) %>% 
  unlist()

# beta regression model
# slope を入れるかどうか検討
beta_model <- betareg(y_transect ~ x_pred, link = "logit")
y_calibrated <- predict(beta_model, newdata = data.frame(y_pred = x_pred), type = "response")

data.frame(x_pred, y_transect, calibrated) %>% 
  head()

get_mae_rmse(calibrated*100, y_transect*100)
plot_slope_flat(df_transect, y_calibrated*100)
mtext(adj=0, line=-1,padj=3, paste0("  MAE: ", round(get_mae_rmse(calibrated*100, y_transect*100)[1],2)))
mtext(adj=0, line=-2,padj=3, paste0("  RMSE: ", round(get_mae_rmse(calibrated*100, y_transect*100)[2],2)))

get_mae_rmse(x_pred*100, y_transect*100)
plot_slope_flat(df_transect, x_pred*100)
mtext(adj=0, line=-1,padj=3, paste0("  MAE: ", round(get_mae_rmse(x_pred*100, y_transect*100)[1],2)))
mtext(adj=0, line=-2,padj=3, paste0("  RMSE: ", round(get_mae_rmse(x_pred*100, y_transect*100)[2],2)))
