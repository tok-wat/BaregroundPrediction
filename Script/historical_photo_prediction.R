library(pls)
library(tidyverse)
library(ggplot2)

############ Load CSV file ##############
df_historical <- read.csv("./CSV/historical_photo_extracted.csv") %>% 
  select(-c(system.index, .geo, vegType)) %>%
  filter(Month > 4) %>% 
  na.omit()  # 31 rows

head(df_historical)
hist(df_historical$Year,breaks=1989:2013,
     xlab="Year",
     xlim=c(1989,2013),
     main=NULL)

# Load a script# LoYearad a script
source("./Script/functions.R")

# Load the model
pls_1 <- readRDS("./Script/models/pls_1.rds")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
predictors_1 <- c(vi)

# Initial model prediction 
x_pred <- predict(pls_1, df_historical[,predictors_1], ncomp=2)[,,1]/100

data.frame(x_pred, df_historical$BareGround) %>% 
  head()

get_mae_rmse(x_pred*100, df_historical$BareGround)
plot_slope_flat(df_historical, x_pred*100)
mtext(adj=0, line=-1,padj=3, paste0("  MAE: ", round(get_mae_rmse(x_pred*100, df_historical$BareGround)[1],2)))
mtext(adj=0, line=-2,padj=3, paste0("  RMSE: ", round(get_mae_rmse(x_pred*100, df_historical$BareGround)[2],2)))

# residuals or prediction error
data.frame(x_pred, df_historical$BareGround) %>% 
  mutate(pred_error = x_pred*100-x_pred-df_historical.BareGround) %>% 
  mutate(Year = df_historical$Year) %>%
  mutate(Month = df_historical$Month) %>% 
  select(pred_error, Year, Month) %>% 
  plot()

# extracting coefficients 
extract_coefficients <- function(plsmodel,ncomp){
  std_coef <- coef(plsmodel, type = "original", ncomp = ncomp, intercept = TRUE)
  unstd_coef <- c(std_coef[1], std_coef[2:nrow(std_coef)]/plsmodel$scale)
  return(unstd_coef)
}
extract_coefficients(pls_1, 2)  
