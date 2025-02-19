############ Functions ##############
# variable selection based on vip value
pls_variable_selection <- function(dataframe, predictors){
  df_temp <- dataframe[,c('BareGround', predictors)]
  pls_initial <- plsr(BareGround ~., data = df_temp, validation="CV", scale = TRUE)
  low_vip_names <-VIP(pls_initial,2) %>% 
    sort() %>% 
    names()
  
  out_df <- data.frame()
  less_than_median <- VIP(pls_initial,2) < median(VIP(pls_initial,2))
  n_loop <- less_than_median %>% 
    unname() %>% 
    as.integer() %>% 
    sum()
  
  all_variables_model <- plsreg1(df_temp[,predictors],df_temp[,'BareGround'], comps = NULL,  crosval = TRUE)
  out_df <- c(all_variables_model$Q2[2,5], sum(all_variables_model$R2), "all_variables") %>% 
    t() %>% 
    as.data.frame()
  colnames(out_df) <- c("Q2cum", "R2cum", "excluded")
  
  for(i in 1:n_loop){
    excluded_name <- low_vip_names[1]
    low_vip_names <- low_vip_names[-1]
    plsBandsVi <- plsreg1(df_temp[,c(low_vip_names)],df_temp[,'BareGround'], comps = NULL, crosval = TRUE)
    temp_Q2cum <- plsBandsVi$Q2[2,5]
    temp_R2sum <- sum(plsBandsVi$R2)
    output <- c(temp_Q2cum, temp_R2sum, excluded_name)
    out_df <- rbind(out_df, output)
  }
  return(out_df)
}
pls_best_model <- function(pls_variable_selection_output, dataframe, predictor){
  max_index <- which.max(pls_variable_selection_output$Q2cum)
  removed_predictor <- pls_variable_selection_output %>% 
    slice(2:max_index) %>% 
    pull(excluded)
  selected_predictor <- predictor[!predictor %in% removed_predictor]
  best_model <- plsr(BareGround ~., data = dataframe[,c('BareGround', selected_predictor)], validation="CV", scale = TRUE)
  return(best_model)
}
# extract q2 using plsreg1
get_q2 <- function(dataframe, pls_model){
  temp_model <- plsreg1(dataframe[,c(colnames(pls_model$model)[-1])],dataframe[,'BareGround'], crosval = TRUE)
  return(temp_model$Q2)
}
get_coefficients <- function(pls_model){
  ncomp <- selectNcomp(pls_model, method = "onesigma")
  std_coef <- coef(pls_model, type="original", ncomp=ncomp, intercept=TRUE)
  unstd_coef <- c(std_coef[1], std_coef[2:nrow(std_coef)]/pls_model$scale)
  return(unstd_coef)
}
# evaluate prediction accuracy
get_mae_rmse <- function(vec_pred, vec_test){
  rmse <- sqrt(mean((vec_pred - vec_test)^2))
  mae <- mean(abs(vec_pred - vec_test))
  out_vec <- c(mae, rmse)
  names(out_vec) <- c("MAE", "RMSE")
  return(out_vec)
}
# plot obs vs preds scatter plot
plot_obs_pred <- function(test_vec, pred_vec){
  plot(test_vec, pred_vec,
       xlim =c(0,100),ylim=c(-20,120),
       xlab = "Photo interpretation (%)", ylab = "Model Prediction (%)" , col=NULL, bg=rgb(0, 0, 1, alpha=0.2), pch=21, cex=1.5
  )
  temp_lm <-lm(pred_vec ~ test_vec)
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  r_sq <- round(summary(temp_lm)$r.squared,3)
  mtext(
    bquote("  R"^"2"~"="~.(r_sq)),
    adj=0,
    line=-2
  )  
}    

# summarise MAE and RMSE as a dataframe
get_rmse_df <- function(pred_list){
  rmse_df <- data.frame(MAE=as.numeric(), RMSE=as.numeric())
  for(pred in pred_list){
    temp_rmse <- get_mae_rmse(pred, df_non_growing_test$BareGround) %>% 
      t() %>% 
      as.data.frame()
    rmse_df <- rbind(rmse_df, temp_rmse)
  }
  return(rmse_df)
}

######### Functions for pls_beta_regression ##################
get_plsbeta_df <- function(dataframe, predictors){
  df_temp <- dataframe %>% 
    select(BareGround, all_of(predictors)) %>% 
    mutate(BareGround = BareGround/100) %>% 
    mutate(BareGround = if_else(BareGround == 0, 1e-8, BareGround)) %>% 
    mutate(BareGround = if_else(BareGround == 1, 1-(1e-8), BareGround))
  return(df_temp) 
}
get_selected_var <- function(ci){
  selected_predictors <- (ci[,7]<0&ci[,8]<0)|(ci[,7]>0&ci[,8]>0)
  temp_pred <- selected_predictors[selected_predictors == TRUE] %>% 
    names()
  return(temp_pred)
}

######### Functions for randomforest_regression ##############
rfe_tune <- function(dataframe, predictors){
  temp_df <- dataframe %>% 
    select(BareGround, all_of(predictors))  
  temp_rfe_tune <- caret::rfe(x = temp_df[,names(temp_df) %in% predictors],
                              y = temp_df$BareGround,
                              sizes = 1:length(predictors),
                              metric = "RMSE",
                              rfeControl = rfeControl(functions =rfFuncs,
                                                      method = "cv",
                                                      repeats = 10, 
                                                      seeds = set.seed(123))
  )
  return (temp_rfe_tune)
}
run_rf <- function(dataframe, predictors, tuned_rfe){
  temp_df <- dataframe %>% 
    select(BareGround, all_of(predictors(tuned_rfe)))
  best <- tuneRF(temp_df[,-1],temp_df$BareGround,doBest=T)
  temp_rf_model <- randomForest(BareGround~., data = temp_df, mtry=best$mtry, ntree=500)
  return(temp_rf_model)
}