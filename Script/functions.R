######### Functions for pls regression with caret #############
# ------------------------------------------------------------------------------
# Function to train a Partial Least Squares (PLS) regression model using caret
# Input:
#   - dataframe: A data frame containing the response and explanatory variables.
#                The response variable must be named "BareGround".
#   - predictors: A character vector of explanatory variable names to use in the model.
#   - blocks: A list of spatial resampling indices, typically from CAST::knndm or similar.
# Output:
#   - A trained caret model object using PLS regression with cross-validation.
func_caret_pls <- function(dataframe, predictors, blocks, n_comp){
  caret_model <- train(
    x = dataframe[, predictors],
    y = dataframe$BareGround,
    method = "pls",                           # Use Partial Least Squares
    preProcess = c('center', 'scale'),        # Standardize predictors before training
    importance = TRUE,                        # Enable variable importance calculation
    metric = "Qsquared",
    trControl = trainControl(
      method = "cv",                          # Cross-validation
      summaryFunction = func_q2,
      selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
      index = blocks$indx_train,               # Cross-validation folds (training)
      indexOut = blocks$indx_test,             # Cross-validation folds (testing)
      savePredictions = "final"               # Save final predictions
    ),
    tuneGrid = data.frame(ncomp = n_comp)
  )
}

# ------------------------------------------------------------------------------
# Function to determine whether log transformation improves correlation with BareGround
# Input:
#   - sf_object: An sf object with predictor variables and a response variable 'BareGround'.
#   - predictors: A character vector of variable names to evaluate for log transformation.
# Output:
#   - A list with three elements:
#       1. A data frame comparing correlations before and after log transformation.
#       2. A character vector of variables selected for log transformation.
#       3. A transformed sf object with selected variables log-transformed.
# func_log_transform <- function(sf_object, predictors) {
#   temp_df <- sf_object %>% 
#     st_drop_geometry() %>% 
#     select(BareGround, all_of(predictors))
#   
#   coef_before <- cor(temp_df)[predictors, "BareGround", drop=FALSE] %>% 
#     as.data.frame() %>% 
#     rename(before_transform = 1)
#   
#   # Identify predictors with negative values
#   add_values <- temp_df %>%
#     summarise(across(everything(), min)) %>%
#     pivot_longer(cols = everything()) %>% 
#     mutate(add=case_when(
#       value <= 0 ~ abs(value)+0.01,
#       TRUE ~ 0
#     )) %>% 
#     select(-value) %>% 
#     deframe()
#   
#   temp_df_log <- temp_df %>% 
#     mutate(across(all_of(predictors), ~.x + add_values[cur_column()])) %>% 
#     mutate(across(all_of(predictors),log)) 
#   
#   coef_after <- cor(temp_df_log)[predictors, "BareGround", drop=FALSE] %>% 
#     as.data.frame() %>% 
#     rename(log_transform = 1)
#   
#   df_log <- bind_cols(coef_before, coef_after)
#   rownames(df_log) <- predictors
#   
#   # Select variables where log-transform improved correlation
#   transform_vars <- df_log %>% 
#     filter(abs(log_transform) > abs(before_transform)) %>% 
#     rownames()
#   
#   transform_neg_var_names <- intersect(names(add_values[add_values>0.01]), transform_vars)
#   
#   if (length(transform_neg_var_names) > 0) {
#     sf_transformed <- sf_object %>%
#       mutate(across(all_of(transform_neg_var_names), ~ .x + add_values[cur_column()])) %>%
#       mutate(across(all_of(transform_vars), log))
#   } else {
#     sf_transformed <- sf_object %>%
#       mutate(across(all_of(transform_vars), log))
#   }
#   
#   return(list(
#     correlation_table=df_log, 
#     transformed_variables=transform_vars, 
#     sf_transformed=sf_transformed, 
#     added_values=add_values[transform_neg_var_names]))
# }

# ------------------------------------------------------------------------------
# Function to split spatial data into training and validation sets with stratified 
# sampling (5-bins) and then create spatial blocks
# Input:
#   - sf_object: An sf object containing spatial point data.
#   - ratio: Proportion of data to assign to training set (default is 0.7).
#   - k: Number of nearest neighbors used in spatial blocking (default is 5).
#   - seed: Random seed for reproducibility (default is 1234).
# Output:
#   - A list containing:
#       1. train_df: A data frame of the training data (geometry dropped).
#       2. valid_df: A data frame of the validation data (geometry dropped).
#       3. spatial_cv_block: A spatial blocking object for cross-validation.
func_data_split <- function(sf_object, ratio = 0.7, k = 5, seed = 1234) {
  set.seed(seed)
  
  sf_bin <- sf_object %>% 
    mutate(y_bin = cut(sf_object$BareGround,
                       breaks =quantile(sf_object$BareGround, probs = seq(0, 1, 0.2), na.rm = TRUE),
                       include.lowest=TRUE,
                       labels=FALSE))
  train_idx <- createDataPartition(sf_bin$y_bin, p=ratio, list=FALSE)
  
  sf_train <- sf_object[train_idx,] 
  df_train <- st_drop_geometry(sf_train) 
  df_valid <- sf_object[-train_idx,] %>%
    st_drop_geometry() 
  
  spatial_blocks <- knndm(sf_train, study_area, k = k)
  
  return(list(train_df = df_train, valid_df = df_valid, spatial_cv_block = spatial_blocks))
}

# ------------------------------------------------------------------------------
# Function to calculate Q² (predictive R²) from a caret-trained PLS model
# Input:
#   - pls_model: A trained PLS model object returned by caret::train.
# Output:
#   - A numeric value representing the Q² (predictive R²) of the model using cross-validated predictions.
# func_calc_q2 <- function(pls_model) {
#   cv_preds <- pls_model$pred
#   best_ncomp <- pls_model$bestTune$ncomp
#   best_preds <- subset(cv_preds, ncomp == best_ncomp)
#   
#   ss_res <- sum((best_preds$obs - best_preds$pred)^2)
#   ss_tot <- sum((best_preds$obs - mean(best_preds$obs))^2)
#   cv_r2  <- 1 - ss_res / ss_tot
#   return(cv_r2)
# }

func_q2 <- function(data, lev =NULL, model=NULL) {
  ss_res <- sum((data$obs - data$pred)^2)
  ss_tot <- sum((data$obs - mean(data$obs))^2)
  cv_r2  <- 1 - ss_res / ss_tot
  default_stats <- defaultSummary(data, lev, model)
  names(cv_r2) <- "Qsquared"
  return(c(default_stats, cv_r2))
}

# ------------------------------------------------------------------------------
# Function to perform variable selection using PLS and VIP, evaluating model performance
# Input:
#   - list_data: A list containing training and validation data frames, and a spatial blocking object.
#   - predictors: A character vector of predictor variable names.
# Output:
#   - A list containing:
#       1. metrics: A data frame with performance metrics for each VIP-based variable subset.
#       2. vip_full_model: A named numeric vector of VIP scores from the full model.

func_pls_var_selection <- function(list_data, predictors, n_comp) {
  full_model <- func_caret_pls(list_data$train_df, predictors, list_data$spatial_cv_block, n_comp)
  vip_scores <- VIP(full_model$finalModel, unlist(full_model$bestTune))
  vip_ordered <- sort(vip_scores) %>% names()
  
  df_metrics <- full_model$results %>%
    mutate(variables = "All")
  
  for (i in 1:(length(vip_ordered) - 2)) {
    sub_model <- func_caret_pls(
      list_data$train_df,
      vip_ordered[-c(1:i)],
      list_data$spatial_cv_block,
      n_comp
    )
    
    metrics <- sub_model$results %>%
      mutate(
        variables = paste(sub_model$finalModel$xNames, collapse = ",")
      )
    
    df_metrics <- bind_rows(df_metrics, metrics)
  }
  
  return(list(metrics = df_metrics, vip_full_model = vip_scores))
}

# ------------------------------------------------------------------------------
# Function to extract unstandardized regression coefficients from a caret PLS model
# Input:
#   - caret_model: A trained PLS model object returned by caret::train.
# Output:
#   - A named numeric vector of unstandardized coefficients, including the intercept.
get_caret_coefficients <- function(caret_model){
  # Extract standardized regression coefficients from the model
  beta_scaled <- coef(caret_model$finalModel, caret_model$bestTune$ncomp)[,,1]
  
  # Retrieve predictor names and standardization parameters from the caret model
  x_names <- names(beta_scaled)
  x_means <- caret_model$preProcess$mean[x_names]  # Predictor means
  x_sds   <- caret_model$preProcess$std[x_names]   # Predictor standard deviations
  y_mean  <- caret_model$finalModel$Ymeans         # Mean of the response variable
  
  # Convert standardized coefficients to unstandardized scale
  beta_unscaled <- beta_scaled / x_sds
  
  # Compute unstandardized intercept
  intercept_unscaled <- y_mean - sum(beta_unscaled * x_means)
  
  # Combine intercept and unstandardized coefficients into a named vector
  coef_caret_unscaled <- c(Intercept = intercept_unscaled, beta_unscaled)
  
  return(coef_caret_unscaled)
}


# ######### Functions for pls_beta_regression ##################
# get_plsbeta_df <- function(dataframe, predictors){
#   df_temp <- dataframe %>% 
#     select(BareGround, all_of(predictors)) %>% 
#     mutate(BareGround = BareGround/100) %>% 
#     mutate(BareGround = if_else(BareGround == 0, 1e-8, BareGround)) %>% 
#     mutate(BareGround = if_else(BareGround == 1, 1-(1e-8), BareGround))
#   return(df_temp) 
# }
# get_selected_var <- function(ci){
#   selected_predictors <- (ci[,7]<0&ci[,8]<0)|(ci[,7]>0&ci[,8]>0)
#   temp_pred <- selected_predictors[selected_predictors == TRUE] %>% 
#     names()
#   return(temp_pred)
# }
# predict_plsRbeta <- function(plsRbeta_model, dataframe){
#   temp_df <- dataframe %>% 
#     mutate(intercept = 1) %>% 
#     select(intercept, all_of(rownames(plsRbeta_model$Coeffs)[-1]
#     ))
#   
#   # Compute predictions using the PLS-beta model
#   temp_logit_values <- as.matrix(temp_df) %*% as.vector(plsRbeta_model$Coeffs)
#   temp_prediction <- exp(temp_logit_values) / (1 + exp(temp_logit_values)) * 100
#   return(temp_prediction)
# }

# ######### Functions for randomforest_regression ##############
# # ------------------------------------------------------------------------------
# # Function to select predictors based on Variance Inflation Factor (VIF)
# # Input: 
# #   - dataframe: A data frame containing the response and explanatory variables. 
# #                The response variable must be named "BareGround".
# #   - predictors: A character vector specifying the names of candidate explanatory variables.
# # Output:
# #   - A character vector of selected predictors with VIF less than the specified threshold.
# func_colinear <- function(dataframe, predictors) {
#   pref_order <- dataframe %>%
#     select(BareGround,all_of(vi)) %>% 
#     cor() %>% 
#     as.data.frame() %>% 
#     select(BareGround) %>% 
#     arrange(desc(abs(BareGround))) %>% 
#     rownames()
#   cor_select(
#     df = dataframe,
#     predictors = predictors,
#     preference_order = pref_order[2:6], # prioritise highly correlated variables
#     max_cor = 0.90,
#     quiet = TRUE
#   ) %>%
#     as.vector()
# }
# 
# # ------------------------------------------------------------------------------
# # Function to perform Recursive Feature Elimination (RFE) using Random Forest
# # Input: 
# #   - dataframe: A data frame containing the response and explanatory variables. 
# #                The response variable must be named "BareGround".
# #   - predictors: A character vector specifying the names of candidate explanatory variables.
# #   - blocks: A list of spatial resampling indices created with CAST::knndm (or similar).
# #   - seed: An integer to set the random seed for reproducibility.
# # Output:
# #   - A list containing:
# #     [[1]]: A data frame of model performance metrics (e.g., RMSE, R²) for each subset size.
# #     [[2]]: The number of predictors in the most parsimonious model within one standard error (1-SE) of the best model.
# #     [[3]]: A character vector of selected predictors from the 1-SE model.
# func_rfe <- function(dataframe, predictors, blocks, seed = 1234) {
#   set.seed(seed)  # Ensure reproducibility
#   
#   output <- rfe(
#     x = dataframe[, predictors],
#     y = dataframe$BareGround,
#     sizes = 1:length(predictors),
#     rfeControl = rfeControl(
#       functions = rfFuncs,
#       method = "cv",
#       index = blocks$indx_train,
#       indexOut = blocks$indx_test,
#     )
#   )
#   
#   nparam <- oneSE(output$results, metric = "RMSE", maximize = FALSE, num=5)
#   selected_pred <- output$optVariables[1:nparam]
#   
#   return(list(
#     results = output$results,
#     n_selected = nparam,
#     selected_predictors = selected_pred
#   ))
# }
# 
# # ------------------------------------------------------------------------------
# # Function to train a Random Forest model using the caret package
# # Input: 
# #   - dataframe: A data frame containing the response and explanatory variables. 
# #                The response variable must be named "BareGround".
# #   - predictors: A character vector specifying the names of explanatory variables to use.
# #   - blocks: A list of spatial resampling indices created with CAST::knndm (or similar).
# #   - seed: An integer to set the random seed for reproducibility.
# # Output:
# #   - A caret model object trained using cross-validation.
# func_caret_rf <- function(dataframe, predictors, blocks, seed = 1234) {
#   set.seed(seed)  # Ensure reproducibility
#   
#   caret_model <- train(
#     x = dataframe[, predictors],
#     y = dataframe$BareGround,
#     method = "rf",
#     preProcess = c("center", "scale"),
#     importance = TRUE,
#     trControl = trainControl(
#       method = "cv",
#       index = blocks$indx_train,
#       indexOut = blocks$indx_test,
#       savePredictions = "final"
#     )
#   )
#   
#   return(caret_model)
# }
# 
# 


########## Functions for validation #############
# plot obs vs preds scatter plot
plot_obs_pred <- function(test_vec, pred_vec){
  plot(test_vec, pred_vec,
       xlim =c(0,100),ylim=c(-20,120),
       xlab = "Photo interpretation (%)", ylab = "Model Prediction (%)" , col=NULL, bg=rgb(0, 0, 1, alpha=0.3), pch=21, cex=1.5
  )
  temp_lm <-lm(pred_vec ~ test_vec)
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  r_sq <- round(summary(temp_lm)$r.squared,3)
  rmse <- mean((test_vec-pred_vec)^2) %>% 
    sqrt() %>% 
    round(3)
  mae <- mean(abs(test_vec-pred_vec)) %>% 
    round(3)
  mtext(
    bquote("R"^"2"~"="~.(r_sq)),
    adj=0.05,
    line=-4.6
  )
  mtext(
    bquote("RMSE = "~.(rmse)),
    adj=0.05,
    line=-2
  )
  mtext(
    bquote("MAE = "~.(mae)),
    adj=0.05,
    line=-3.3
  )
}

plot_slope_flat <- function(dataframe, valid_vec){
  df_temp <- dataframe %>% 
    select(BareGround,Notes) %>% 
    mutate(pred = valid_vec) 
  df_flat_temp <- df_temp %>% 
    filter(!str_detect(Notes, "_Slope"))
  df_slope_temp <- df_temp %>% 
    filter(str_detect(Notes, "_Slope"))
  # plot flat points with blue circles
  plot(df_flat_temp$BareGround, df_flat_temp$pred,
       xlim =c(0,100),ylim=c(-20,120),
       xlab = "Transect Measurement (%)", ylab = "Model Prediction (%)" , col=NULL, bg=rgb(0, 0, 1, alpha=0.4), pch=21, cex=1.5
  )
  # plot slope points with red triangles
  points(df_slope_temp$BareGround, df_slope_temp$pred, col=NULL, bg=rgb(1, 0, 0, alpha=0.4), pch=24, cex=1.5)
  
  temp_lm <-lm(df_temp$pred ~ df_temp$BareGround)
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  r_sq <- round(summary(temp_lm)$r.squared,3)
  rmse <- mean((df_temp$pred-df_temp$BareGround)^2) %>% 
    sqrt() %>% 
    round(3)
  mae <- mean(abs(df_temp$pred-df_temp$BareGround)) %>% 
    round(3)
  mtext(
    bquote("R"^"2"~"="~.(r_sq)),
    adj=0.05,
    line=-4.6
  )
  mtext(
    bquote("RMSE = "~.(rmse)),
    adj=0.05,
    line=-2
  )
  mtext(
    bquote("MAE = "~.(mae)),
    adj=0.05,
    line=-3.3
  )
  legend("bottomright",
         legend = c("Flat", "Slope"), 
         col=c(rgb(0, 0, 1, alpha=0.2), rgb(1, 0, 0, alpha=0.2)),
         pt.bg = c(rgb(0, 0, 1, alpha=0.2), rgb(1, 0, 0, alpha=0.2)), 
         pch = c(21, 24), 
         cex = 1.2,
         bty = "n"
         )
} 
# 4 plots at the same time
plot_4plots_slope_flat <- function(pred_list, dataframe){
  par(mfrow = c(2,2))
  for(pred in pred_list){
    plot_slope_flat(dataframe,pred)
  }
}

# ######### Functions for bootstrapping ##################
# run_boot <- function(vec_pred,vec_test){
#   temp_df <- data.frame(vec_pred,vec_test)
#   colnames(temp_df) <- c("pred", "obs")
#   temp_func <- function(data,i){
#     x <- data[i,"obs"]
#     y <- data[i,"pred"]
#     rmse <- sqrt(mean((x - y)^2))
#     mae <- mean(abs(x - y))
#     temp_lm <-lm(y ~ x)
#     r_sq <- summary(temp_lm)$r.squared
#     out_vec <- c(mae, rmse,r_sq)
#     names(out_vec) <- c("MAE", "RMSE","R^2")
#     return(out_vec)
#   }
#   temp_boot <- boot(data=temp_df, 
#                     statistic=temp_func,
#                     R=1000,
#                     stype="i")
#   temp_mae_ci <- boot.ci(temp_boot, type="bca", index=1) %>% 
#     get_ci(type="bca")
#   temp_rmse_ci <- boot.ci(temp_boot, type="bca", index=2) %>% 
#     get_ci(type="bca")
#   temp_r2_ci <- boot.ci(temp_boot, type="bca", index=3) %>% 
#     get_ci(type="bca")
#   
#   vec_mae <- c(temp_boot$t0[1], mean(temp_boot$t[,1]), temp_mae_ci) %>% 
#     unlist()
#   vec_rmse <- c(temp_boot$t0[2], mean(temp_boot$t[,2]), temp_rmse_ci) %>% 
#     unlist()
#   vec_r2 <- c(temp_boot$t0[3], mean(temp_boot$t[,3]), temp_r2_ci) %>% 
#     unlist()
#   
#   out_df <- data.frame(vec_mae, vec_rmse, vec_r2) %>% 
#     t()
#   colnames(out_df) <- c("Original", "Mean_boot", "bca_Lower","bca_Upper")
#   rownames(out_df) <- c("MAE", "RMSE","R^2")
#   return(out_df)
# }
# # create a dataframe suitable to draw a plot
# get_dotplot_df <- function(pred_list, vec_test){
#   out_df <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
#   colnames(out_df) <- c("Original", "Mean_boot", "bca_Lower","bca_Upper", "model")
#   for(i in 1:length(pls_valid_list)){
#     temp_df <- run_boot(pred_list[[i]], df_transect$BareGround) %>% 
#       as.data.frame() %>% 
#       mutate(model=paste0("model_",i))
#     out_df <- rbind(out_df,temp_df)
#   }
#   return(out_df)
# }
