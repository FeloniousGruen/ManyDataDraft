
source("~/Simulation Function.R")
source("Functions.R")
# ----------------------------------------------------------------------------
# Parallel implementation of simulations for different values of bias
# ----------------------------------------------------------------------------
# Create empty list to store summary dataframes
result_dfs <- list()
bias_list <- c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,2.25,2.5,3,4,5,7.5,10,15,20)


n_boot <- 100



for (bias in bias_list) {
  
  # Run sim_func once to determine structure (using i=1 as a test)
  # If sim_func can fail on the first call, you might want to wrap this in a tryCatch as well,
  # or choose a different representative i. This assumes it won't fail here.
  test_result <- sim_func(bias, 1, confounding_fn_inc = T, atmle_inc = T,procova_inc = T,
                          RBesT_inc = F)
  col_names <- colnames(test_result)
  n_cols <- length(col_names)
  
  # Setup parallel backend
  num_cores <- detectCores() - 1  # Leave two cores for other tasks (1 is often too few)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Run simulations in parallel for the current bias
  result <- foreach(i = 1:n_boot,
                    .combine = 'rbind', 
                    .packages = c("data.table", "causl", "ManyData", "mvtnorm", "survey", "loo",
                                  "MASS", "magrittr", "numDeriv", "dplyr", "ggplot2", "matrixcalc",
                                  "IntegrativeHTEcf", "tidyr", "RBesT", "sl3", "atmle", "SuperLearner", "purrr",
                                  "MatchIt")) %dopar% {
                                    tryCatch(
                                      {
                                        sim_func(bias, i, confounding_fn_inc = T, atmle_inc = T,procova_inc = T,
                                                 RBesT_inc = F)
                                      },
                                      error = function(e) {
                                        # Return a zero-filled row with the same column structure
                                        zero_df <- data.frame(matrix(0, nrow = 1, ncol = n_cols))
                                        colnames(zero_df) <- col_names
                                        zero_df
                                      }
                                    )
                                    
                                  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Store results and save to CSV
  result_dfs[[as.character(bias)]] <- result
  run_name <- paste("Bias", bias)
  write.csv(result, paste(run_name, "Save_name.csv"))
}



true_ATE <- 0.5
result_dfs <- lapply(result_dfs, function(df) {
  
  
  # Remove rows where all values are zero
  df <- df[rowSums(df != 0) > 0, ]
  
  return(df)
})

compute_mse_opt <- function(df) {
  df$MSE_new_opt <- (df$new_opt_ATE_est - true_ATE)^2 
  return(df)
}
compute_mse_new<- function(df) {
  df$MSE_new <- (df$new_ATE_est - true_ATE)^2 
  return(df)
}
compute_mse_old <- function(df) {
  df$MSE_old <- (df$old_ATE_est - true_ATE)^2 
  return(df)
}

compute_mse_meta <- function(df) {
  df$MSE_simple_meta <- (df$simple_meta_ATE_est - true_ATE)^2 
  return(df)
}
compute_mse_rct <- function(df) {
  df$MSE_rct <- (df$rct_ATE_est - true_ATE)^2 
  return(df)
}
compute_mse_cf <- function(df) {
  df$MSE_cf <- (df$cf_ATE - true_ATE)^2 
  return(df)
}
compute_mse_meta_prior <- function(df) {
  df$MSE_meta_prior <- (df$posterior_mean - true_ATE)^2 
  return(df)
}

compute_mse_atmle <- function(df) {
  df$MSE_atmle <- (df$atmle_ATE - true_ATE)^2 
  return(df)
}
compute_mse_procova <- function(df) {
  df$MSE_procova <- (df$procova_ATE - true_ATE)^2 
  return(df)
}
compute_mse_oberst <- function(df) {
  df$MSE_oberst <- (df$oberst_ATE_est - true_ATE)^2 
  return(df)
} 
compute_mse_shrink <- function(df) {
  df$MSE_shrink_gsar <- (df$green_strawd_shrink_ATE_est - true_ATE)^2 
  df$MSE_shrink_rosenman <- (df$rosenman_shrink_ATE_est - true_ATE)^2 
  return(df)
} 
compute_mse_aipw <- function(df) {
  df$MSE_aipw <- (df$AIPW_ATE_est - true_ATE)^2 
  return(df)
}


result_dfs <- lapply(result_dfs, compute_mse_old)
result_dfs <- lapply(result_dfs, compute_mse_opt)
result_dfs <- lapply(result_dfs, compute_mse_meta)
result_dfs <- lapply(result_dfs, compute_mse_rct)
result_dfs <- lapply(result_dfs, compute_mse_cf)
result_dfs <- lapply(result_dfs, compute_mse_meta_prior)
result_dfs <- lapply(result_dfs, compute_mse_atmle)
result_dfs <- lapply(result_dfs, compute_mse_procova)
result_dfs <- lapply(result_dfs, compute_mse_oberst)
result_dfs <- lapply(result_dfs, compute_mse_shrink)
result_dfs <- lapply(result_dfs, compute_mse_aipw)









