##########################################################
# Use of multiple datasets for LEADER ATE estimation
##########################################################

##########################################################
# Package Installation
##########################################################
cran_packages <- c(
  "data.table","causl","mvtnorm","survey","doParallel","loo","coda",
  "parallel","foreach","doSNOW","MASS", "magrittr",
  "numDeriv","dplyr","ggplot2","matrixcalc","tidyr", "RBesT", "SuperLearner", 
  "MatchIt", "purrr", "remotes"
)

# Install CRAN packages if missing
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing CRAN package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
}
# Handle GitHub packages
# remotes::install_github("XiLinStats/ManyData")
# remotes::install_github("IntegrativeStats/IntegrativeHTEcf")
# remotes::install_github("tlverse/sl3")
# remotes::install_github("tq21/atmle")


# Load all packages
# - CRAN packages
lapply(cran_packages, library, character.only = TRUE)
# - GitHub packages
library(ManyData)            
library(IntegrativeHTEcf)     
library(atmle)
library(sl3)


# Load required functions
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
  result_dfs <- add_all_mse(result_dfs, true_ATE = 0.5)
  write.csv(result, paste(run_name, "Save_name.csv"))
}









