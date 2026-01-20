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
#    (a) ManyData from XiLinStats
#    (b) IntegrativeHTEcf from IntegrativeStats

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



##########################################################
# Required specifications
##########################################################
#dataframe <- (data - with covariates we discussed, outcome coded as "Y", treatment as "ARM", and region as region) 
#region_name <- (this should refer to prefered region)


# Parameters
source("~/Functions.R")


# HERE ####
sim_func <- function(omega,seed, RBesT_inc = T, atmle_inc = T, 
                     procova_inc = T, confounding_fn_inc = T,
                     n_o = 500, n_e = 500, n_datasets = 10){
  
  ########################################################################################################
  # Simulation
  ########################################################################################################


    set.seed(seed+1)
  

    prop_list <- seq_len(n_datasets) * omega / n_datasets
  
    data.list <- generate.data.list.HTE(n_e,n_o,seed,coeff.U.list = prop_list, X.diff = T, 
                                        U.diff = F, U1.mean.shift = 1, n_datasets = n_datasets,
                                        HTE = T)
  
    data_exp_b <- data.list[[1]][,c(1:6)]                  # randomised arm (S = 1)
    
    for (i in seq_len(n_datasets)) {              # observational arms (S = 2,3,…)
      assign(paste0("data_obs_b", i), data.list[[i + 1]][,c(1:6)])
    }
    
    ## Combine every element except the first (that one is dat.r)
    data_obs_b <- rbindlist(data.list[-1], use.names = TRUE, fill = TRUE)[,c(1:6)]
    
    
    
    ########################################################################################################
    # Fit_causl estimation
    ########################################################################################################
    # Update families to remove hidden confounder U for estimation
    
    family <- list(1, c(5, 5,5), 1, 1)        # used for *all* arms
    
    ## ----- 2.  Build the single randomised arm ----------------------------- ##
    forms_obs <- forms_exp <- list(
      c(X1 ~ 1),
      list(A ~ X1 + X2 , X3 ~ 1, X2 ~ 1),                # treatment randomised
      Y ~ A + X3 + A:X3 + X2,
      ~ 1
    )
    
    forms_exp[[2]] <- list(A ~ 1, X3 ~ 1, X2 ~ 1)
    
    
    forms2 <- list(obs = causl:::tidy_formulas(unlist(forms_obs[-2]), kwd = "cop"),
                   exp = causl:::tidy_formulas(unlist(forms_obs[-2]), kwd = "cop"))
    
    
    full_form <- list(obs =  causl:::merge_formulas(forms2$obs),
                      exp =  causl:::merge_formulas(forms2$exp))
    
    msks <- list(obs = causl:::masks(forms2$obs,family = unlist(family[-2]),full_form$obs$wh),
                 exp = causl:::masks(forms2$exp,family = unlist(family[-2]),full_form$exp$wh))
    
    vars <- causl:::lhs(unlist(forms2$obs[-length(forms2$obs)]))

    # Model with copula fitting
    
    fit_exp <- fit_causl(dat = data_exp_b, formulas = unlist(forms2$exp), family = unlist(family[-2])
                         ,control = list(sandwich = F)
    )
    
    
    fit_obs <- fit_causl(dat = data_obs_b, formulas = unlist(forms2$obs), family = unlist(family[-2])
                         ,control = list(sandwich = F)
    )
    
    
    # List to store fit results
    fit_obs_list <- lapply(1:n_datasets, function(i) {
      # Retrieve the data (data_obs_b1, data_obs_b2, ...)
      current_data <- get(paste0("data_obs_b", i))
      # Fit the causl model
      fit_result <- fit_causl(dat = current_data, 
                              formulas = unlist(forms2$obs), 
                              family = unlist(family[-2])
                               , control = list(sandwich = F) # Uncomment if needed
      )
      
      return(fit_result)
    })
    
    # Optionally, assign each fit result to variables fit_obs1, fit_obs2, ..., fit_obs10
    for (i in 1:n_datasets) {
      assign(paste0("fit_obs", i), fit_obs_list[[i]])
    }
    
    
    mm <- list(obs = fit_obs$mm, exp = fit_exp$mm)

    
  
  ########################################################################################################
  # Estimating E(A) to implement new estimation method
  ########################################################################################################
  A_empirical <- -fit_exp$FI
  
  ########################################################################################################
  # Using Matrices A, V and W to calculate ELPD using one observational and one RCT
  ########################################################################################################

  res_single_dataset <- new_method_optim(fit_exp,list(fit_obs),  A_empirical)
  eta_max <- res_single_dataset$par
  
  # Get estimate
  new_samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = eta_max, n_sample = 10000,
                                  use_sandwich = F)
  
  indices_for_A <- find_indices_for_var(fit_exp, "A")
  data <- calculate_summary_flexible(new_samples, indices_for_A, data_exp_b, main_var = "A",
                                     name = "new")
  
  data <- cbind(data, eta_max = eta_max)
  
  
  ########################################################################################################
  # New method using numerical optimisation
  ########################################################################################################

  res <- new_method_optim(fit_exp,fit_obs_list,  A_empirical)
  
  # Get estimate
  new_opt_samples <- approx_posterior_multiple_datasets(fit_exp = fit_exp,  
                                                        fit_obs_list = fit_obs_list, 
                                                        eta_list = res$par, 
                                                        n_sample = 10000,
                                                        use_sandwich = F)
  

  new_method_data <- calculate_summary_flexible(new_opt_samples, indices_for_A, 
                              data_exp_b, main_var = "A", name = "new_opt")
  print(new_method_data)
  data <- cbind(data, new_method_data)
  
  
  
  
  
  ########################################################################################################
  # Old power likelihood method with one random and one observational dataset
  ########################################################################################################
  

  
  old_result <- old_power_likelihood(data_exp_b, fit_exp, fit_obs, mm, msks, family, vars)
  old_eta <- old_result[[1]]
  old_samples <- old_result[[2]]
  
  old_method_data <- calculate_summary_flexible(old_samples, indices_for_A, 
                                   data_exp_b, main_var = "A",name = "old")
  data <- cbind(data, old_method_data)

  
  
  
  ########################################################################################################
  # RCT and simple meta estimate
  ########################################################################################################
  
  rct_samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = 0, n_sample = 10000,
                                  use_sandwich = F)
  
  rct_method_data <- calculate_summary_flexible(rct_samples, indices_for_A, 
                                 data_exp_b, main_var = "A",name = "rct")
  data <- cbind(data, rct_method_data)
  
  
  
  
  
  simple_meta_samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = 1, n_sample = 10000,
                                  use_sandwich = F)
  
  simple_meta_method_data <- calculate_summary_flexible(simple_meta_samples, indices_for_A, 
                                data_exp_b, main_var = "A",name = "simple_meta")
  data <- cbind(data, simple_meta_method_data)
  
  
  obs_samples <- approx_posterior(fit_exp = fit_obs,  fit_exp, eta = 0, n_sample = 10000,
                                  use_sandwich = F)
  
  obs_only <- calculate_summary_flexible(obs_samples, indices_for_A, 
                                                data_exp_b, main_var = "A",name = "obs_only")
  data <- cbind(data, obs_only)
  
  
  obs_samples_most_cnf <- approx_posterior(fit_exp = fit_obs_list[[n_datasets]],  fit_exp, eta = 0, n_sample = 10000,
                                  use_sandwich = F)
  
  obs_only_most_cnf <- calculate_summary_flexible(obs_samples_most_cnf, indices_for_A, 
                                         data_exp_b, main_var = "A",name = "obs_only_most_cnf")
  data <- cbind(data, obs_only_most_cnf)
  
  
  obs_samples_least_cnf <- approx_posterior(fit_exp = fit_obs_list[[1]],  fit_exp, eta = 0, n_sample = 10000,
                                           use_sandwich = F)
  
  obs_only_least_cnf <- calculate_summary_flexible(obs_samples_most_cnf, indices_for_A, 
                                                  data_exp_b, main_var = "A",name = "obs_only_least_cnf")
  data <- cbind(data, obs_only_least_cnf)
  
  

  
  
  ##############################################################################
  # AIPW method
  ##############################################################################
  data_exp_b$S <- 1
  data_obs_b$S <- 0
  confidence.level = 0.05
  reweight.obs = FALSE
  aipw <- aipw.est(rbind(data_exp_b, data_obs_b), covariates = c("X1","X2","X3"),
                       Q.SL.library =  c("SL.mean","SL.speedglm"),
                       g.SL.library =  c("SL.mean","SL.speedglm"),
                       prop.bound = 0.025,
                       estimates = c("obs","RCT"),
                       reweight.obs = reweight.obs,
                       confidence.level = confidence.level)
  
  aipw_data <- data.frame(AIPW_ATE_est = aipw$rct$ate,
                          AIPW_ATE_est_low_CI = aipw$rct$ate.l,
                          AIPW_ATE_est_high_CI = aipw$rct$ate.u,
                          AIPW_obs_ATE_est = aipw$obs$ate,
                          AIPW_obs_ATE_est_low_CI = aipw$obs$ate.l,
                          AIPW_obs_ATE_est_high_CI = aipw$obs$ate.u)
  data <- cbind(data, aipw_data)
  
  ##############################################################################
  # Shrinkage method
  ##############################################################################
  # Create Strata for shrinkage methods
  
  dat.strat <- create_strata(4,dat = rbind(data_exp_b[,-c("U1")], data_obs_b[,-c("U1")]))
  
  shrinkage_est <- shrinkage.est(dat.strat,
                                 Q.SL.library =  c("SL.speedglm","SL.gam"),
                                 g.SL.library =  c("SL.mean","SL.speedglm"),
                                 prop.bound = 0.025,
                                 reweight.obs = reweight.obs)
  shrinkage_est
  
  shrinkage_data <- data.frame(
   green_strawd_shrink_ATE_est = shrinkage_est$gsar2,
   rosenman_shrink_ATE_est = shrinkage_est$kappa2
  )
  data <- cbind(data, shrinkage_data)
  
  
  
  ##############################################################################
  # Oberst method
  ##############################################################################
  
  
  oberst.est <- function(aipw.est, output_lambda = FALSE){
    
    if (any(!(c("rct","obs") %in% names(aipw.est))))
      stop(paste0("Need both rct and obs estimates in aipw.est"))
    
    out <- list()
    
    rct <- aipw.est$rct
    obs <- aipw.est$obs
    
    oberst <- obserst(rct$ate, obs$ate, rct$ate.se^2, obs$ate.se^2)
    
    out$ate <- oberst$theta_oberst
    
    if (output_lambda){
      out$lambda <- oberst$lambda
      
    }
    
    return(out)
  }
  
  oberst_est <- oberst.est(aipw.est = aipw, output_lambda = TRUE)

  oberst_data <- data.frame(
    oberst_ATE_est = oberst_est$ate,
    oberst_lambda = oberst_est$lambda
  )
  data <- cbind(data, oberst_data)
  
  
  ##############################################################################
  # RBesT meta-analysis method
  ##############################################################################
  # Goal: Summarize each study's difference in means (treatment - control).
  

  if (RBesT_inc) {
    # -----------------------------------------------------
    # Complete Script: Separate Priors for Untreated and Treated,
    # and Comparison via pmixdiff and Simulation
    # -----------------------------------------------------
    
    data_obs_list <- lapply(1:n_datasets, function(i) {
      get(paste0("data_obs_b", i))
    })
    
    result <- meta_analysis_summary(data_obs_list, data_exp_b, 
                                    outcome_type = "continuous",
                                    treatment_col = "A",
                                    outcome_col = "Y",
                                    robust_weight = 0.2,
                                    robust_mean = 0)
    
    
    
    # Create a small data frame of results
    RBesT_df <- data.frame(
      posterior_mean    = result[[1]],
      posterior_sd      = result[[2]],
      posterior_2.5pc   = result[[3]],
      posterior_97.5pc  = result[[4]]
    )
    data <- cbind(data, RBesT_df)
    
  }
  
  # ########################################################################################################
  # # ATMLE
  # ########################################################################################################

  if (atmle_inc) {
      # Implement matching then run A-TMLE
      data_exp_b$S <- 1
      data_obs_b$S <- 0
      matched_controls <- matching_controls_fn(data_obs_b, data_exp_b,
                                               treatment = "A", 
                                               match_formula = S ~ X1+X2,
                                               ratio = 1) 
      
      data_combined <- bind_rows(matched_controls, data_exp_b) %>%
        ungroup() %>%
        dplyr::select(-U1) %>%
        as.data.frame()
      
      
      atmle_res <- atmle(
        data          = data_combined,
        S             = "S",
        W             = c("X1", "X2", "X3"),
        A             = "A",
        Y             = "Y",
        controls_only = TRUE,
        family        = "gaussian",
        atmle_pooled  = TRUE,
        verbose       = FALSE
      )
      
      # Create a small data frame of results
      atmle_df <- data.frame(
        atmle_ATE         = atmle_res$est,
        atmle_ATE_low_CI  = atmle_res$lower,
        atmle_ATE_high_CI = atmle_res$upper
      )
      
      # Combine with data
      data <- cbind(data, atmle_df)
      
    }
  
  
  
  
  # ########################################################################################################
  # # PROCOVA
  # ########################################################################################################
  

  if (procova_inc) {
    procova_est <- procova.est(data_exp_b, data_obs_b, 
                               covariates = c("X1", "X2", "X3"),
                               outcome = "Y",
                               treatment = "A",
                               continuous = T,
                               interactions = NULL)
    procova_df <- data.frame(
      procova_ATE                = procova_est$ate,
      procova_ATE_low_CI         = procova_est$ate.l,
      procova_ATE_high_CI        = procova_est$ate.u
    )
    data <- cbind(data, procova_df)
    
  }
  
  ########################################################################################################
  # # Confounding Function
  ########################################################################################################
  
  if (confounding_fn_inc){
    # Just note here which columns are being used, as well as the propensity and outcome models
    data_exp_b <- as.data.table(data_exp_b)
    data_obs_b <- as.data.table(data_obs_b)
    
    rct_data <- dataInput(data_exp_b[, .SD, .SDcols = c("X1", "X2","X3", "A", "Y")], 
                          outcome.model = Y ~ A +X1+ X3 + A:X3,
                          ps.model = A~1)
    
    obs_data <- dataInput(data_obs_b[, .SD, .SDcols = c("X1", "X2","X3", "A", "Y")], 
                          outcome.model = Y ~ A +X1+ X3 + A:X3,
                          ps.model = A ~ X1+X2)
    
    intHTE_obj <- IntegrativeHTEcf::IntHTEcf(data.rct = rct_data,
                                              data.rwe = obs_data,
                                              outcome.method = "gam",
                                              cfName = c("X1", "X2", "X3"),
                                              ps.method = "glm",
                                              n.boot = 100)
    
  cf_df <- data.frame(
    cf_A = unname(intHTE_obj$est.int[1]),
    cf_AX3 = unname(intHTE_obj$est.int[2]),
    cf_A_var = unname(intHTE_obj$cov.est.int[1,1]),
    cf_AX3_var = unname(intHTE_obj$cov.est.int[2,2]),
    cf_A_AX3_cov = unname(intHTE_obj$cov.est.int[1,2]),
    cf_ATE = unname(intHTE_obj$est.int[1]) + mean(data_exp_b$X3) * unname(intHTE_obj$est.int[2])
  )
  data <- cbind(data, cf_df)

  }
  
  ########################################################################################################
  # Add eta estimates
  ########################################################################################################
  
  # Dynamically add columns for res$par
  for (i in seq_along(res$par)) {
    col_name <- paste0("eta_", i)  # Create dynamic column names: eta_1, eta_2, ..., eta_i
    data[[col_name]] <- res$par[i] # Assign res$par[i] to the respective column
  }
  
  for (i in seq_along(res_not_log$par)) {
    col_name <- paste0("eta_not_log", i)  # Create dynamic column names: eta_1, eta_2, ..., eta_i
    data[[col_name]] <- res_not_log$par[i] # Assign res$par[i] to the respective column
  }
  
  data <- cbind(data, as.data.frame(t(egs)))
  
  # View the resulting data frame
  print(data)
  
  return(data)
  
}




