##########################################################
# Function definitions
########################################################## 

approx_posterior <- function(fit_exp, # model fit objective from fit_causl: MLE and its variance from the experimental data 
                             fit_obs, # model fit objective from fit_causl: MLE and its variance from the observational data 
                             eta, # a given eta 
                             n_sample,  # number of posterior samples to draw
                             use_sandwich = F){ # whether to use the sandwich estimator of the variance of MLEs. Default is FALSE - use the inverse of FI.
  if (eta < 0) {
    stop("Must provide a positive eta")
  }
  if (eta > 1) {
    warning("We recommend setting eta between 0 and 1")
  }
  
  theta_e <- fit_exp$par # MLE of the experimental data parameters
  theta_o <- fit_obs$par # MLE of the observational data parameters
  
  n_e <- nrow(fit_exp$mm) 
  n_o <- nrow(fit_obs$mm)
  
  FI <- fit_exp$FI + eta * fit_obs$FI # fisher information
  
  if (!any(is.na(FI)) && rcond(FI) > 1e-16) {
    invFI <- solve.default(FI)
  }
  else {
    invFI <- tryCatch(MASS::ginv(FI), error = function(e) {
      cat("ERROR : FI is singular", "\n")
    })
  }
  
  if (use_sandwich == TRUE){ # If using the sandwich estimator of variance
    V <- n_e * fit_exp$sandwich
    if (!any(is.na(V)) && rcond(V) > 1e-16) {
      invV <- solve.default(V)
    }
    else {
      invV <- tryCatch(MASS::ginv(V), error = function(e) {
        cat("ERROR : V is singular", "\n")
      })
    }
    W <- n_o * fit_obs$sandwich
    if (!any(is.na(W)) && rcond(W) > 1e-16) {
      invW <- solve.default(W)
    }
    else {
      invW <- tryCatch(MASS::ginv(W), error = function(e) {
        cat("ERROR : W is singular", "\n")
      })
    }
    
    VW <- n_e * invV + eta * n_o * invW
    if (!any(is.na(VW)) && rcond(VW) > 1e-16) {
      invvW <- solve.default(VW)
    }
    else {
      invvW <- tryCatch(MASS::ginv(VW), error = function(e) {
        cat("ERROR : VW is singular", "\n")
      })
    }
    
    theta_hat <- invvW %*% (n_e * invV %*% theta_e + eta * n_o *
                              invW %*% theta_o)
    
  } else{ # If using the inverse of Fisher Information matrix
    theta_hat <- invFI %*% (fit_exp$FI %*% theta_e + eta *fit_obs$FI %*% theta_o)
  }
  
  samples <- MASS::mvrnorm(n_sample, theta_hat, invFI)
  return(samples)
}

approx_posterior_multiple_datasets <- function(fit_exp, fit_obs_list, eta_list, 
                                               n_sample, use_sandwich = F) {
  
  library(MASS) # Ensure MASS is loaded for mvrnorm and ginv
  
  # -----------------------
  # Argument normalization
  # -----------------------
  
  
  # If eta_list is a single numeric, turn it into a vector
  if (!is.vector(eta_list)) {
    eta_list <- c(eta_list)
  }
  
  # If fit_obs_list is not a list, make it a list
  if (!is.list(fit_obs_list)) {
    fit_obs_list <- list(fit_obs_list)
  }
  
  # -----------------------
  # Basic argument checks
  # -----------------------
  
  # Check that eta_list and fit_obs_list have the same length
  if (!(length(eta_list) ==  length(fit_obs_list))) {
    stop("eta_list and fit_obs_list must all have the same length.")
  }
  
  # If any values are below 0, set them to 0:
  if (any(eta_list < 0)) {
    message("Some elements of 'eta_list' were below 0 and have been set to 0.")
    eta_list[eta_list < 0] <- 0
  }
  
  # If any values are above 1, set them to 1:
  if (any(eta_list > 1)) {
    message("Some elements of 'eta_list' were above 1 and have been set to 1.")
    eta_list[eta_list > 1] <- 1
  }
  
  # Check fit_exp structure
  required_fields <- c("par", "FI")
  for (field in required_fields) {
    if (!is.list(fit_exp) || is.null(fit_exp[[field]])) {
      stop(paste("fit_exp must have a '", field, "' component.", sep = ""))
    }
  }
  
  # Check fit_obs_list structure and dimensions
  for (i in seq_along(fit_obs_list)) {
    for (field in required_fields) {
      if (!is.list(fit_obs_list[[i]]) || is.null(fit_obs_list[[i]][[field]])) {
        stop(paste("Each element of fit_obs_list must have a '", field, "' component. ",
                   "Element at index ", i, " is missing it.", sep = ""))
      }
    }
  }
  
  
  # -----------------------
  # Dimension checks
  # -----------------------
  
  
  # Check that par is numeric
  if (!is.numeric(fit_exp$par)) {
    stop("fit_exp$par must be numeric.")
  }
  
  # Check dimensions of FI and sandwich in fit_exp
  if (!is.matrix(fit_exp$FI)) {
    stop("fit_exp$FI must be a matrix")
  }
  
  # Check each fit_obs_i
  for (i in seq_along(fit_obs_list)) {
    fit_obs_i <- fit_obs_list[[i]]
    if (!is.numeric(fit_obs_i$par)) {
      stop(paste("fit_obs_list[[", i, "]]$par must be numeric.", sep = ""))
    }
    if (length(fit_obs_i$par) != length(fit_exp$par)) {
      stop(paste("Parameter dimension mismatch: fit_obs_list[[", i, "]]$par does not match ",
                 "fit_exp$par in length.", sep = ""))
    }
    if (!is.matrix(fit_obs_i$FI) ) {
      stop(paste("fit_obs_list[[", i, "]]$FI  must be a matrix.", sep = ""))
    }
    
    # Check matrix dimensions match fit_exp
    if (!all(dim(fit_obs_i$FI) == dim(fit_exp$FI))) {
      stop(paste("Dimension mismatch: fit_obs_list[[", i, "]]$FI does not match fit_exp$FI.", sep = ""))
    }
  }
  
  if (use_sandwich == T){
    
    # Check dimensions of FI and sandwich in fit_exp
    if (!is.matrix(fit_exp$sandwich)) {
      stop("fit_exp$sandwich must be a matrix.")
    }
    
    # Check each fit_obs_i
    for (i in seq_along(fit_obs_list)) {
      fit_obs_i <- fit_obs_list[[i]]
      
      if (!is.matrix(fit_obs_i$sandwich)) {
        stop(paste("fit_obs_list[[", i, "]]$sandwich must be a matrix.", sep = ""))
      }
      
      # Check matrix dimensions match fit_exp
      if (!all(dim(fit_obs_i$sandwich) == dim(fit_exp$sandwich))) {
        stop(paste("Dimension mismatch: fit_obs_list[[", i, "]]$sandwich does not match fit_exp$sandwich.", sep = ""))
      }
    }
    
  }
  
  
  # Basic definitions
  
  theta_e <- fit_exp$par
  n_e <- nrow(fit_exp$mm) 
  n_o_list <- lapply(1:length(fit_obs_list), function(i) {
    # Retrieve the fit_obs_i object
    fit_obj <- fit_obs_list[[i]]
    
    # Extract the mm component
    nrow(fit_obj$mm)
  })
  
  n_o_list <- unlist(n_o_list)
  
  # If n_o_list is a single numeric, turn it into a vector
  if (!is.vector(n_o_list)) {
    n_o_list <- c(n_o_list)
  }
  
  
  # -----------------------
  # Compute FI
  # -----------------------
  
  FI <- fit_exp$FI
  for (i in seq_along(fit_obs_list)) {
    eta_i <- eta_list[i]
    FI <- FI + eta_i * fit_obs_list[[i]]$FI
  }
  
  if (!is.matrix(FI) || nrow(FI) != ncol(FI)) {
    stop("FI must be a square matrix.")
  }
  
  if (!any(is.na(FI)) && rcond(FI) > 1e-16) {
    invFI <- solve.default(FI)
  } else {
    invFI <- tryCatch({
      ginv(FI)
    }, error = function(e) {
      stop("ERROR: FI is singular and cannot be inverted or pseudo-inverted.")
    })
  }
  
  # Sandwiching method
  if (use_sandwich == T){
    # Calculate matrices
    V <- n_e * fit_exp$sandwich
    
    if (!is.matrix(V) || nrow(V) != ncol(V)) {
      stop("V must be a square matrix.")
    }
    
    if (!any(is.na(V)) && rcond(V) > 1e-16) {
      invV <- solve.default(V)
    } else {
      invV <- tryCatch({
        ginv(V)
      }, error = function(e) {
        stop("ERROR: V is singular and cannot be inverted or pseudo-inverted.")
      })
    }
    
    # For each fit_obs_i in fit_obs_list, compute W_i and invW_i
    W_list <- vector("list", length(fit_obs_list))
    invW_list <- vector("list", length(fit_obs_list))
    theta_inf_list <- vector("list", length(fit_obs_list))
    
    for (i in seq_along(fit_obs_list)) {
      fit_obs_i <- fit_obs_list[[i]]
      eta_i <- eta_list[i]
      n_o_i <- n_o_list[i]
      
      W_i <- n_o_i * fit_obs_i$sandwich
      if (!is.matrix(W_i) || nrow(W_i) != ncol(W_i)) {
        stop(paste("W_i for i=", i, " is not a square matrix.", sep=""))
      }
      
      if (!any(is.na(W_i)) && rcond(W_i) > 1e-16) {
        invW_i <- solve.default(W_i)
      } else {
        invW_i <- tryCatch({
          ginv(W_i)
        }, error = function(e) {
          stop(paste("ERROR: W_i is singular and cannot be inverted or pseudo-inverted for i=", i, ".", sep = ""))
        })
      }
      
      W_list[[i]] <- W_i
      invW_list[[i]] <- invW_i
      theta_inf_list[[i]] <- fit_obs_i$par
    }
    #Compute VW and theta hat
    VW <- n_e * invV
    for (i in seq_along(fit_obs_list)) {
      eta_i <- eta_list[i]
      n_o_i <- n_o_list[i]
      # Add to VW
      VW <- VW + eta_i * n_o_i * invW_list[[i]]
    }
    
    if (!is.matrix(VW) || nrow(VW) != ncol(VW)) {
      stop("VW must be a square matrix.")
    }
    
    if (!any(is.na(VW)) && rcond(VW) > 1e-16) {
      invvW <- solve.default(VW)
    } else {
      invvW <- tryCatch({
        ginv(VW)
      }, error = function(e) {
        stop("ERROR: VW is singular and cannot be inverted or pseudo-inverted.")
      })
    }
    
    
    
    theta_sum <- n_e * (invV %*% theta_e)
    for (i in seq_along(fit_obs_list)) {
      eta_i <- eta_list[i]
      n_o_i <- n_o_list[i]
      theta_inf_i <- theta_inf_list[[i]]
      theta_sum <- theta_sum + eta_i * n_o_i * (invW_list[[i]] %*% theta_inf_i)
    }
    
    if (ncol(invvW) != nrow(theta_sum)) {
      stop("Dimension mismatch in computing theta_hat: invvW and theta_sum are not conformable.")
    }
    
    theta_hat <- invvW %*% theta_sum
  } else {
    # Calculation without sandwiching 
    theta_hat <- invFI %*% (fit_exp$FI %*% theta_e) 
    for (i in seq_along(fit_obs_list)) {
      eta_i <- eta_list[i]
      fit_obs_i <- fit_obs_list[[i]]
      theta_o_i <- fit_obs_i$par
      theta_hat <- theta_hat + invFI %*% (eta_i * fit_obs_i$FI %*% theta_o_i)
    }
  }
  
  
  if (length(theta_hat) != nrow(invFI)) {
    stop("Dimension mismatch: length(theta_hat) and nrow(invFI) must match.")
  }
  
  # -----------------------
  # Generate samples
  # ----------------------
  
  
  samples <- tryCatch({
    MASS::mvrnorm(n_sample, theta_hat, invFI)
  }, error = function(e) {
    stop("Error in mvrnorm: ", e$message)
  })
  
  return(samples)
}




make_pd <- function(M, tol = 1e-12) {
  # Step 1: Apply nearPD to get a positive-semidefinite matrix
  M_pd <- nearPD(M, corr = FALSE)$mat
  
  # Step 2: Check if the eigenvalues are strictly above 'tol'.
  #         If not, add a ridge to push them above 'tol'.
  eigvals <- eigen(M_pd, only.values = TRUE)$values
  min_eig <- min(eigvals)
  
  if (min_eig < tol) {
    # Increase the diagonal by enough to lift the minimum eigenvalue to 'tol'
    bump <- (tol - min_eig)
    diag(M_pd) <- diag(M_pd) + bump
  }
  
  return(as.matrix(M_pd))
}

calculate_elpd <- function(samples, data, mm, msks,
                           method = "WAIC", 
                           inCop, family, fam_cop) {
  if (!(method %in% c("WAIC", "LOO"))) {
    stop("Method has to be either WAIC or LOO.")
  }
  
  # For each row of 'samples', compute the log-likelihood vector:
  lst <- vector("list", length = nrow(samples))
  for (i in seq_len(nrow(samples))) {
    theta2 <- samples[i, ]
    msks2  <- data.table::copy(msks)
    np     <- sum(msks$beta_m > 0)
    
    # Update beta and phi from the parameter draw 'theta2'
    msks2$beta_m[msks$beta_m > 0] <- theta2[seq_len(np)]
    msks2$phi_m[msks$phi_m > 0]   <- theta2[-seq_len(np)]
    
    # Calculate log-likelihood for this draw
    ll_i <- causl:::ll(
      data,
      mm,
      msks2$beta_m,
      phi = msks2$phi_m,
      inCop   = inCop,
      family  = family,
      fam_cop = fam_cop
    )
    
    lst[[i]] <- as.vector(ll_i)
  }
  
  # Combine all draws into a single matrix: (n_draws x n_data_points)
  mtrx <- do.call(rbind, lst)
  
  # --- NEW: remove any rows that contain NA ---
  # i.e. only keep draws that yield valid (non-NA) log-likelihoods
  idx_ok  <- stats::complete.cases(mtrx)
  mtrx_ok <- mtrx[idx_ok, , drop = FALSE]
  
  # If all rows had NA, you'll end up with zero rows
  if (nrow(mtrx_ok) == 0) {
    stop("No valid draws left after removing NA rows.")
  }
  
  # Finally compute WAIC or LOO on the filtered matrix
  if (method == "WAIC") {
    elpd <- loo::waic(mtrx_ok)
  } else {
    elpd <- loo::loo(mtrx_ok)
  }
  
  return(elpd)
}


estimate_A_empirical_new <- function(data, msks, fit_exp, mm, family,
                                     dat_cols = c("C1","Z", "Y"),
                                     progress_every = 100, eps = 1e-8, inCop = 1:2) {
  # Create msks_ll from the 'exp' element of msks and extract initial parameters
  msks_ll <- msks$exp
  theta2 <- fit_exp$par
  N <- nrow(data)
  
  # Create a copy of msks_ll to avoid modifying the original
  # If data.table is available, use its copy function; otherwise, a simple assignment may suffice.
  msks2 <- if (exists("copy", mode = "function")) copy(msks_ll) else msks_ll
  
  # Determine positions where beta and phi are active
  np <- sum(msks_ll$beta_m > 0)
  msks2$beta_m[msks_ll$beta_m > 0] <- theta2[seq_len(np)]
  msks2$phi_m[msks_ll$phi_m > 0] <- theta2[-seq_len(np)]
  
  # Subset the data to only include the relevant columns
  dat_i <- as.data.frame(data)[dat_cols]

  # Define a wrapper for the log-likelihood function that takes parameters and an index i
  loglik_wrapper <- function(params, i) {
    # Use current values from msks2 and update where necessary
    beta_current <- msks2$beta_m
    phi_current <- msks2$phi_m
    beta_current[msks_ll$beta_m > 0] <- params[seq_len(np)]
    phi_current[msks_ll$phi_m > 0] <- params[-seq_len(np)]
    
    # Extract the i-th row from the data
    row_i <- dat_i[i, ]
    
    # Call the causl:::ll function with all required inputs
    ll_value <- causl:::ll(
      dat = row_i,
      mm = mm$exp[i, ],
      beta = beta_current,
      phi = phi_current,
      inCop = inCop,
      fam_cop = 1,
      family = unlist(family[c(1, 3)]),
      use_cpp = TRUE,
      exclude_Z = FALSE
    )
    
    # Return the exponentiated log-likelihood value
    return(ll_value)
  }
  
  # Use the initial parameter estimate as phi_hat and determine its length
  phi_hat <- fit_exp$par
  p <- length(phi_hat)
  
  # Initialize a matrix to accumulate the Hessians
  A_sum <- matrix(0, nrow = p, ncol = p)

  # Loop over each data point to compute and sum the Hessian matrices
  for (i in 1:N) {
    loglik_i_phi <- function(phi) {
      loglik_wrapper(phi, i)
    }
    A_i <- tryCatch({
      numDeriv::hessian(func = loglik_i_phi, x = phi_hat, method.args = list(eps = eps))
    }, error = function(e) {
      warning(sprintf("Hessian computation failed at data point %d: %s", i, e$message))
      return(matrix(0, nrow = p, ncol = p))
    })
    
    # Check for NaN or Inf values and replace with a zero matrix if needed
    if (any(is.nan(A_i)) || any(is.infinite(A_i))) {
      warning(sprintf("Hessian computation returned NaN/Inf at data point %d. Replacing with zero matrix.", i))
      A_i <- matrix(0, nrow = p, ncol = p)
    }
    
    A_sum <- A_sum + A_i
    
    # Optional: print progress every 'progress_every' iterations
    if (i %% progress_every == 0) {
      cat("Processed", i, "out of", N, "data points.\n")
    }
  }
  
  # Compute the empirical A matrix by averaging the summed Hessians
  A_empirical <- A_sum / N
  
  return(A_empirical)
}

estimate_A_empirical_old <- function(data, msks, fit_exp, mm, family,
                                     dat_cols = c("C1","Z", "Y"),inCop = 1:2,
                                     progress_every = 100, eps = 1e-8) {
  # Create msks_ll from the 'exp' element of msks and extract initial parameters
  msks_ll <- msks$exp
  theta2 <- fit_exp$par
  N <- nrow(data)
  
  # Create a copy of msks_ll to avoid modifying the original
  # If data.table is available, use its copy function; otherwise, a simple assignment may suffice.
  msks2 <- if (exists("copy", mode = "function")) copy(msks_ll) else msks_ll
  
  # Determine positions where beta and phi are active
  np <- sum(msks_ll$beta_m > 0)
  msks2$beta_m[msks_ll$beta_m > 0] <- theta2[seq_len(np)]
  msks2$phi_m[msks_ll$phi_m > 0] <- theta2[-seq_len(np)]
  
  # Subset the data to only include the relevant columns
  dat_i <- as.data.frame(data)[dat_cols]
  
  # Define a wrapper for the log-likelihood function that takes parameters and an index i
  loglik_wrapper <- function(params, i) {
    # Use current values from msks2 and update where necessary
    beta_current <- msks2$beta_m
    phi_current <- msks2$phi_m
    beta_current[msks_ll$beta_m > 0] <- params[seq_len(np)]
    phi_current[msks_ll$phi_m > 0] <- params[-seq_len(np)]
    
    # Extract the i-th row from the data
    row_i <- dat_i[i, ]
    
    # Call the causl:::ll function with all required inputs
    ll_value <- causl:::ll(
      dat = row_i,
      mm = mm$exp[i, ],
      beta = beta_current,
      phi = phi_current,
      inCop = inCop,
      fam_cop = 1,
      family = unlist(family[c(1, 3)]),
      use_cpp = TRUE,
      exclude_Z = FALSE
    )
    
    # Return the exponentiated log-likelihood value
    return(exp(ll_value))
  }
  
  # Use the initial parameter estimate as phi_hat and determine its length
  phi_hat <- fit_exp$par
  p <- length(phi_hat)
  
  # Initialize a matrix to accumulate the Hessians
  A_sum <- matrix(0, nrow = p, ncol = p)
  
  # Loop over each data point to compute and sum the Hessian matrices
  for (i in 1:N) {
    loglik_i_phi <- function(phi) {
      loglik_wrapper(phi, i)
    }
    A_i <- tryCatch({
      numDeriv::hessian(func = loglik_i_phi, x = phi_hat, method.args = list(eps = eps))
    }, error = function(e) {
      warning(sprintf("Hessian computation failed at data point %d: %s", i, e$message))
      return(matrix(0, nrow = p, ncol = p))
    })
    
    # Check for NaN or Inf values and replace with a zero matrix if needed
    if (any(is.nan(A_i)) || any(is.infinite(A_i))) {
      warning(sprintf("Hessian computation returned NaN/Inf at data point %d. Replacing with zero matrix.", i))
      A_i <- matrix(0, nrow = p, ncol = p)
    }
    
    A_sum <- A_sum + A_i
    
    # Optional: print progress every 'progress_every' iterations
    if (i %% progress_every == 0) {
      cat("Processed", i, "out of", N, "data points.\n")
    }
  }
  
  # Compute the empirical A matrix by averaging the summed Hessians
  A_empirical <- A_sum / N
  
  return(A_empirical)
}


new_method_optim <- function(fit_exp, fit_obs_list, A_empirical) {
  # Early definitions
  K <- length(fit_obs_list)
  phi_e <- matrix(fit_exp$par)
  n_e <- nrow(fit_exp$mm)
  V <-  matrix(solve.default(fit_exp$FI), nrow = nrow(fit_exp$FI), ncol = ncol(fit_exp$FI)) *n_e
  Vinv <- solve.default(V)
  
  
  
  # Compute the number of observations per dataset from fit_obs_list
  # (each element should have a model matrix 'mm')
  n_o_list <- sapply(fit_obs_list, function(fit_obj) {
    if (!is.null(fit_obj$mm)) nrow(fit_obj$mm)
    else stop("Each fit_obs object must have a 'mm' component.")
  })
  
  # Initialize lists for W, its inverse, and phi_o (the parameter vector for each dataset)
  W_list <- vector("list", K)
  Winv_list <- vector("list", K)
  phi_o_list <- vector("list", K)
  
  # Loop through each dataset and compute W, Winv, and phi_o
  for (i in seq_along(fit_obs_list)) {
    fit_obs_i <- fit_obs_list[[i]]
    n_o_i <- n_o_list[i]
    
    # Check that FI is a square matrix
    if (!is.matrix(fit_obs_i$FI) || nrow(fit_obs_i$FI) != ncol(fit_obs_i$FI)) {
      stop(sprintf("fit_obs_list[[%d]]$FI is not a square matrix.", i))
    }
    
    # Compute W_i as the inverse of FI scaled by n_o_i
    W_i <- solve(fit_obs_i$FI) * n_o_i
    
    # Get phi_o_i as a column vector
    phi_o_i <- matrix(fit_obs_i$par, ncol = 1)
    
    # Compute the inverse of W_i, using a pseudo-inverse if needed
    if (rcond(W_i) > 1e-16) {
      Winv_i <- solve(W_i)
    } else {
      Winv_i <- MASS::ginv(W_i)
    }
    
    # Store the computed matrices and vector
    W_list[[i]] <- W_i
    Winv_list[[i]] <- Winv_i
    phi_o_list[[i]] <- phi_o_i
  }
  
  # Define the objective function f_fun which depends on the vector x (length K)
  f_fun <- function(x) {
    # Construct S = n_e * Vinv + sum_{k=1}^{K} [x[k] * n_o_list[k] * Winv_list[[k]]]
    S <- n_e * Vinv
    for (k in 1:K) {
      S <- S + x[k] * n_o_list[k] * Winv_list[[k]]
    }
    Sinv <- solve(S)
    
    # term1 is a double sum over datasets k and l
    term1 <- 0
    for (k in 1:K) {
      for (l in 1:K) {
        mat_k <- x[k] * n_o_list[k] * Winv_list[[k]]
        mat_l <- x[l] * n_o_list[l] * Winv_list[[l]]
        delta_phi_k <- phi_o_list[[k]] - phi_e
        delta_phi_l <- phi_o_list[[l]] - phi_e
        term <- t(delta_phi_k) %*% t(mat_k) %*% Sinv %*% A_empirical %*% Sinv %*% mat_l %*% delta_phi_l
        term1 <- term1 + term
      }
    }
    
    # term2 is the trace term
    term2 <- sum(diag(A_empirical %*% Sinv))
    
    result <- term1 + term2
    return(as.numeric(result))
  }
  
  # Set up the initial guess and bounds for the optimizer
  x_start <- rep(0, K)
  lower_bounds <- rep(0, K)
  upper_bounds <- rep(1, K)
  
  # Run the optimization using L-BFGS-B
  res <- optim(
    par = x_start,
    fn = f_fun,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(maxit = 1000, trace = 1, factr = 1e-7, pgtol = 1e-8, fnscale = -1)
  )
  
  # Print diagnostic output
  cat("Optimized parameters:\n")
  print(res$par)
  cat("Value of function at the optimum:\n")
  print(res$value)
  cat("Convergence code (0 means successful):\n")
  print(res$convergence)
  
  # Return the optimization result
  return(res)
}



old_power_likelihood <- function(data, fit_exp, fit_obs, mm, msks, family, vars,
                                 eta_seq = c(seq(0,0.05,0.01),seq(0, 1, 0.05)),
                                 n_sample = 2000) {
  # Initialize tracking variables
  best_elpd <- -Inf
  opt_eta_old <- NA
  old_samples <- NA
  

  # Loop over the sequence of eta values
  for (eta in eta_seq) {
    tryCatch({
      # Compute approximate posterior samples for the current eta
      samples <- approx_posterior(fit_exp = fit_exp, fit_obs = fit_obs,
                                  eta = eta, n_sample = n_sample, use_sandwich = FALSE)
      
      # Calculate elpd for these samples. Note: convert data to data.table and subset to the variables in 'vars'.
      elpd_sample <- calculate_elpd(samples,
                                    as.data.table(data)[, vars, with = FALSE],
                                    mm$exp,
                                    msks = msks$exp,
                                    method = "WAIC",
                                    family = unlist(family[-2][-length(family[-2])]),
                                    fam_cop = family[length(family)])
      
      # Update best estimate if the current elpd is higher
      if (best_elpd < elpd_sample[["estimates"]]["elpd_waic", "Estimate"]) {
        best_elpd <- elpd_sample[["estimates"]]["elpd_waic", "Estimate"]
        old_samples <- samples
        opt_eta_old <- eta
      }
    }, error = function(e) {
      cat("ERROR: For eta =", eta, "skipped due to:", e$message, "\n")
    })
  }
  
  cat("Eta under previous estimates is:", opt_eta_old, "\n")
  results <- list(opt_eta_old, old_samples)
  return(results)
}

meta_analysis_summary <- function(data_obs_list, data_new, 
                                  outcome_type = c("binary", "continuous"),
                                  treatment_col,
                                  outcome_col,
                                  robust_weight = 0.2,
                                  robust_mean = 0,
                                  robust_sigma = 3) {
  # Load required packages
  require(dplyr)
  require(data.table)
  require(tidyr)
  
  outcome_type <- match.arg(outcome_type)
  
  if (outcome_type == "binary") {
    ## ===== Binary Outcome Method ===== ##
    # 1. Helper: Summarize subject-level binary data into binomial counts.
    summarize_binomial_data <- function(df, study_label, treatment_col, outcome_col) {
      df %>%
        group_by(!!sym(treatment_col)) %>%
        summarise(
          r = sum(!!sym(outcome_col)),   # number of successes
          n = n(),                       # total subjects in this arm
          .groups = "drop"
        ) %>%
        mutate(study = study_label, .before = 1)
    }
    
    # 2. Summarize each historical study.
    hist_summaries <- lapply(seq_along(data_obs_list), function(i) {
      summarize_binomial_data(data_obs_list[[i]], paste0("obs_b", i),
                              treatment_col, outcome_col)
    })
    hist_data <- bind_rows(hist_summaries)
    
    # 3. Split historical data into placebo (treatment == 0) and treated (treatment == 1)
    hist_placebo <- filter(hist_data, !!sym(treatment_col) == 0)
    hist_treat   <- filter(hist_data, !!sym(treatment_col) == 1)
    
    # 4. Fit separate gMAP models for each group.
    set.seed(123)
    map_hist_placebo <- gMAP(
      formula   = cbind(r, n) ~ 1 | study,
      data      = hist_placebo,
      family    = binomial(link = "logit"),
      tau.dist  = "HalfNormal",
      tau.prior = c(0, 5),
      chains    = 4,
      iter      = 4000,
      warmup    = 1000
    )
    set.seed(123)
    map_hist_treat <- gMAP(
      formula   = cbind(r, n) ~ 1 | study,
      data      = hist_treat,
      family    = binomial(link = "logit"),
      tau.dist  = "HalfNormal",
      tau.prior = c(0, 5),
      chains    = 4,
      iter      = 4000,
      warmup    = 1000
    )
    
    # 5. Convert each historical model to a mixture prior and robustify.
    prior_placebo <- automixfit(map_hist_placebo)
    robust_placebo <- robustify(
      priormix = prior_placebo,
      weight   = robust_weight,
      mean     = robust_mean,
      sigma    = robust_sigma
    )
    prior_treat <- automixfit(map_hist_treat)
    robust_treat <- robustify(
      priormix = prior_treat,
      weight   = robust_weight,
      mean     = robust_mean,
      sigma    = robust_sigma
    )
    
    # 6. Summarize the new trial data for each arm.
    new_placebo <- data_new %>%
      filter(!!sym(treatment_col) == 0) %>%
      summarise(r = sum(!!sym(outcome_col)), n = n())
    new_treat <- data_new %>%
      filter(!!sym(treatment_col) == 1) %>%
      summarise(r = sum(!!sym(outcome_col)), n = n())
    
    cat("New study placebo: r =", new_placebo$r, "n =", new_placebo$n, "\n")
    cat("New study treated: r =", new_treat$r, "n =", new_treat$n, "\n\n")
    
    # 7. Update each prior with the new trial data using postmix().
    post_placebo <- postmix(robust_placebo, r = new_placebo$r, n = new_placebo$n)
    post_treat   <- postmix(robust_treat,   r = new_treat$r,   n = new_treat$n)
    
    # 8. Compute the probability that (treated - placebo) < 0.
    prob_smaller <- pmixdiff(post_treat, post_placebo, 0, lower.tail = FALSE)
    cat("Probability that (treated - placebo) < 0:", prob_smaller, "\n\n")
    
    # 9. Simulate from each updated mixture to obtain the difference distribution.
    n_sim <- 10000
    draws_treat   <- rmix(post_treat, n = n_sim)
    draws_placebo <- rmix(post_placebo, n = n_sim)
    draws_diff    <- draws_treat - draws_placebo
    
    # 10. Compute and return summary statistics.
    posterior_mean  <- unname(mean(draws_diff))
    posterior_sd    <- unname(sd(draws_diff))
    posterior_low   <- unname(quantile(draws_diff, 0.025))
    posterior_high  <- unname(quantile(draws_diff, 0.975))
    
    return(list(posterior_mean = posterior_mean,
                posterior_sd   = posterior_sd,
                Low            = posterior_low,
                High           = posterior_high))
    
  } else if (outcome_type == "continuous") {
    ## ===== Continuous Outcome Method ===== ##
    # 1. Helper: Summarize each two-arm study.
    summarize_two_arm <- function(df, study_label, treatment_col, outcome_col) {
      df_summary <- df %>%
        group_by(!!sym(treatment_col)) %>%
        summarise(
          mean = mean(!!sym(outcome_col)),
          sd   = sd(!!sym(outcome_col)),
          n    = n(),
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = !!sym(treatment_col),
          values_from = c(mean, sd, n),
          names_glue = paste0("{.value}_", treatment_col, "{", treatment_col, "}")
        ) %>%
        mutate(
          diff = .data[[paste0("mean_", treatment_col, "1")]] - .data[[paste0("mean_", treatment_col, "0")]],
          se   = sqrt((.data[[paste0("sd_", treatment_col, "1")]]^2 / .data[[paste0("n_", treatment_col, "1")]]) +
                        (.data[[paste0("sd_", treatment_col, "0")]]^2 / .data[[paste0("n_", treatment_col, "0")]]))
        )
      df_summary$study <- study_label
      df_summary
    }
    
    # 2. Summarize historical studies.
    data_prior_2X <- bind_rows(lapply(seq_along(data_obs_list), function(i) {
      summarize_two_arm(data_obs_list[[i]], paste0("obs_b", i), treatment_col, outcome_col)
    }))
    
    # 3. Fit a MAP prior on differences (with known within‐study variance).
    set.seed(123)
    map_prior_fit_2X <- gMAP(
      formula   = cbind(diff, se) ~ 1 | study,
      data      = data_prior_2X,
      family    = gaussian,
      tau.dist  = "HalfNormal",
      tau.prior = c(0, 5),
      chains    = 4,
      iter      = 10000,
      warmup    = 2000,
      thin      = 2
    )
    
    # 4. Mixture approximation & robustification.
    map_mix_2X <- automixfit(map_prior_fit_2X)
    map_robust_2X <- robustify(
      priormix = map_mix_2X,
      weight   = robust_weight,
      mean     = robust_mean,
      sigma    = robust_sigma
    )
    
    # 5. Summarize the new trial data.
    data_exp_b_summ <- data_new %>%
      group_by(!!sym(treatment_col)) %>%
      summarise(
        mean = mean(!!sym(outcome_col)),
        sd   = sd(!!sym(outcome_col)),
        n    = n(),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = !!sym(treatment_col),
        values_from = c(mean, sd, n),
        names_glue = paste0("{.value}_", treatment_col, "{", treatment_col, "}")
      ) %>%
      mutate(
        new_diff = .data[[paste0("mean_", treatment_col, "1")]] - .data[[paste0("mean_", treatment_col, "0")]],
        new_se   = sqrt((.data[[paste0("sd_", treatment_col, "1")]]^2 / .data[[paste0("n_", treatment_col, "1")]]) +
                          (.data[[paste0("sd_", treatment_col, "0")]]^2 / .data[[paste0("n_", treatment_col, "0")]]))
      )
    
    new_diff <- data_exp_b_summ$new_diff
    new_se   <- data_exp_b_summ$new_se
    
    # 6. Update the prior with the new two-arm data.
    post_exp_2X <- postmix(map_robust_2X, m = new_diff, se = new_se)
    
    # 7. Extract summary statistics from the updated mixture.
    s_post_exp_2X <- summary(post_exp_2X)
    posterior_mean <- unname(s_post_exp_2X[1])
    posterior_sd   <- unname(s_post_exp_2X[2])
    posterior_low  <- unname(s_post_exp_2X[3])
    posterior_high <- unname(s_post_exp_2X[5])
    
    return(list(posterior_mean = posterior_mean,
                posterior_sd   = posterior_sd,
                posterior_low  = posterior_low,
                posterior_high = posterior_high))
  }
}

procova.est <- function(dat.r, dat.o, covariates, 
                        treatment = "ARM", outcome = "Y",
                        prog.SL.library = c("SL.mean", "SL.glmnet", "SL.gam"),
                        Q.SL.library = c("SL.glm", "SL.gam"),
                        g.SL.library = c("SL.glm", "SL.glmnet"),
                        interactions = NULL,    # vector of variable names to interact with treatment
                        continuous = TRUE,
                        confidence.level = 0.05) {
  
  # Convert inputs to data.table if they are not already
  if (!data.table::is.data.table(dat.r)) {
    dat.r <- data.table::as.data.table(dat.r)
  }
  if (!data.table::is.data.table(dat.o)) {
    dat.o <- data.table::as.data.table(dat.o)
  }
  
  # Ensure that covariates exist in dat.r
  if (any(!(covariates %in% names(dat.r))))
    stop("Argument covariates contains values outside the column names of input data.")
  
  # If interactions are provided, ensure they are in dat.r
  if (!is.null(interactions) && length(interactions) > 0) {
    if (any(!(interactions %in% names(dat.r))))
      stop("Some specified interaction variables are not found in the data.")
  }
  
  out <- list()
  
  # Set families based on outcome type
  if (continuous) {
    sl_family <- gaussian()
    glm_family <- gaussian()
  } else {
    sl_family <- binomial(link = "logit")
    glm_family <- binomial(link = "logit")
  }
  
  ## 1. Prognostic SuperLearner Fit
  # Fit the prognostic model using all covariates from controls (treatment == 0)
  sl <- SuperLearner(
    Y = dat.o[dat.o[[treatment]] == 0, get(outcome)],
    X = as.data.frame(scale(dat.o[dat.o[[treatment]] == 0, covariates, with = FALSE],
                            scale = FALSE)),
    family = sl_family,
    SL.library = prog.SL.library
  )
  
  # Get prognostic predictions and center them
  prog_raw <- predict(sl, dat.r[, covariates, with = FALSE], onlySL = TRUE)$pred
  dat.r[, prog := prog_raw - mean(prog_raw)]
  
  ## 2. Prepare Data for the Final Regression Model
  if (!is.null(interactions) && length(interactions) > 0) {
    # Remove interaction variables from scaling; scale the remaining covariates.
    other_covariates <- setdiff(covariates, interactions)
    if (length(other_covariates) > 0) {
      scaled_covs <- scale(dat.r[, other_covariates, with = FALSE], scale = FALSE)
      dat.r.slim <- cbind(
        dat.r[, c(treatment, outcome, "prog"), with = FALSE],
        as.data.frame(scaled_covs)
      )
    } else {
      dat.r.slim <- dat.r[, c(treatment, outcome, "prog"), with = FALSE]
    }
    
    # Add the interaction variables in their original (unscaled) form
    dat.r.slim <- cbind(dat.r.slim, dat.r[, interactions, with = FALSE])
    
    # For each interaction variable, create a new column for the product with treatment.
    # Name these columns "trt_int_<var>".
    for (intvar in interactions) {
      int_colname <- paste("trt_int", intvar, sep = "_")
      dat.r.slim[, (int_colname) := dat.r[[treatment]] * dat.r[[intvar]] ]
    }
  } else {
    # If no interactions are provided, scale all covariates.
    dat.r.slim <- cbind(
      dat.r[, c(treatment, outcome, "prog"), with = FALSE],
      as.data.frame(scale(dat.r[, covariates, with = FALSE], scale = FALSE))
    )
  }
  
  ## 3. Fit the Final Regression Model
  form <- as.formula(paste(outcome, "~ ."))
  procova_m <- glm(form, data = dat.r.slim, family = glm_family)
  
  # Compute sandwich covariance matrix for robust standard errors
  procova_m_sandwich <- sandwich::sandwich(procova_m)
  
  ## 4. Calculate the Average Treatment Effect (ATE)
  if (is.null(interactions) || length(interactions) == 0) {
    # No interactions: ATE is simply the coefficient for treatment.
    procova_se <- as.numeric(sqrt(procova_m_sandwich[treatment, treatment]))
    ate <- as.numeric(procova_m$coefficients[treatment])
  } else {
    # With interactions: ATE = beta_t + sum_j mean(X_j) * beta_{t*X_j}
    mean_interactions <- sapply(interactions, function(var) mean(dat.r[[var]]))
    interaction_coef_names <- paste("trt_int", interactions, sep = "_")
    
    ate <- as.numeric(procova_m$coefficients[treatment] +
                        sum(procova_m$coefficients[interaction_coef_names] * mean_interactions))
    
    # Delta method variance calculation:
    # Let m_vec = [1, mean(X_1), ..., mean(X_k)]
    # and let v be the vector of coefficients for treatment and the interactions.
    m_vec <- c(1, mean_interactions)
    coef_names <- c(treatment, interaction_coef_names)
    cov_sub <- procova_m_sandwich[coef_names, coef_names]
    var_ate <- as.numeric(t(m_vec) %*% cov_sub %*% m_vec)
    procova_se <- sqrt(var_ate)
  }
  
  out$ate   <- ate
  out$ate.u <- ate + qnorm(1 - confidence.level / 2) * procova_se
  out$ate.l <- ate + qnorm(confidence.level / 2) * procova_se
  
  return(out)
}

calculate_summary_flexible <- function(new_opt_samples, indices, data, main_var = "X", name = "") {
  # new_opt_samples: matrix/data frame with samples for the global parameter vector.
  # indices: named list from find_indices_for_var(fit_exp, main_var), e.g. c("Y_X", "Y_X:C1", "cop_X", …)
  # data: data frame containing the covariates used in interactions (e.g., "C1")
  # main_var: the variable of interest (e.g., "X")
  # name: a character string to prepend to every column name (if empty, names remain unchanged)
  
  # 1. Compute estimates and sds for every coefficient that involves the main variable.
  coef_names <- names(indices)
  coef_est <- numeric(length(coef_names))
  coef_sd  <- numeric(length(coef_names))
  
  for (i in seq_along(indices)) {
    col_idx <- indices[[i]]
    samples_col <- new_opt_samples[, col_idx]
    coef_est[i] <- mean(samples_col)
    coef_sd[i]  <- sd(samples_col)
  }
  
  # Build a data frame with one column per coefficient (estimate and sd)
  coef_df <- data.frame(matrix(ncol = 0, nrow = 1))
  for (i in seq_along(indices)) {
    key <- coef_names[i]
    # Clean up the key for naming (replace ":" with "_" for valid column names)
    clean_key <- gsub(":", "_", key)
    coef_df[[paste0(clean_key, "_est")]] <- coef_est[i]
    coef_df[[paste0(clean_key, "_sd")]]  <- coef_sd[i]
  }
  
  # 2. Compute the ATE using the Y-model coefficients.
  # Look at all indices for keys that begin with "Y_"
  y_keys <- coef_names[grep("^Y_", coef_names)]
  
  # Look for the main effect for X in the Y-model.
  main_effect_key <- paste0("Y_", main_var)
  if (! main_effect_key %in% y_keys) {
    warning("Main effect ", main_effect_key, " not found among Y coefficients. ATE will be NA.")
    ate_est <- NA
    ate_sd  <- NA
    ate_low_CI  <- NA
    ate_high_CI <- NA
  } else {
    # Get the simulation draws for the main effect.
    main_idx <- indices[[main_effect_key]]
    main_samples <- new_opt_samples[, main_idx]
    
    # Initialize ATE samples as the main effect.
    ate_samples <- main_samples
    
    # Identify interaction terms among the Y-model coefficients (those with a colon in the name).
    interaction_keys <- y_keys[grep(":", y_keys)]
    if (length(interaction_keys) > 0) {
      for (ikey in interaction_keys) {
        # Remove the "Y_" prefix. We assume the remaining string is of the form "X:Other" or "Other:X".
        remainder <- sub("^Y_", "", ikey)
        # Split by colon.
        parts <- unlist(strsplit(remainder, ":"))
        # Determine the other variable (the one that is not main_var).
        other_var <- setdiff(parts, main_var)
        if (length(other_var) == 0) {
          warning("Could not determine other variable for interaction ", ikey)
          next
        }
        other_var <- other_var[1]
        # Check that the other variable exists in the data.
        if (!(other_var %in% names(data))) {
          warning("Covariate ", other_var, " not found in data. Skipping interaction ", ikey)
          next
        }
        # Compute the mean of the other covariate.
        other_mean <- mean(data[[other_var]])
        # Get the simulation draws for the interaction effect.
        inter_idx <- indices[[ikey]]
        inter_samples <- new_opt_samples[, inter_idx]
        # Add to the ATE samples: interaction effect * mean(other variable)
        ate_samples <- ate_samples + inter_samples * other_mean
      }
    }
    
    # Summarize the ATE simulation draws.
    ate_est <- mean(ate_samples)
    ate_sd  <- sd(ate_samples)
    ate_low_CI  <- unname(quantile(ate_samples, 0.025))
    ate_high_CI <- unname(quantile(ate_samples, 0.975))
  }
  
  # 3. Combine the coefficient summaries and the ATE summary into one data frame row.
  out_df <- coef_df
  out_df$ATE_est       <- ate_est
  out_df$ATE_sd        <- ate_sd
  out_df$ATE_low_CI    <- ate_low_CI
  out_df$ATE_high_CI   <- ate_high_CI
  
  # If a prefix is provided, prepend it to every column name.
  if(nzchar(name)) {
    colnames(out_df) <- paste(name, colnames(out_df), sep = "_")
  }
  
  return(out_df)
}

find_indices_for_var <- function(fit_obj, var, tol = 1e-8) {
  indices <- list()
  # Loop over each component in fit_obj$pars (e.g., "Y", "cop", "Z")
  for(comp in names(fit_obj$pars)) {
    comp_obj <- fit_obj$pars[[comp]]
    if (!is.null(comp_obj$beta)) {
      local_beta <- comp_obj$beta
      # If beta is a matrix, use its rownames; if it's a vector, use names()
      if (is.matrix(local_beta)) {
        beta_names <- rownames(local_beta)
        n_local <- nrow(local_beta)
      } else if (!is.null(names(local_beta))) {
        beta_names <- names(local_beta)
        n_local <- length(local_beta)
      } else {
        next  # skip if no names available
      }
      
      # Loop over the local coefficients
      for(j in seq_len(n_local)) {
        coeff_name <- beta_names[j]
        # If this coefficient's name mentions the variable of interest:
        if (grepl(var, coeff_name)) {
          # Get the coefficient value (from a matrix or a vector)
          local_val <- if (is.matrix(local_beta)) local_beta[j, 1] else local_beta[j]
          # Find the matching global index in fit_obj$par (using tolerance)
          match_idx <- which(abs(fit_obj$par - local_val) < tol)
          if (length(match_idx) >= 1) {
            global_index <- match_idx[1]
            indices[[paste(comp, coeff_name, sep = "_")]] <- global_index
          } else {
            warning("No global match found for ", comp, " coefficient ", coeff_name)
          }
        }
      }
    }
  }
  return(indices)
}

matching_controls_fn <- function(data_obs_b, data_exp_b, 
                                 treatment = "X", match_formula = S ~ C1 + Z,
                                 ratio = 1){
  data_obs_b$S <- 0
  data_exp_b$S <- 1
  
  matching_externals <- data_obs_b[data_obs_b[[treatment]] == 0,]
  matching_internals <- data_exp_b[data_exp_b[[treatment]] == 0,]
  
  
  # Combine the datasets (assuming they have similar covariates)
  combined_data <- rbind(matching_internals, matching_externals)
  
  # Perform nearest neighbor matching
  match_out <- matchit(match_formula, data = combined_data, method = "nearest",
                       ratio = ratio)
  
  # Extract the matched dataset
  matched_data <- match.data(match_out)
  
  # Identify the original columns
  original_cols <- names(combined_data)
  
  # Subset matched_data to only include the original columns
  matched_data <- matched_data[, ..original_cols]
  
  # Select the matched observational controls
  matched_controls <- subset(matched_data, S == 0)
  
  return(matched_controls)
}


generate.data.list.HTE <- function(n.r, n.o, seed,
                                   coeff.U.list,
                                   n_datasets = length(coeff.U.list),  # <- NEW
                                   HTE = FALSE,
                                   rct.treat.prob = 0.5,
                                   X.diff = FALSE,
                                   U.diff = FALSE,
                                   U1.mean.shift = 0) {
  
  ## ----- 1.  House-keeping ------------------------------------------------ ##
  # sensible default
  if (length(coeff.U.list) != n_datasets)
    coeff.U.list <- rep(coeff.U.list, length.out = n_datasets)    # recycle
  
  data.list <- vector("list", n_datasets + 1)   # +1 for the single dat.r
  
  family <- list(c(1), c(5, 5,5), 1, 1)        # used for *all* arms
  
  ## ----- 2.  Build the single randomised arm ----------------------------- ##
  forms_exp_d <- list(
    c(X1 ~ 1),
    list(A ~ 1, X3 ~ 1, X2 ~ 1),                # treatment randomised
    Y ~ A + X3 + A:X3 +X2+ U1,
    ~ 1
  )
  
  pars_exp_d <- list(
    X1 = list(beta = 0,          phi = 1),
    X2 = list(beta = 0,          phi = 1),
    X3 = list(beta = 0),
    A  = list(beta = logit(rct.treat.prob)),
    Y  = list(beta = c(0.5, 0.2, 1, 1, 1, 0), phi = 1),
    cop = list(beta = matrix(c(0.5), nrow = 1))
  )
  
  ## ... covariate‐shift options for X3
  if (X.diff) {
    pars_exp_d$X1 <- list(beta = 0.5, phi = 1)
    pars_exp_d$X3 <- list(beta = logit(0.6))
  }
  if (HTE) pars_exp_d$Y$beta[length(pars_exp_d$Y$beta)] <- 0.5
  
  set.seed(seed)
  dat.r <- data.frame(U1 = rnorm(n.r, 0, 1))
  dat.r <- as.data.table(
    causl:::rfrugalParam(n        = n.r,
                         formulas = forms_exp_d,
                         family   = family,
                         pars     = pars_exp_d,
                         dat      = dat.r))
  dat.r[, `:=`(
    NCO.strong = 1 + X2 + X3 + U1 + rnorm(n.r, 0, 2),
    S          = 1,                     # <- flag for the randomised arm
    coeff.U    = NA_real_,
    seed       = seed
  )]
  
  data.list[[1]] <- dat.r               # first list element is always dat.r
  
  
  ## ----- 3.  Loop over observational arms ------------------------------- ##
  for (k in seq_len(n_datasets)) {
    
    coeff.U <- coeff.U.list[k]
    
    forms_obs_d <- list(
      c(X1 ~ 1),
      list(A ~ X1 + X2 + U1, X3 ~ 1, X2 ~ 1),
      Y ~ A + X3 + A:X3 + X2 + U1,
      ~ 1
    )
    
    pars_obs_d <- list(
      X1 = list(beta = 0,          phi = 1),
      X2 = list(beta = 0,          phi = 1),
      X3 = list(beta = 0),
      A  = list(beta = c(0, 1, 1, coeff.U)),
      Y  = list(beta = c(0.5, 0.2, 1, 1,1, 0), phi = 1),
      cop = list(beta = matrix(c(0.5), nrow = 1))
    )
    
    if (X.diff) {
      pars_obs_d$X1 <- list(beta = 0.5, phi = 1)
      pars_obs_d$X3 <- list(beta = logit(0.6))
    }
    if (HTE) pars_obs_d$Y$beta[length(pars_obs_d$Y$beta)] <- 0.5
    
    ## --- Generate U1 with / without mean-shift
    if (U.diff) {
      dat.o <- data.frame(U1 = rnorm(n.o, U1.mean.shift, 1))
    } else {
      dat.o <- data.frame(U1 = rnorm(n.o, 0, 1))
    }
    
    dat.o <- as.data.table(
      causl:::rfrugalParam(n        = n.o,
                           formulas = forms_obs_d,
                           family   = family,
                           pars     = pars_obs_d,
                           dat      = dat.o))
    dat.o[, `:=`(
      NCO.strong = 1 + X2 + X3 + U1 + rnorm(n.o, 0, 2),
      S          = k + 1,               # <- 2, 3, 4, … for the obs arms
      coeff.U    = coeff.U,
      seed       = seed
    )]
    
    data.list[[k + 1]] <- dat.o
  }
  
  return(data.list)                     # list[[1]]=dat.r, list[[2]]…=dat.o
}


aipw_calc <- function(dat,covariates, g.SL.library, Q.SL.library, g.bound = 0.025,
                      X3.reweight = FALSE){
  
  ps.sl <-  SuperLearner(Y = dat$A, X = as.data.frame(dat[,covariates,with = F]),
                         family = binomial(),
                         SL.library = g.SL.library)
  ps <- predict(ps.sl, dat[,covariates, with = F], onlySL = TRUE)$pred
  ps.trunc <- pmin(pmax(ps,g.bound), 1 - g.bound)
  
  mu.1.sl <- SuperLearner(Y = dat[A == 1,]$Y, X = as.data.frame(dat[A == 1,covariates,with = F]),
                          family = gaussian(),
                          SL.library = Q.SL.library)
  
  mu.0.sl <- SuperLearner(Y = dat[A == 0,]$Y, X = as.data.frame(dat[A == 0,covariates,with = F]),
                          family = gaussian(),
                          SL.library = Q.SL.library)
  
  dat[, ':='(ps = ps.trunc,
             mu.1 = predict(mu.1.sl, dat[,covariates, with = F], onlySL = TRUE)$pred,
             mu.0  = predict(mu.0.sl, dat[,covariates, with = F], onlySL = TRUE)$pred)]
  
  dat[,aipw := A*Y/ps - (1 - A)*Y/(1 - ps) - ((A - ps)/(ps*(1 - ps))) * ((1-ps)*mu.1 + ps * mu.0) ]
  dat[,X3.wt := X3 * 0.5/mean(dat$X3) + (1-X3) * 0.5/(1-mean(dat$X3)) ]
  
  if (X3.reweight){
    dat[,aipw := aipw *  X3.wt]
  }
  
  ate <- mean(dat$aipw)
  ate.se <- sqrt(var(dat$aipw)/nrow(dat))
  return(list(ate = ate,
              se = ate.se,
              dat.ps = dat$ps,
              mu.1.sl = mu.1.sl,
              mu.0.sl = mu.0.sl))
  
}



aipw.est <- function(dat,covariates,
                     Q.SL.library =  c("SL.mean","SL.speedglm"),
                     g.SL.library =  c("SL.mean","SL.speedglm"),
                     prop.bound = 0.025,
                     estimates = c("RCT","obs","pooled"),
                     reweight.obs = FALSE,
                     confidence.level = 0.05){
  
  if (any(!(estimates %in% c("RCT","obs","pooled"))))
    stop(paste0("Argument estimates should take values within RCT, obs or pooled"))
  if (any(!(covariates %in% names(dat))))
    stop(paste0("Argument covariates contains values outside the column names of input data."))
  
  dat.r <- dat[S == 1]
  dat.o <- dat[S == 0]
  
  out <- list()
  
  if ("RCT" %in% estimates) { 
    rct <- list()
    
    # RCT
    aipw.r <- aipw_calc(dat = dat.r,
                        covariates = covariates,
                        Q.SL.library = Q.SL.library,
                        g.SL.library = g.SL.library,
                        g.bound = prop.bound)
    
    rct$ate <- aipw.r$ate
    rct$ate.se <- aipw.r$se
    rct$ate.l <- rct$ate - qnorm(1 - confidence.level/2) * aipw.r$se
    rct$ate.u <- rct$ate + qnorm(1 - confidence.level/2) * aipw.r$se
    
    out$rct <- rct
    
  }
  
  # Observational data 
  
  if ("obs" %in% estimates) { 
    
    obs <- list()
    
    aipw.o <- aipw_calc(dat = dat.o,
                        covariates = covariates,
                        Q.SL.library = Q.SL.library,
                        g.SL.library = g.SL.library,
                        g.bound = prop.bound,
                        X3.reweight = reweight.obs)
    
    obs$ate <- aipw.o$ate
    obs$ate.se <- aipw.o$se
    obs$ate.l <- obs$ate - qnorm(1 - confidence.level/2) * aipw.o$se
    obs$ate.u <- obs$ate + qnorm(1 - confidence.level/2) * aipw.o$se
    
    
    out$obs <- obs
  }
  
  if ("pooled" %in% estimates) {
    
    pooled <- list()
    
    aipw.pooled <- aipw_calc(dat = dat,
                             covariates = covariates,
                             Q.SL.library = Q.SL.library,
                             g.SL.library = g.SL.library,
                             g.bound = prop.bound)
    
    pooled$ate <- aipw.pooled$ate
    pooled$ate.se <- aipw.pooled$se
    pooled$ate.l <- pooled$ate - qnorm(1 - confidence.level/2) * aipw.pooled$se
    pooled$ate.u <- pooled$ate + qnorm(1 - confidence.level/2) * aipw.pooled$se
    
    out$pooled <- pooled
  }
  
  return(out)
}

create_strata <- function(K,dat){
  
  dat.r <- dat[S == 1]
  dat.o <- dat[S == 0]
  
  q <- seq(0,1,2/K)[-1]
  table <- as.data.table(rbind(as.data.table(list(X3 = 0,
                                                  bin = q,
                                                  X1 = quantile(dat.o[X3 == 0]$X1,q))),
                               as.data.table(list(X3 = 1,
                                                  bin = q,
                                                  X1 = quantile(dat.o[X3 == 1]$X1,q)))))
  table[bin == 1, X1 := Inf][,K := 1:.N]
  
  # stratify data
  setkey(table,X3,X1)
  setkey(dat.r,X3,X1)
  setkey(dat.o,X3,X1)
  
  dat.r <- table[dat.r, roll = -Inf]
  dat.o <- table[dat.o, roll = -Inf]
  
  out <- as.data.table(rbind(dat.r,dat.o))
}

shrinkage.est <- function(stratified.dat,
                          Q.SL.library =  c("SL.mean","SL.speedglm"),
                          g.SL.library =  c("SL.mean","SL.speedglm"),
                          # splits = 3,
                          prop.bound = 0.025,
                          reweight.obs = FALSE){
  
  if (!("K" %in% names(stratified.dat)))
    stop(paste0("Input data needs to be stratified with strata indicator K"))
  
  K.list <- unique(stratified.dat$K)
  
  K.count <- length(K.list)
  
  if (K.count < 4)
    stop(paste0("Need more than three strata."))
  
  dat.r <- stratified.dat[S == 1]
  
  
  out <- list()
  
  ate.r.vec <- vector()
  ate.se.r.vec <- vector()
  ate.o.vec <- vector()
  ate.se.o.vec <- vector()
  
  for (i in 1:K.count) {
    k <- K.list[i]
    aipw.k <- aipw.est(stratified.dat[K==k,],
                       covariates = c("X1","X2"),
                       Q.SL.library =  Q.SL.library,
                       g.SL.library =  g.SL.library,
                       prop.bound = prop.bound,
                       estimates = c("RCT","obs"),
                       reweight.obs = reweight.obs)
    ate.r.vec[i] <- aipw.k$rct$ate
    ate.se.r.vec[i] <- aipw.k$rct$ate.se
    ate.o.vec[i] <- aipw.k$obs$ate
    ate.se.o.vec[i] <- aipw.k$obs$ate.se
    
    
  }
  
  # variance matrix
  shringkage.var <- as.vector(ate.se.r.vec^2) * diag(K.count)
  # shrinkage
  gsar2.vec <- gsar2(ate.r.vec, ate.o.vec, Sigma_e = shringkage.var)
  kappa2.vec <- rbob_kappa2(ate.r.vec, ate.o.vec, Sigma_e = shringkage.var)
  # reweight 
  wt <- dat.r[,.N,K]$N
  out$gsar2 <- as.vector(gsar2.vec) %*% wt/nrow(dat.r)
  out$kappa2 <- as.vector(kappa2.vec) %*% wt/nrow(dat.r)
  
  return(out)
}




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
  # New method (MDPL) using numerical optimisation
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
  # Old single dataset power likelihood method 
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
  
  # View the resulting data frame
  print(data)
  
  return(data)
  
}

add_all_mse <- function(result_dfs, true_ATE = 0.5) {
  
  # helper to compute squared error
  se <- function(x) (x - true_ATE)^2
  
  # map from existing estimate columns → new MSE column names
  col_map <- c(
    new_opt_ATE_est            = "MSE_new_opt",
    new_ATE_est                = "MSE_new",
    old_ATE_est                = "MSE_old",
    simple_meta_ATE_est        = "MSE_simple_meta",
    rct_ATE_est                = "MSE_rct",
    cf_ATE                     = "MSE_cf",
    posterior_mean             = "MSE_meta_prior",
    atmle_ATE                  = "MSE_atmle",
    procova_ATE                = "MSE_procova",
    oberst_ATE_est             = "MSE_oberst",
    green_strawd_shrink_ATE_est= "MSE_shrink_gsar",
    rosenman_shrink_ATE_est    = "MSE_shrink_rosenman",
    AIPW_ATE_est               = "MSE_aipw"
  )
  
  lapply(result_dfs, function(df) {
    
    ## 1. drop rows where *every* entry is zero
    df <- df[rowSums(df != 0) > 0, , drop = FALSE]
    
    ## 2. add MSE columns wherever the corresponding estimate exists
    for (orig in names(col_map)) {
      if (orig %in% names(df)) {
        df[[col_map[[orig]]]] <- se(df[[orig]])
      }
    }
    
    df
  })
}
