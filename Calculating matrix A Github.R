rm(list = ls())


# List of required packages
required_packages <- c("data.table", "causl", "ManyData", "mvtnorm", "survey", "loo", "MASS", "magrittr", "numDeriv", "dplyr")

# Install any missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
lapply(required_packages, library, character.only = TRUE)



# Causal model parameterization

n_o <- 2500
n_e <- 500


family <- list(1,c(5, 5),1,1)
forms_obs <- forms_exp <- list(Z~1,
                               list(X ~ U, U ~ 1),
                               Y ~ X + U,
                               ~ 1)
forms_exp[[2]][1] <- list(X ~ 1)

pars_exp <- list(Z = list(beta = 0, phi = 1),
                 X = list(beta = c(0, 0.2)),
                 U = list(beta = 0),
                 Y = list(beta = c(0,1, 0.2),phi = 1),
                 cop = list(beta = 0))
pars_obs <- pars_exp
pars_exp$X <- list(beta = 0)

data_obs_b <- as.data.table(causl:::rfrugalParam(n = n_o, formulas = forms_obs, family = family, pars = pars_obs))
data_exp_b <- as.data.table(causl:::rfrugalParam(n = n_e, formulas = forms_exp, family = family, pars = pars_exp))


# Check rough outline of datasets
lm(Y~X*Z, data_obs_b)
lm(Y~X*Z, data_exp_b)


# Update families to remove hidden confounder U for estimation

family <- list(1,c(5),1,1)

forms_obs <- forms_exp <- list(Z ~ 1,
                               list(X ~ Z ),
                               Y ~ X,
                               ~ 1)
forms_exp[[2]][1] <- list(X ~ 1)


forms2 <- list(obs = causl:::tidy_formulas(unlist(forms_obs[-2]), kwd = "cop"),
               exp = causl:::tidy_formulas(unlist(forms_obs[-2]), kwd = "cop"))


full_form <- list(obs =  causl:::merge_formulas(forms2$obs),
                  exp =  causl:::merge_formulas(forms2$exp))

msks <- list(obs = causl:::masks(forms2$obs,family = c(rep(1,6),1),full_form$obs$wh),
             exp = causl:::masks(forms2$exp,family = c(rep(1,6),1),full_form$exp$wh))

vars <- causl:::lhs(unlist(forms2$obs[-length(forms2$obs)]))



# Model with copula fitting

fit_exp <- fit_causl(dat = data_exp_b, formulas = unlist(forms2$exp), family = unlist(family[-2]))

fit_obs <- fit_causl(dat = data_obs_b, formulas = unlist(forms2$obs), family = unlist(family[-2]))

mm <- list(obs = fit_obs$mm, exp = fit_exp$mm)






########################################################################################################
# Define density functions
########################################################################################################


# Define your dataset
data <- data_exp_b  
N <- n_e


# Function to claculate copula density
gaussian_copula_density <- function(u_Z, u_Y, theta_intercept, theta_slope, X) {
  # Compute the linear predictor
  eta <- theta_intercept + theta_slope * X
  
  # Apply the logistic transformation to get rho
  rho <- 2 * (1 / (1 + exp(-eta))) - 1  # Equivalent to 2 * expit(eta) - 1
  
  # Ensure rho is strictly within (-1, 1)
  epsilon <- 1e-10
  rho <- pmin(pmax(rho, -1 + epsilon), 1 - epsilon)
  
  # Adjust u_Z and u_Y to avoid issues with qnorm
  epsilon_cdf <- 1e-10
  u_Z <- pmin(pmax(u_Z, epsilon_cdf), 1 - epsilon_cdf)
  u_Y <- pmin(pmax(u_Y, epsilon_cdf), 1 - epsilon_cdf)
  
  # Compute determinant
  det_Sigma <- 1 - rho^2
  if (any(det_Sigma <= 0)) {
    warning("det_Sigma is non-positive, returning zero density.")
    return(0)
  }
  
  # Convert u_Z and u_Y to quantiles
  q_Z <- qnorm(u_Z)
  q_Y <- qnorm(u_Y)
  
  # Compute the exponent
  exponent <- -(q_Z^2 - 2 * rho * q_Z * q_Y + q_Y^2) / (2 * det_Sigma)
  
  # Compute the density
  density <- (1 / (2 * pi * sqrt(det_Sigma))) * exp(exponent)
  
  # Marginal standard normal densities
  f_Z <- dnorm(q_Z)
  f_Y <- dnorm(q_Y)
  
  # Compute the copula density
  copula_density <- density / (f_Z * f_Y)
  
  return(copula_density)
}




# Define the density function for a single data point
density_function_i <- function(phi, x_i) {
  # Extract parameters
  intercept_Z     <- phi[1]  # (Intercept) for Z ~ 1
  intercept_Y     <- phi[2]  # (Intercept) for Y ~ X
  slope_Y         <- phi[3]  # Coefficient for X in Y ~ X
  theta_intercept <- phi[4]  # (Intercept) for copula ~ X
  theta_slope     <- 0  # Coefficient for X in copula ~ X
  sigma_Z         <- phi[5]  # Residual s.e. for Z
  sigma_Y         <- phi[6]  # Residual s.e. for Y ~ X
  
  # Extract data point values
  z <- x_i$Z
  y <- x_i$Y
  X <- x_i$X
  
  # Marginal density of Z
  f_Z <- dnorm(z, mean = intercept_Z, sd = sigma_Z)
  
  # Conditional density of Y given X
  mu_Y_given_X <- intercept_Y + slope_Y * X
  f_Y_given_X <- dnorm(y, mean = mu_Y_given_X, sd = sigma_Y)
  
  # Compute CDFs for copula
  u_Z <- pnorm(z, mean = intercept_Z, sd = sigma_Z)
  u_Y <- pnorm(y, mean = mu_Y_given_X, sd = sigma_Y)
  
  # Copula density
  copula_dens <- gaussian_copula_density(u_Z, u_Y, theta_intercept, theta_slope, X)
  
  # Joint density at x_i
  density_value <- f_Z * f_Y_given_X * copula_dens
  
  return(density_value)
}




########################################################################################################
# Estimating E(A)
########################################################################################################


# Estimated parameters from fit_causl
phi_hat <- fit_exp$par

# Number of parameters
p <- length(phi_hat)

# Initialize a zero matrix for summing Hessians
A_sum <- matrix(0, nrow = p, ncol = p)



# Loop over all data points
for (i in 1:N) {
  # Extract the ith data point
  x_i <- data_exp_b[i, ]
  
  # Define the density function for the ith data point
  density_function_i_phi <- function(phi) {
    density_function_i(phi, x_i)
  }
  
  # Compute the Hessian matrix at phi_hat for x_i
  A_i <- numDeriv::hessian(func = density_function_i_phi, x = phi_hat, method.args = list(eps = 1e-4))
  
  # Handle potential NaN or Inf values
  if (any(is.nan(A_i)) || any(is.infinite(A_i))) {
    warning(paste("Hessian computation failed at data point", i))
    A_i <- matrix(0, nrow = p, ncol = p)
  }
  
  # Sum the Hessian matrices
  A_sum <- A_sum + A_i
  
  # Optional: Print progress every 100 iterations
  if (i %% 100 == 0) {
    cat("Processed", i, "out of", N, "data points.\n")
  }
}

# Compute the empirical expectation
A_empirical <- matrix(A_sum / N, nrow = nrow(A_sum), ncol = ncol(A_sum))






########################################################################################################
# Using Matrices A, V and W to calculate ELPD 
########################################################################################################
V <-  matrix(fit_exp$sandwich, nrow = nrow(fit_exp$sandwich), ncol = ncol(fit_exp$sandwich))
W <-  matrix(fit_obs$sandwich, nrow = nrow(fit_obs$sandwich), ncol = ncol(fit_obs$sandwich))

phi_o <- matrix(fit_obs$par)
phi_e <- matrix(fit_exp$par)
phi_o_t <- t(phi_o)
phi_e_t <- t(phi_e)


Vinv <- solve.default(V)
Winv <- solve.default(W)





# New method implementation
etas <- seq(0, 1, 0.04)
results <- numeric(length(etas))
term1s <- numeric(length(etas))
term2s <- numeric(length(etas))

for (i in seq_along(etas)) {
  eta <- etas[i]
  
  # If eta_i depends on eta, define it here:
  eta_i <- eta  # Modify this line based on how eta_i is calculated
  
  S <- n_e * Vinv + eta_i * n_o *Winv
  Sinv <- solve.default(S)

  
  # Compute the expression:
  term1 <- (phi_o_t - phi_e_t) %*% t(eta_i * n_o * Winv) %*% Sinv %*% A_empirical %*%
    Sinv %*% (eta_i * n_o * Winv) %*% (phi_o - phi_e)
  # Double check that Sinv is right thing here
  term2 <- sum(diag(A_empirical %*% Sinv))

  results[i] <- term1 + term2
  term1s[i] <- term1
  term2s[i] <- term2
}


# Find the eta with the highest value
max_index <- which.max(results)

eta_max <- etas[max_index]

# Report the eta with the highest value
cat("Eta with the highest value is:", eta_max, "\n")
results

term1s
term2s



# Estimation of ELPD using older method


approx_posterior <- function (fit_exp, fit_obs, eta, n_sample, n_e, n_o) 
{
  if (eta < 0) {
    stop("Must provide a positive eta")
  }
  if (eta > 1) {
    warning("We recommend setting eta between 0 and 1")
  }
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
  theta_0 <- fit_exp$par
  theta_inf <- fit_obs$par
  VW <- n_e * invV + eta * n_o * invW
  if (!any(is.na(VW)) && rcond(VW) > 1e-16) {
    invvW <- solve.default(VW)
  }
  else {
    invvW <- tryCatch(MASS::ginv(VW), error = function(e) {
      cat("ERROR : VW is singular", "\n")
    })
  }
  theta_hat <- invvW %*% (n_e * invV %*% theta_0 + eta * n_o * 
                            invW %*% theta_inf)
  FI <- fit_exp$FI + eta * fit_obs$FI
  if (!any(is.na(FI)) && rcond(FI) > 1e-16) {
    invFI <- solve.default(FI)
  }
  else {
    invFI <- tryCatch(MASS::ginv(FI), error = function(e) {
      cat("ERROR : FI is singular", "\n")
    })
  }
  samples <- mvrnorm(n_sample, theta_hat, invFI)
  return(samples)
}

elpd <- -Inf
opt_samples <- NA
opt_eta <- NA

for (eta in etas) {
  tryCatch(
    {samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = eta, n_sample = 2000,
                                  n_e = n_e, n_o = n_o)
      elpd_sample <- calculate_elpd(samples, data_exp_b[,vars,with = F],mm$exp,
                                    msks = msks$exp,method = "WAIC",
                                    family = unlist(family[-2][-length(family[-2])]),
                                    fam_cop = family[length(family)])
      
      if (elpd < elpd_sample[["estimates"]]["elpd_waic", "Estimate"]) {
        elpd <- elpd_sample[["estimates"]]["elpd_waic", "Estimate"]
        opt_samples <- samples
        opt_eta <- eta
        
      }
    }
    ,error = function(e) {cat("ERROR : RCT_size = ",eta, "eta = ", eta, " skipped", "\n")})
}




opt_eta



