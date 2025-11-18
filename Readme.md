# Data-Fusion Simulation Framework 

This repository contains the code used to **simulate and compare methods for Average Treatment Effect (ATE) estimation when an RCT is supplemented by multiple observational datasets**.  The workflow is designed around the *causl* package and implements both classical and cutting-edge approaches (AIPW, PROCOVA, shrinkage, power-likelihood, Bayesian mixture priors, etc.).

## Repository layout

| File | Role |
|------|------|
| **`Simulation code.R`** | Main driver script. 1) Generates synthetic RCT + observational data and calls a suite of estimators through `sim_func()`, 2) runs large simulation loops over bias scenarios and saves the results to CSV. |
| **`Functions.R`** | Library of helper functions sourced by the main script—e.g. posterior samplers (`approx_posterior*`), optimisation helpers (`new_method_optim`, `old_power_likelihood`), AIPW machinery, PROCOVA, shrinkage utilities, etc. Note: new_opt_method corresponds to MDPL.  |
## Quick start

```r
# 0.  Clone the repo and open an R session in the repo root

# 1.  Source the function library
source("Functions.R")

# 2.  Run a single simulation instance
out <- sim_func(omega = 0.5,      # strength of unmeasured confounding
                seed  = 123,      # RNG seed
                n_o   = 500,      # obs sample size per dataset
                n_e   = 500,      # RCT sample size
                n_datasets = 10)  # number of observational datasets
# 3.  Inspect results (Note: new_opt_method corresponds to MDPL)
head(out)

# 4. Run all code in 'Simulation code.R' to implement method for a variety of different methods and bias values in parallel.




