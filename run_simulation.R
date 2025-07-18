####################
# SIMULATIONS
####################
# In the following, the performance of the estimators are investigated based
# on simulations. Samples are generated from the Fréchet copula with parameter
# alpha.

source("simulations.R")
source("main_estimators.R")

# Define parameters
N <- 250 # Sample sizes
p1 <- 0.4 # P(X = 0)
p2 <- 0.9 # P(Y = 0)
alpha <- 0 # Copula parameter
sim <- 1000 # Number of simulations

# Perform simulations
results <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)

# Compute statistics
stat <- corr_bzid_simulation_stats(results, alpha, p1, p2)
