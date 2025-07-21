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
p1 <- 0.2 # P(X = 0)
p2 <- 0.6 # P(Y = 0)
alpha <- 0 # Copula parameter
sim <- 1000 # Number of simulations

# Perform simulations
results <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)

# Compute statistics
stat <- corr_bzid_simulation_stats(results, alpha, p1, p2)

# Plot results
library(ggplot2)
library(gridExtra)

grid.arrange(
  # Gini's gamma
  ggplot(results) +
    geom_boxplot(aes(y=Gini_gamma), outlier.size=.2, linewidth=.1) +
    geom_hline(yintercept=as.double(stat$True[1]),
               col='red', linetype='dashed') +
    theme_classic(),
  
  # Spearman's footrule
  ggplot(results) +
    geom_boxplot(aes(y=Spm_foot), outlier.size=.2, linewidth=.1) +
    geom_hline(yintercept=as.double(stat$True[2]),
               col='red', linetype='dashed') +
    theme_classic(),
  
  # Spearman's rho
  ggplot(results) +
    geom_boxplot(aes(y=Spm_rho), outlier.size=.2, linewidth=.1) +
    geom_hline(yintercept=as.double(stat$True[3]),
               col='red', linetype='dashed') +
    theme_classic()
)
