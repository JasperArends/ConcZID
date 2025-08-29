########################################
# PERFORM SIMULATIONS
########################################
# The performance of the estimators are investigated based on simulations.
# Samples are generated from the Fréchet copula with parameter alpha.

# Load dependencies
source("simulation_functions.R")
source("main_estimators.R")

# Reproducability
set.seed(1)

#########################
# Global constants
#########################
N <- 150 # Number of observations
sim <- 1000 # Number of runs

########################################
# MAIN ESTIMATORS
########################################

#########################
# INFLATION COMBINATION 1
#########################

p1 <- .2
p2 <- .2

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

#########################
# INFLATION COMBINATION 2
#########################

p1 <- .2
p2 <- .8

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

#########################
# INFLATION COMBINATION 3
#########################

p1 <- .8
p2 <- .8

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)
corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

########################################
# BOUNDS ESTIMATION
########################################

#########################
# INFLATION COMBINATION 1
#########################

# Inflation probabilities
p1 <- .2
p2 <- .2

# Simulations
bnd_res <- bnd_bzid_simulation_main(N, sim, p1, p2)

# Statistics
bnd_bzid_simulation_stats(bnd_res, p1, p2)        

#########################
# INFLATION COMBINATION 2
#########################

# Inflation probabilities
p1 <- .2
p2 <- .8

# Simulations
bnd_res <- bnd_bzid_simulation_main(N, sim, p1, p2)

# Statistics
bnd_bzid_simulation_stats(bnd_res, p1, p2)  

#########################
# INFLATION COMBINATION 3
#########################

# Inflation probabilities
p1 <- .8
p2 <- .8

# Simulations
bnd_res <- bnd_bzid_simulation_main(N, sim, p1, p2)

# Statistics
bnd_bzid_simulation_stats(bnd_res, p1, p2)  
