########################################
# PERFORM SIMULATIONS
########################################
# The performance of the estimators are investigated based on simulations.
# Samples are generated from the Fréchet copula with parameter alpha.

# Load dependencies
source("simulation_functions.R")
source("main_estimators.R")

library(doParallel)
registerDoParallel(19)

# Reproducability
set.seed(1)

#########################
# Global constants
#########################

sim <- 1000 # Number of runs

########################################
# MAIN ESTIMATORS
########################################

N <- 150 # Sample size
bsim <- 1000 # Number of resamples in bootstrap

# The following simulatiosn are repeated for
# N = 250, 500, 1000 with bsim = 0 without
# resetting the seed

#########################
# INFLATION COMBINATION 1
#########################

p1 <- 0
p2 <- 0

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)


#########################
# INFLATION COMBINATION 2
#########################

p1 <- .2
p2 <- .2

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

#########################
# INFLATION COMBINATION 3
#########################

p1 <- .2
p2 <- .8

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

#########################
# INFLATION COMBINATION 4
#########################

p1 <- .8
p2 <- .8

# Low dependence
alpha <- .2
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# Mixed dependence
alpha <- .5
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)

# High dependence
alpha <- .8
df.sim <- corr_bzid_simulation_main(N, sim, bsim, alpha, p1, p2)
df.tot <- corr_bzid_simulation_stats(df.sim, alpha, p1, p2)



########################################
# COMPARISON OF SPEARMAN'S RHO
########################################

set.seed(1)

N <- 150 # Sample size
sim <- 1000 # Number of simulations

#########################
# INFLATION COMBINATION 1
#########################

p1 <- .2
p2 <- .2

# Low dependence
alpha <- .2
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# Mixed dependence
alpha <- .5
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# High dependence
alpha <- .8
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

#########################
# INFLATION COMBINATION 2
#########################

p1 <- .2
p2 <- .8

# Low dependence
alpha <- .2
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# Mixed dependence
alpha <- .5
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# High dependence
alpha <- .8
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

#########################
# INFLATION COMBINATION 3
#########################

p1 <- .5
p2 <- .2

# Low dependence
alpha <- .2
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# Mixed dependence
alpha <- .5
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# High dependence
alpha <- .8
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)


#########################
# INFLATION COMBINATION 4
#########################

p1 <- .5
p2 <- .8

# Low dependence
alpha <- .2
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# Mixed dependence
alpha <- .5
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# High dependence
alpha <- .8
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

#########################
# INFLATION COMBINATION 5
#########################

p1 <- .8
p2 <- .8

# Low dependence
alpha <- .2
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# Mixed dependence
alpha <- .5
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)

# High dependence
alpha <- .8
dt.rho <- corr_bzid_spm_rho_main(N, sim, alpha, p1, p2)
corr_bzid_spm_rho_stats(dt.rho, alpha, p1, p2)
