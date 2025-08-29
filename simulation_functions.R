########################################
# SIMULATION FUNCTIONS
########################################

####################
# IMPORT LIBRARIES
####################
library(Metrics)

####################
# GENERATE SAMPLE
####################
gen_sample <- function(N, alpha, p1, p2) {
  # Generates a sample from the Fréchet-copula
  U <- matrix(nrow=N, ncol=2) # Sample
  U[,1] <- runif(N) # X
  
  Z <- runif(N, min=0, max=1) # Decision variable
  V <- runif(N, min=0, max=1) # Independent sample
  
  U[,2] <- (Z <= alpha & alpha > 0) * U[,1] + # Upper Fréchet-Hoeffding copula
    (Z <= -alpha & alpha < 0) * (1 - U[,1]) + # Lower Fréchet-Hoeffding copula
    (Z > abs(alpha)) * V # Independence copula
  
  # Induce inflation
  U[U[,1] <= p1, 1] <- 0
  U[U[,2] <= p2, 2] <- 0
  
  return ( U )
}


####################
# MAIN ESTIMATORS
####################

corr_bzid_simulation_main <- function(N, sim, alpha, p1, p2) {
  # Ensure N and sim and alpha are numeric
  if (!is.numeric(N) || !is.numeric(sim) || !is.numeric(alpha)) {
    stop("N, sim and alpha must be numeric.")
  }
  
  # Ensure p1 and p2 are numeric
  if (!is.numeric(p1) || !is.numeric(p2)) {
    stop("p1 and p2 must be numeric.")
  }
  
  # Create data frames
  results <- data.frame(matrix(nrow=0, ncol=4))
  colnames(results) <- c("Simulation", "Gini_gamma", "Spm_foot", "Spm_rho")
  
  # Perform simulation
  for (ix in 1:sim) {
    # Generate sample
    X <- gen_sample(N, alpha, p1, p2)
    
    # Assume p1 <= p2
    if (sum(X[,2] == 0) < sum(X[,1] == 0))
      X <- X[,c(2, 1)]
    
    # Compute estimates
    results[ix,] <- c(ix,
                      gini_gamma.est(X[,1], X[,2]),
                      spm_foot.est(X[,1], X[,2]),
                      spm_rho.est(X[,1], X[,2]))
  }
  
  return ( results )
}

corr_bzid_simulation_stats <- function(results, alpha, p1, p2) {
  # Compute true values
  conc_true <- list(gini_gamma = gini_gamma.bound(p1,p2)$up * alpha * (alpha > 0) -
                      gini_gamma.bound(p1,p2)$lw * alpha * (alpha < 0),
                    spm_foot = spm_foot.bound(p1,p2)$up * alpha * (alpha > 0) -
                      spm_foot.bound(p1,p2)$lw * alpha * (alpha < 0),
                    spm_rho = spm_rho.bound(p1,p2)$up * alpha * (alpha > 0) -
                      spm_rho.bound(p1,p2)$lw * alpha * (alpha < 0))
  
  # Statistics
  corr_stat <- data.frame(matrix(nrow=3, ncol=7))
  colnames(corr_stat) <- c("Variable", "True", "Mean", "Variance", "MSE", "Left_95CI", "Right_95CI")
  
  # Gini's gamma
  corr_stat[1,] <- c("Gini's gamma",
                     conc_true$gini_gamma,
                     mean(results[,2]),
                     var(results[,2]),
                     mse(conc_true$gini_gamma, results[,2]),
                     quantile(results[,2], 0.025)[[1]],
                     quantile(results[,2], 0.975)[[1]]
                     )
  
  # Spearman's footrule
  corr_stat[2,] <- c("Spearman's footrule",
                     conc_true$spm_foot,
                     mean(results[,3]),
                     var(results[,3]),
                     mse(conc_true$spm_foot, results[,3]),
                     quantile(results[,3], 0.025)[[1]],
                     quantile(results[,3], 0.975)[[1]]
                     )
  
  # Spearman's rho
  corr_stat[3,] <- c("Spearman's rho",
                     conc_true$spm_rho,
                     mean(results[,4]),
                     var(results[,4]),
                     mse(conc_true$spm_rho, results[,4]),
                     quantile(results[,4], 0.025)[[1]],
                     quantile(results[,4], 0.975)[[1]]
                     )
  
  return ( corr_stat )
}

####################
# BOUND ESTIMATORS
####################

bnd_bzid_simulation_main <- function(N, sim, p1, p2) {
  # Ensure N and sim and alpha are numeric
  if (!is.numeric(N) || !is.numeric(sim) || !is.numeric(alpha)) {
    stop("N, sim and alpha must be numeric.")
  }
  
  # Ensure p1 and p2 are numeric
  if (!is.numeric(p1) || !is.numeric(p2)) {
    stop("p1 and p2 must be numeric.")
  }
  
  # Create data frames
  results <- data.frame(matrix(nrow=0, ncol=7))
  colnames(results) <- c("Simulation",
                         "gamma_lw", "gamma_up",
                         "foot_lw", "foot_up",
                         "rho_lw", "rho_up")
  
  # Perform simulation
  for (ix in 1:sim) {
    # Generate sample
    X <- gen_sample(N, 1, p1, p2)
    
    # Assume p1 <= p2
    if (sum(X[,2] == 0) < sum(X[,1] == 0))
      X <- X[,c(2, 1)]
    
    # Estimate inflation
    p1.est <- sum(X[,1] == 0)/N
    p2.est <- sum(X[,2] == 0)/N
    
    # Compute estimates
    results[ix,] <- c(ix,
                      gini_gamma.bound_est(X[,1], X[,2])$lw,
                      gini_gamma.bound_est(X[,1], X[,2])$up,
                      spm_foot.bound_est(X[,1], X[,2])$lw,
                      spm_foot.bound_est(X[,1], X[,2])$up,
                      spm_rho.bound_est(X[,1], X[,2])$lw,
                      spm_rho.bound_est(X[,1], X[,2])$up)
    
  }
  
  return ( results )
}

bnd_bzid_simulation_stats <- function(results, p1, p2) {
  # Compute true values
  bnd_true <- list(gini_gamma.lw = gini_gamma.bound(p1,p2)$lw,
                    gini_gamma.up = gini_gamma.bound(p1,p2)$up,
                    spm_foot.lw = spm_foot.bound(p1,p2)$lw,
                    spm_foot.up = spm_foot.bound(p1,p2)$up,
                    spm_rho.lw = spm_rho.bound(p1,p2)$lw,
                    spm_rho.up = spm_rho.bound(p1,p2)$up)
  
  # Statistics
  bnd_stat <- data.frame(matrix(nrow=3, ncol=7))
  colnames(bnd_stat) <- c("Variable", "True lower", "True upper", "Mean lower", "Mean upper",
                          "MSE lower", "MSE upper")
  
  # Gini's gamma
  bnd_stat[1,] <- c("Gini's gamma",
                     bnd_true$gini_gamma.lw,
                     bnd_true$gini_gamma.up,
                     mean(results$gamma_lw),
                     mean(results$gamma_up),
                     mse(bnd_true$gini_gamma.lw, results$gamma_lw),
                     mse(bnd_true$gini_gamma.up, results$gamma_up)
  )
  
  # Spearman's footrule
  bnd_stat[2,] <- c("Spearman's footrule",
                     bnd_true$spm_foot.lw,
                     bnd_true$spm_foot.up,
                     mean(results$foot_lw),
                     mean(results$foot_up),
                     mse(bnd_true$spm_foot.lw, results$foot_lw),
                     mse(bnd_true$spm_foot.up, results$foot_up)
  )
  
  # Spearman's rho
  bnd_stat[3,] <- c("Spearman's rho",
                     bnd_true$spm_rho.lw,
                     bnd_true$spm_rho.up,
                     mean(results$rho_lw),
                     mean(results$rho_up),
                     mse(bnd_true$spm_rho.lw, results$rho_lw),
                     mse(bnd_true$spm_rho.up, results$rho_up)
  )
  
  return ( bnd_stat )
}
