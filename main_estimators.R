########################################
# ESTIMATORS AND BOUNDS OF CONCORDANCE
# MEASURES FOR ZERO-INFLATED DATA
########################################


########################################
# MAIN ESTIMATOR FUNCTIONS
########################################

####################
# SPEARMAN'S
# FOOTRULE
####################
spm_foot.est <- function(X, Y) {
  Q.HM <- QM.est(X, Y)
  Q.MM <- (1 - max(sum(X==0)/length(X),sum(Y==0)/length(Y))^2)
  
  return ( (3 * Q.HM - Q.MM)/2 )
}

####################
# GINI'S GAMMA
####################
gini_gamma.est <- function(X, Y) {
  Q.HM <- QM.est(X, Y)
  Q.HW <- QW.est(X, Y)
  
  return ( Q.HM + Q.HW )
}

####################
# SPEARMAN'S RHO
####################
spm_rho.est <- function(X, Y) {
  # Ensure X and Y are numeric vectors
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }
  
  # Number of observations
  N <- length(X)
  if (length(Y) != N) {
    stop("X and Y must have the same length.")
  }
  
  # Probability masses
  p00 <- sum(X == 0 & Y == 0) / N # P(X = 0, Y = 0)
  p01 <- sum(X == 0 & Y >  0) / N  # P(X = 0, Y > 0)
  p10 <- sum(X >  0 & Y == 0) / N  # P(X > 0, Y = 0)
  p11 <- sum(X >  0 & Y >  0) / N   # P(x > 0, Y > 0)
  
  ######################################
  # Split data
  D10 <- which(X >  0 & Y == 0)
  D01 <- which(X == 0 & Y >  0)
  D11 <- which(X >  0 & Y >  0)
  
  ######################################
  # Compute probabilities
  p1s <- ifelse(length(D10) > 0 & length(D11) > 0,  # P[X10 > X11]
                mean(outer(X[D10], X[D11], FUN = ">")),
                0)
  p2s <- ifelse(length(D01) > 0 & length(D11) > 0,  # P[Y01 > Y11]
                mean(outer(Y[D01], Y[D11], FUN = ">")),
                0)
  
  ######################################
  # Spearman's rho for positive values
  rho11 <- ifelse(length(D11) > 0, cor(X[D11], Y[D11], method="spearman"), 0)
  rho10 <- ifelse(length(D10) > 0 & length(D11) > 0,
                  12 * mean( colMeans(outer(X[D10], X[D11], FUN = "<")) *
                               colMeans(outer(Y[D11], Y[D11], FUN = "<"))) - 6 * (1 - p1s),
                  0)
  rho01 <- ifelse(length(D01) > 0 & length(D11) > 0,
                  12 * mean( colMeans(outer(X[D11], X[D11], FUN = "<")) *
                               colMeans(outer(Y[D01], Y[D11], FUN = "<"))) - 6 * (1 - p2s),
                  0)
  rho00 <- ifelse(length(D10) > 0 & length(D01) > 0,
                  12 * mean( colMeans(outer(X[D10], X[D11], FUN = ">")) *
                               colMeans(outer(Y[D01], Y[D11], FUN = ">")) ) + 3 - 6 * (p1s + p2s),
                  0)
  
  ######################################
  # Estimate Spearman's rho
  rho <- p11^3 * rho11 + p11^2 * p10 * rho10 + p11^2 * p01 * rho01 + p11 * p10 * p01 * rho00 +
    3 * (p11 * p00 - p10 * p01) + 3 * p11 * (p10 * (1 - 2 * p1s) + p01 * (1 - 2 * p2s))
  
  return ( rho )
}

########################################
# CONCORDANCE W.R.T. UPPER
# FRÉCHET-HOEFFDING COPULA
########################################
QM.est <- function(X, Y) {
  # Ensure X and Y are numeric vectors
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }
  
  N <- length(X)
  if (length(Y) != N) {
    stop("X and Y must have the same length.")
  }

  # Assume p1 <= p2
  if (sum(X == 0) > sum(Y == 0)) {
    Z <- X
    X <- Y
    Y <- Z
  }
  
  # Probability masses
  p00 <- sum(X == 0 & Y == 0)/N
  p10 <- sum(X >  0 & Y == 0)/N
  p01 <- sum(X == 0 & Y >  0)/N
  p11 <- sum(X >  0 & Y >  0)/N
  
  p1 <- sum(X == 0)/N
  p2 <- sum(Y == 0)/N # Assume p2 >= p1
  
  # Divide into p11(I), p11(II), p10(I) and p01(I)
  Q11I <- (X > 0 & rank(X, ties.method="max")/N <= p2 & Y >  0) # X > 0, Y > 0 and F(X) <= p2
  Q11II <- (X > 0 & rank(X, ties.method="max")/N > p2 & Y >  0) # X > 0, Y > 0 and F(X) >  p2
  Q10I <- (X > 0 & rank(X, ties.method="max")/N <= p2 & Y == 0) # X > 0, Y = 0 and F(X) <= p2
  
  # Normalised ranks of X | X > 0 and Y | Y > 0 (we simply set R(0) = 0)
  Rp <- rep(0, N) # Initialise ranks
  bf <- rank(X, ties.method="max")/N > p2
  Rp[bf] <- (rank(X[bf])/sum(bf)) # Ranks of X | F(X) > p2
  Sp <- rep(0, N) # Initialise ranks
  Sp[Y > 0] <- rank(Y[Y > 0])/sum(Y>0) # Ranks of Y | Y > 0
  
  Rp2 <- rep(0, N) # Initialise ranks
  bf <- X > 0 & rank(X, ties.method="max")/N <= p2 # Filter X > 0 and F(X) < p2
  Rp2[bf] <- rank(X[bf])/sum(bf)
  
  # Relevant statistics
  stat <- c(mean(pmin(Rp[Q11II], Sp[Q11II])),
            mean(Rp2[Q11I]))
  
  stat[is.na(stat)] <- 0
  
  # Compute concordance
  QHM <- (1 - p2) * (sum(Q11II)/N * (4 * stat[1] - 1) + sum(Q10I)/N) +
    (p2 - p1) * (sum(Q11I)/N * (2 * stat[2] - 1) + sum(Q11II)/N) +
    p00 * (1 - p2) + p11 * p1 - p01 * (p2 - p1)
  
  return ( QHM )
}

########################################
# CONCORDANCE W.R.T. LOWER
# FRÉCHET-HOEFFDING COPULA
########################################
QW.est <- function(X, Y) {
  # Ensure X and Y are numeric vectors
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }
  
  N <- length(X)
  if (length(Y) != N) {
    stop("X and Y must have the same length.")
  }

  # Assume p1 <= p2
  if (sum(X == 0) > sum(Y == 0)) {
    Z <- X
    X <- Y
    Y <- Z
  }
  
  # Probability masses
  p00 <- sum(X == 0 & Y == 0)/N
  p10 <- sum(X >  0 & Y == 0)/N
  p01 <- sum(X == 0 & Y >  0)/N
  p11 <- sum(X >  0 & Y >  0)/N
  
  p1 <- sum(X == 0)/N
  p2 <- sum(Y == 0)/N
  
  if ( p1 + p2 < 1 ) {
    
    QI <- X > 0 & rank(X, ties.method="max")/N <= 1 - p2
    QII <- rank(X, ties.method="max")/N > 1 - p2
    QA <- rank(Y, ties.method="max")/N > 1 - p1
    QB <- Y > 0 & rank(Y, ties.method="max")/N <= 1 - p1
    
    Q11I <- QI & Y > 0
    Q11II <- QII & Y > 0
    Q11A <- QA & X > 0
    Q11B <- QB & X > 0
    
    Q11BI <- Q11B & Q11I
    
    # Estimate margins
    #   Ranks of X
    RpI <- rep(0, N) # R(X | X > 0 and F(X) <= 1 - p2)
    RpI[QI] <- rank(X[QI])/sum(QI)
    RpII <- rep(0, N) # R(X | F(X) > 1 - p2)
    RpII[QII] <- rank(X[QII])/sum(QII)
    
    #   Ranks of Y
    SpA <- rep(0, N) # R(Y | Y > 0 and G(Y) <= 1 - p1)
    SpA[QA] <- rank(Y[QA])/sum(QA)
    SpB <- rep(0, N) # R(Y | G(Y) > 1 - p1)
    SpB[QB] <- rank(Y[QB])/sum(QB)
    
    # Compute relevant statistics
    stat <- c(mean(pmax(RpI[Q11BI] + SpB[Q11BI] - 1, 0)),
              mean(SpA[Q11A]), # P(Y1 > Y2 | C1101, G(Y1) > 1 - p1))
              mean(RpII[Q11II]), # P(X1 > X2 | C1110, F(X1) > 1 - p2)
              mean(SpB[Q11B & Q11II]), # P(Y1 > Y2 | C1111, F(X1) > 1 - p2, G(Y1) <= 1 - p1)
              mean(RpI[Q11A & Q11I]) ) # P(X1 > X2 | C1111, F(X1) <= 1 - p2, G(Y1) > 1 - p1)
    stat[is.na(stat)] <- 0
    
    # Compute concordance
    QHW <- (1 - p1 - p2) * (sum(Q11BI)/N * 4 * stat[1] - 1 + 4 * sum(Q11A & Q11II)/N +
                              4 * sum(Q11A * Q11I)/N * stat[5] +
                              4 * sum(Q11B & Q11II)/N * stat[4]) +
      p1 * (sum(Q11A)/N * 2 * stat[2] - (1 - p1)) +
      p2 * (sum(Q11II)/N * 2 * stat[3] - (1 - p2))
  }
  
  else { # p1 + p2 > 1
    Q1m <- X > 0
    Qm1 <- Y > 0
    Q11 <- Q1m & Qm1
    
    # Compute margins
    Rp <- rep(0, N) # F(X) | X > 0
    Rp[Q1m] <- rank(X[Q1m])/sum(Q1m)
    Sp <- rep(0, N) # G(Y) | Y > 0
    Sp[Qm1] <- rank(Y[Qm1])/sum(Qm1)
    
    # Relevant statistics
    stat <- c(mean(Sp[Q11]), mean(Rp[Q11]))
    stat[is.na(stat)] <- 0
    
    QHW <- p11 * (1 - p1) * (2 * stat[2] - 1) +
      p11 * (1 - p2) * (2 * stat[1] - 1) +
      p11 * (p1 + p2 - 1) - p10 * (1 - p2) - p01 * (1 - p1)
  }
  
  return ( QHW )
}

########################################
# BOUNDS
########################################
gini_gamma.bound_est <- function(X, Y) {
  # Ensure X and Y are numeric vectors
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }
  
  # Estimate probabilities of inflation
  p1 <- sum(X == 0)/length(X) # P[X = 0]
  p2 <- sum(Y == 0)/length(Y) # P[Y = 0]
  
  # Upper bound
  gamma.up <- ifelse(p1 + p2 >= 1,
                     (1 - max(p1,p2)) * (p1 + p2 + max(p1,p2)),
                     (1 - max(p1,p2)^2))
  
  # Lower bound
  gamma.lw <- ifelse(p1 + p2 >= 1,
                     max(p1,p2) - max(p1,p2)^2 - 3 * (1 - p1)*(1 - p2),
                     (1 - p1 - p2)^2 - 2 * (1 - p1) * (1 - p2))
  
  return (list(up=gamma.up, lw=gammal.w))
}

spm_foot.bound_est <- function(X, Y) {
  # Ensure X and Y are numeric vectors
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }
  
  # Estimate probabilities of inflation
  p1 <- sum(X == 0)/length(X) # P[X = 0]
  p2 <- sum(Y == 0)/length(Y) # P[Y = 0]
  
  # Upper bound
  spm_foot.up <- 1 - max(p1,p2)^2
  
  # Lower bound
  spm_foot.lw <- ifelse(p1 + p2 >= 1,
                        (3 * max(p1,p2) - 2 * max(p1,p2)^2 - 3 * (1 - p1) * (1 - p2) - 1)/2,
                        - (1 - max(p1,p2)^2)/2)
  
  return (list(up=spm_foot.up, lw=spm_foot.lw))
}

spm_rho.bound_est <- function(X, Y) {
  # Ensure X and Y are numeric vectors
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }
  
  # Estimate probabilities of inflation
  p1 <- sum(X == 0)/length(X) # P[X = 0]
  p2 <- sum(Y == 0)/length(Y) # P[Y = 0]
  
  # Upper bound
  spm_rho.up <- 1 - max(p1,p2)^3
  
  # Lower bound
  spm_rho.lw <- ifelse(p1 + p2 <= 1,
                       p1^3 + p2^2 - 1,
                       -3 * (1 - p1) * (1 - p2))
  
  return (list(up=spm_rho.up, lw=spm_rho.lw))
}
