########################################
# APPLICATION EXAMPLE
########################################

# Load data
library(fitdistrplus)
data(danishmulti)
danish <- data.frame("Contents" = danishmulti$Contents, "Profits"=danishmulti$Profits)

# Estimate inflation
sum(danish$Contents == 0)/nrow(danish)
sum(danish$Profits == 0)/nrow(danish)

# Gini's gamma
gini_gamma.est(danish$Contents, danish$Profits)
gini_gamma.bound_est(danish$Contents, danish$Profits)

# Spearman's footrule
spm_foot.est(danish$Contents, danish$Profits)
spm_foot.bound_est(danish$Contents, danish$Profits)

# Spearman's rho
spm_rho.est(danish$Contents, danish$Profits)
spm_rho.bound_est(danish$Contents, danish$Profits)
