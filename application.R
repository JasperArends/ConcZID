########################################
# APPLICATION
########################################
# The code below illustrates how to apply the estimators of Gini's gamma, 
# Spearman's footrule and Spearman's rho to two data sets as investigated in
# the main manuscript.

# Load dependencies and libraries
source("main_estimators.R")
library(dplyr)

########################################
# DANISH INSURANCE DATASET
########################################

# Load data
library(fitdistrplus)
data(danishmulti)
danish <- data.frame("Contents" = danishmulti$Contents, "Profits"=danishmulti$Profits)

# Probabilities of inflation
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

########################################
# PRECIPITATION DATASET
########################################
# The data is freely available from the Royal Netherlands Meteorological
# Institute (KNMI) and can be downloaded at
# https://www.knmi.nl/nederland-nu/klimatologie/uurgegevens. The following part
# uses the precipitation (RH, given in 10^-1 mm) from Schiphol (240) and
# De Bilt (260) during the time period of 2021 until 2025.

# Load data
schiphol <- read.csv('./Data sets/Precipitation/uurgeg_240_2021-2030.txt', skip=30) %>%
  mutate(YYYYMMDD = as.Date(paste0(YYYYMMDD), format='%Y%m%d'))
debilt <- read.csv('./Data sets/Precipitation/uurgeg_260_2021-2030.txt', skip=30) %>%
  mutate(YYYYMMDD = as.Date(paste0(YYYYMMDD), format='%Y%m%d'))

# Filter data from 2020 - 2025 and 12H00 - 13H00
data <- data.frame(date = schiphol$YYYYMMDD,
                   hour = schiphol$HH,
                   schiphol = schiphol$RH,
                   debilt = debilt$RH) %>%
  filter(date < as.Date("20260101", format='%Y%m%d'), hour == 12,
         weekdays(date) == "Monday")

# Dataset entries with no rain are indicated by -1
data[data < 0] <- 0

# Probabilities of inflation
sum(data$schiphol == 0)/nrow(data)
sum(data$debilt == 0)/nrow(data)

# Add random noise to positive tied observations
data$schiphol[data$schiphol > 0] <- data$schiphol[data$schiphol > 0] +
  runif(sum(data$schiphol > 0))
data$debilt[data$debilt > 0] <- data$debilt[data$debilt > 0] +
  runif(sum(data$debilt > 0))


# Gini's gamma
gini_gamma.est(data$schiphol, data$debilt)
gini_gamma.bound_est(data$schiphol, data$debilt)

# Spearman's footrule
spm_foot.est(data$schiphol, data$debilt)
spm_foot.bound_est(data$schiphol, data$debilt)

# Spearman's rho
spm_rho.est(data$schiphol, data$debilt)
spm_rho.bound_est(data$schiphol, data$debilt)
