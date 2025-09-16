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

# Spearman's rho
spm_rho.est(danish$Contents, danish$Profits)
spm_rho.bound_est(danish$Contents, danish$Profits)

# Spearman's footrule
spm_foot.est(danish$Contents, danish$Profits)
spm_foot.bound_est(danish$Contents, danish$Profits)

# Gini's gamma
gini_gamma.est(danish$Contents, danish$Profits)
gini_gamma.bound_est(danish$Contents, danish$Profits)

########################################
# PRECIPITATION DATASET
########################################
# The data is freely available from the Royal Netherlands Meteorological
# Institute (KNMI) and can be downloaded at
# https://www.knmi.nl/nederland-nu/klimatologie/uurgegevens. The following part
# uses the precipitation (RH) from Schiphol (240) and De Bilt (260) during the
# time period of 2020 until 2023.

# Load data
schiphol <- read.csv('uurgeg_240_2021-2030.txt', skip=30)
debilt <- read.csv('uurgeg_260_2021-2030.txt', skip=30)

# Filter data from 2020 - 2023 and 12H00 - 13H00
data <- data.frame(date = schiphol$YYYYMMDD,
                   hour = schiphol$HH,
                   schiphol = schiphol$RH,
                   debilt = debilt$RH) %>%
  filter(date < 20240000, hour == 12)

# Dataset entries with no rain are indicated by -1
data[data < 0] <- 0

# Probabilities of inflation
sum(data$schiphol == 0)/nrow(data)
sum(data$debilt == 0)/nrow(data)

# Add random noise to positive tied observations
data$schiphol[data$schiphol > 0] <- data$schiphol[data$schiphol > 0] 
  runif(sum(data$schiphol > 0))
data$debilt[data$debilt] <- data$debilt[data$debilt > 0] +
  runif(sum(data$debilt > 0))

# Spearman's rho
spm_rho.est(dt$schiphol, dt$debilt)
spm_rho.bound_est(dt$schiphol, dt$debilt)

# Spearman's footrule
spm_foot.est(dt$schiphol, dt$debilt)
spm_foot.bound_est(dt$schiphol, dt$debilt)

# Gini's gamma
gini_gamma.est(dt$schiphol, dt$debilt)
gini_gamma.bound_est(dt$schiphol, dt$debilt)


# Illustrative figures
ggplot(data/10, aes(x=debilt, y=schiphol)) +
  geom_point(size=.4) +
  labs(x="De Bilt", y="Schiphol") +
  coord_fixed(ratio = 1, xlim=c(0, 11), ylim=c(0, 5)) +
  scale_x_continuous(breaks=0:11, expand=expansion(mult = c(0.02, .05))) +
  scale_y_continuous(breaks=0:5, expand=expansion(mult = c(0.05, 0.05))) +
  theme_classic() +
  theme(axis.line = element_line(linewidth=.2),
        axis.ticks.x = element_line(linewidth=.2),
        axis.ticks.y = element_line(linewidth=.2))

ggplot(data/10, aes(x=debilt)) +
  geom_histogram(aes(y=after_stat(count) / nrow(data)), breaks=c(0, seq(0.1, 11, 0.1))) +
  labs(x="", y="Margin") +
  xlim(0, 11) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, .5, 1)) +
  theme_classic() +
  theme(axis.line = element_line(linewidth=.2),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=.2))

ggplot(data/10, aes(x=schiphol)) +
  geom_histogram(aes(y=after_stat(count) / nrow(data)), breaks=c(0, seq(0.1, 5, 0.1))) +
  labs(x="", y="Margin") +
  scale_y_continuous(position="right", limits=c(0, 1), breaks=c(0, .5, 1)) +
  scale_x_reverse() +
  theme_classic() +
  theme(axis.line = element_line(linewidth=.2),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=.2))
