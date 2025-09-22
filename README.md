# Rank-based concordance measures for zero-inflated continuous distributions
Estimators of Gini's gamma, Spearman's footrule and Spearman's rho for zero-inflated continuous distributions and their corresponding upper and lower bounds.
Based on Arends, J.R.M., Lyu, G., Mesfioui, M., Perrone, E., and Trufin, J. (2025). *Rank-based concordance for zero-inflated data: New representations, estimators, and sharp bounds.* In preparation.

## Required libraries
`dplyr` `fitdistrplus` `ggplot2` `Metrics` `pals`

## Files

**application.R**
Example of how the estimators of the measures and their bounds can be computed for the Danish insurance dataset in `fitdistrplus`.

**figures.R**
Code used to create the figures from the paper.

**main_estimators.R**
Estimators for Gini's gamma, Spearman's footrule and Spearman's rho for zero-inflated continuous distributions and their corresponding bounds.

**run_simulation.R**
Simulation study as analysed in the manuscript.

**simulation_functions.R**
Functions used to perform the simulation study.
