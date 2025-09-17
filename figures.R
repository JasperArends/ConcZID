########################################
# FIGURES
########################################
# The code snippets below can be used to reproduce the recreate the figures
# from the paper.

# Load libraries
library(ggplot2)
library(pals)

########################################
# BOUNDS
########################################

# Create a grid over [0, 1] x [0, 1]
N <- 200 # Level of detail
p1 <- seq(0, 1, length.out=N)
p2 <- seq(0, 1, length.out=N)
grid <- expand.grid(p1=p1, p2=p2)


# GINI'S GAMMA

# Compute bounds
grid$gini_gamma.up <- sapply(1:nrow(grid), function(ix) gini_gamma.bound(grid$p1[ix], grid$p2[ix])$up)
grid$gini_gamma.lw <- sapply(1:nrow(grid), function(ix) gini_gamma.bound(grid$p1[ix], grid$p2[ix])$lw)

# Create plots
ggplot(grid, aes(x = p1, y = p2)) +
  geom_raster(aes(fill = gini_gamma.lw)) +
  geom_contour(aes(z = gini_gamma.lw), color = "tan", bins = 30, linewidth = 0.5) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    name = expression(gamma[min]),
    limits = c(-1, 0)
  ) +
  labs(
    title = expression("Numerical evaluation of " * gamma[min]),
    x = expression(p[1]),
    y = expression(p[2])
  ) +
  coord_fixed() +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  theme(plot.margin=margin(t=10, r=10),
        plot.title = element_text(size=11, hjust=.5),
        axis.title.x = element_text(hjust=1),
        axis.title.y = element_text(hjust=1, angle=0),
        axis.text = element_text(colour="black"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2),
        panel.border = element_rect(linewidth=.2, fill=NA))

ggplot(grid, aes(x = p1, y = p2)) +
  geom_raster(aes(fill = gini_gamma.up)) +
  geom_contour(aes(z = gini_gamma.up), color = "tan", bins = 30, linewidth = 0.5) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    name = expression(gamma[max]),
    limits = c(0, 1)
  ) +
  labs(
    title = expression("Numerical evaluation of " * gamma[max]),
    x = expression(p[1]),
    y = expression(p[2])
  ) +
  coord_fixed() +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  theme(plot.margin=margin(t=10, r=10),
        plot.title = element_text(size=11, hjust=.5),
        axis.title.x = element_text(hjust=1),
        axis.title.y = element_text(hjust=1, angle=0),
        axis.text = element_text(colour="black"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2),
        panel.border = element_rect(linewidth=.2, fill=NA))


# SPEARMAN'S FOOTRULE

# Compute bounds
grid$spm_foot.up <- sapply(1:nrow(grid), function(ix) spm_foot.bound(grid$p1[ix], grid$p2[ix])$up)
grid$spm_foot.lw <- sapply(1:nrow(grid), function(ix) spm_foot.bound(grid$p1[ix], grid$p2[ix])$lw)

# Create plots
ggplot(grid, aes(x = p1, y = p2)) +
  geom_raster(aes(fill = spm_foot.lw)) +
  geom_contour(aes(z = spm_foot.lw), color = "tan", bins = 30, linewidth = 0.5) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    name = expression(phi[min]),
    limits = c(-.5, 0),
    labels=function(x) paste0(format(round(x, 2), nsmall = 2))
  ) +
  labs(
    title = expression("Numerical evaluation of " * phi[min]),
    x = expression(p[1]),
    y = expression(p[2])
  ) +
  coord_fixed() +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  theme(plot.margin=margin(t=10, r=10),
        plot.title = element_text(size=11, hjust=.5),
        axis.title.x = element_text(hjust=1),
        axis.title.y = element_text(hjust=1, angle=0),
        axis.text = element_text(colour="black"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2),
        panel.border = element_rect(linewidth=.2, fill=NA))

ggplot(grid, aes(x = p1, y = p2)) +
  geom_raster(aes(fill = spm_foot.up)) +
  geom_contour(aes(z = spm_foot.up), color = "tan", bins = 30, linewidth = 0.5) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    name = expression(phi[max]),
    limits = c(0, 1)
  ) +
  labs(
    title = expression("Numerical evaluation of " * phi[max]),
    x = expression(p[1]),
    y = expression(p[2])
  ) +
  coord_fixed() +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  theme(plot.margin=margin(t=10, r=10),
        plot.title = element_text(size=11, hjust=.5),
        axis.title.x = element_text(hjust=1),
        axis.title.y = element_text(hjust=1, angle=0),
        axis.text = element_text(colour="black"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2),
        panel.border = element_rect(linewidth=.2, fill=NA))

# SPEARMAN'S RHO

# Compute bounds
grid$spm_rho.up <- sapply(1:nrow(grid), function(ix) spm_rho.bound(grid$p1[ix], grid$p2[ix])$up)
grid$spm_rho.lw <- sapply(1:nrow(grid), function(ix) spm_rho.bound(grid$p1[ix], grid$p2[ix])$lw)

# Create plots
ggplot(grid, aes(x = p1, y = p2)) +
  geom_raster(aes(fill = spm_rho.lw)) +
  geom_contour(aes(z = spm_rho.lw), color = "tan", bins = 30, linewidth = 0.5) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    name = expression(rho[min]),
    limits = c(-1, 0)
  ) +
  labs(
    title = expression("Numerical evaluation of " * rho[min]),
    x = expression(p[1]),
    y = expression(p[2])
  ) +
  coord_fixed() +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  theme(plot.margin=margin(t=10, r=10),
        plot.title = element_text(size=11, hjust=.5),
        axis.title.x = element_text(hjust=1),
        axis.title.y = element_text(hjust=1, angle=0),
        axis.text = element_text(colour="black"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2),
        panel.border = element_rect(linewidth=.2, fill=NA))

ggplot(grid, aes(x = p1, y = p2)) +
  geom_raster(aes(fill = spm_rho.up)) +
  geom_contour(aes(z = spm_rho.up), color = "tan", bins = 30, linewidth = 0.5) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    name = expression(rho[max]),
    limits = c(0, 1)
  ) +
  labs(
    title = expression("Numerical evaluation of " * rho[max]),
    x = expression(p[1]),
    y = expression(p[2])
  ) +
  coord_fixed() +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  theme(plot.margin=margin(t=10, r=10),
        plot.title = element_text(size=11, hjust=.5),
        axis.title.x = element_text(hjust=1),
        axis.title.y = element_text(hjust=1, angle=0),
        axis.text = element_text(colour="black"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2),
        panel.border = element_rect(linewidth=.2, fill=NA))


########################################
# ESTIMATOR BOXPLOTS
########################################
# The following snippet can be used to recreate the boxplots from the simulation
# study. The plot below visualizes an example of simulated data, the reported
# figures in the manuscript contain the data from the actual simulation study in
# 'run_simulation.R'.

# Load dependencies
source("simulation_functions.R")
source("main_estimators.R")

# Parameters
p1 <- 0.2
p2 <- 0.2
alpha <- 0.5 # Copula parameter
N <- 150 # Sample size
sim <- 1000 # Number of simulations

# Perform simulations
df.sim <- corr_bzid_simulation_main(N, sim, alpha, p1, p2)

# True values
true_val <- data.frame(variable=c("Gini_gamma", "Spm_foot", "Spm_rho"),
                       value=c(gini_gamma.bound(p1,p2)$up * alpha,
                               spm_foot.bound(p1,p2)$up * alpha,
                               spm_rho.bound(p1,p2)$up * alpha))

# Make plot
df.sim[,2:4] %>%
  reshape2::melt(id.var=c()) %>%
  mutate_at("variable", as.factor) %>%
  ggplot() +
  # Boxplot
  geom_boxplot(aes(x=variable, y=value, fill=variable), outlier.size=.01, linewidth=.2) +
  scale_x_discrete(breaks=c("Gini_gamma", "Spm_foot", "Spm_rho"),
                   labels=c(expression("Gini's " * gamma),
                            expression("Spearman's " * phi),
                            expression("Spearman's " * rho))) +
  scale_fill_manual(breaks=c("Gini_gamma", "Spm_foot", "Spm_rho"),
                    values=c("#b3cde3", "#fbb4ae", "#ccebc5")) +
  # True value reference lines
  geom_hline(aes(yintercept=value), true_val, linetype='dashed', linewidth=.2) +
  facet_grid(~variable, scale="free_x") +
  # Labels and titles
  ggtitle(expression(p[1] * "= 0.2, " * p[2] * " = 0.2")) +
  # Theme settings
  ylim(0, .8) +
  theme_classic() +
  theme(panel.spacing.x=unit(0, "lines"),
        strip.background=element_blank(),
        strip.text.x=element_blank(),
        legend.position="na",
        axis.title = element_blank(),
        axis.text = element_text(size=7),
        plot.title = element_text(size=8),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth=.2))
