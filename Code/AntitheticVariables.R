##########################################################
# Antithetic variables                  
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Stochastic Simulation 
##########################################################

##########################################################
# Libraries
library(ggplot2)      # Version 4.0.2
library(ggthemes)     # Version 5.2.0
library(dplyr)        # Version 1.2.0
##########################################################

##########################################################
# Estimate integral of exp(x) from 0 to 1 using independent
# uniform (0,1) distributions
# exact probability: exp(1)-1
set.seed(314)
I_true = exp(1) - 1

# We consider m independent Monte Carlo experiments
# each one of size n
m = 100
n = 500

# Ordinary Monte Carlo
U = matrix(runif(n * m),nrow = m, ncol = n, byrow = TRUE)
theta_MC = t(apply(exp(U), 1, cummean))

# Point estimate (the first row) and empirical range
df_theta_MC = data.frame(
  n = 1:n,
  estimate = theta_MC[1, ],
  lower = apply(theta_MC, 2, min),
  upper = apply(theta_MC, 2, max)
)

ggplot(df_theta_MC, aes(x = n)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "lightblue",
    alpha = 0.75
  ) +
  geom_line(
    aes(y = estimate),
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = I_true,
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )
##########################################################

##########################################################
# Antithetic variates

# Generate n independent uniforms per experiment
U_AV = matrix(runif(n*m),nrow=m, ncol=n, byrow=TRUE)

# Construct antithetic pairs
Y = exp(U_AV)
Z = exp(1-U_AV)

# Average within each pair
W = (Y+Z)/2

# Sequential antithetic estimates
theta_AV = t(apply(W, 1, cummean))

# Point estimate and empirical range
df_theta_AV = data.frame(
  n = 1:n,
  estimate = theta_AV[1, ],
  lower = apply(theta_AV, 2, min),
  upper = apply(theta_AV, 2, max)
)

# Plot antithetic variates
ggplot(df_theta_AV, aes(x=n)) +
  geom_ribbon(
    aes(ymin=lower, ymax=upper),
    fill="lightblue",
    alpha=0.75
  ) +
  geom_line(
    aes(y=estimate),
    linewidth=0.7
  ) +
  geom_hline(
    yintercept=I_true,
    col="darkred"
  ) +
  theme_minimal() +
  labs(
    x="Number of independent antithetic pairs",
    y="Monte Carlo estimate"
  )

# Comparison of empirical ranges
ggplot() +
  geom_ribbon(
    data=df_theta_MC,
    aes(x=n, ymin=lower, ymax=upper,
        fill="Ordinary MC"),
    alpha=0.45
  ) +
  geom_ribbon(
    data=df_theta_AV,
    aes(x=n, ymin=lower, ymax=upper,
        fill="Antithetic variates"),
    alpha=0.45
  ) +
  geom_line(
    data=df_theta_AV,
    aes(x=n, y=estimate),
    color="steelblue4",
    linewidth=0.7
  ) +
  geom_hline(
    yintercept=I_true,
    color="darkred"
  ) +
  scale_fill_manual(
    name="Empirical range",
    values=c(
      "Ordinary MC"="grey65",
      "Antithetic variates"="steelblue"
    )
  ) +
  theme_minimal() +
  theme(legend.position="top") +
  labs(
    x="Number of independent simulation units",
    y="Monte Carlo estimate"
  )

# Empirical variance comparison


# Variance across the 100 independent experiments
var_MC = apply(theta_MC, 2, var)
var_AV = apply(theta_AV, 2, var)

df_variance = data.frame(
  n = 1:n,
  MC = var_MC,
  AV = var_AV
)

ggplot(df_variance, aes(x=n)) +
  geom_line(
    aes(y=MC, color="Ordinary MC"),
    linewidth=0.8
  ) +
  geom_line(
    aes(y=AV, color="Antithetic variates"),
    linewidth=0.8
  ) +
  scale_color_manual(
    name="Method",
    values=c(
      "Ordinary MC"="grey50",
      "Antithetic variates"="steelblue4"
    )
  ) +
  theme_minimal() +
  theme(legend.position="top") +
  labs(
    x="Number of independent simulation units",
    y="Empirical variance"
  )

# Empirical variance comparison
var_MC = var(theta_MC[, n])
var_AV = var(theta_AV[, n])

c(
  MC = var_MC,
  AV = var_AV,
  Variance_reduction = var_MC / var_AV
)
#######################################################