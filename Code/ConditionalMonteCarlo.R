##########################################################
# Conditional Monte Carlo                   
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
# Estimate P(X+Y > a) where X and Y are independent
# exponential random variables with mean 1
# Exact probability exp(-a) * (1 + a)
set.seed(314)
a = 5
p_true = exp(-a) * (1 + a)

# We consider m independent Monte Carlo experiments
# each one of size n
m = 100
n = 500
X = matrix(rexp(n * m), nrow = m, byrow = TRUE)
Y = matrix(rexp(n * m), nrow = m, byrow = TRUE)

# We count the number of draws that are in the
# objective region and create the matrix of point 
# estimates
Z_MC = (X + Y > a)
theta_MC = t(apply(Z_MC, 1, cummean))

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
    yintercept = p_true,
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )
##########################################################

##########################################################
# Conditional Monte Carlo

# Reuse X to make a paired comparison and obtain the
# conditional Monte Carlo estimate
Z_CMC = exp(-pmax(a - X, 0))
theta_CMC = t(apply(Z_CMC, 1, cummean))

# Point estimate (the first row) and empirical range
df_theta_CMC = data.frame(
  n = 1:n,
  estimate = theta_CMC[1, ],
  lower = apply(theta_CMC, 2, min),
  upper = apply(theta_CMC, 2, max)
)

ggplot(df_theta_CMC, aes(x = n)) +
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
    yintercept = p_true,
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )

# Comparison of empirical ranges
ggplot() +
  geom_ribbon(
    data = df_theta_MC,
    aes(x = n, ymin = lower, ymax = upper,
        fill = "Ordinary MC"),
    alpha = 0.45
  ) +
  geom_ribbon(
    data = df_theta_CMC,
    aes(x = n, ymin = lower, ymax = upper,
        fill = "Conditional MC"),
    alpha = 0.45
  ) +
  geom_line(
    data = df_theta_CMC,
    aes(x = n, y = estimate),
    color = "steelblue4"
  ) +
  geom_hline(
    yintercept = p_true,
    color = "darkred"
  ) +
  scale_fill_manual(
    name = "Empirical range",
    values = c(
      "Ordinary MC" = "grey65",
      "Conditional MC" = "steelblue"
    )
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
)

# Empirical variance comparison
var_MC = var(theta_MC[, n])
var_CMC = var(theta_CMC[, n])

c(MC = var_MC,CMC = var_CMC,
  Variance_reduction = var_MC / var_CMC
)
##########################################################
