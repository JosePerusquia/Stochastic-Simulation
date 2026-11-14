##########################################################
# Control variate                    
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
# Control variate using Z=U

# Optimal coefficient
beta = 6 * (3 - exp(1))
Y_CV = exp(U) - beta * (U - 0.5)

# Point estimate and empirical range
theta_CV = t(apply(Y_CV, 1, cummean))

df_theta_CV = data.frame(
  n = 1:n,
  estimate = theta_CV[1, ],
  lower = apply(theta_CV, 2, min),
  upper = apply(theta_CV, 2, max)
)

ggplot(df_theta_CV, aes(x = n)) +
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

# Comparison of empirical ranges
ggplot() +
  geom_ribbon(
    data = df_theta_MC,
    aes(x = n, ymin = lower, ymax = upper,
        fill = "Ordinary MC"),
    alpha = 0.45
  ) +
  geom_ribbon(
    data = df_theta_CV,
    aes(x = n, ymin = lower, ymax = upper,
        fill = "Control variates"),
    alpha = 0.45
  ) +
  geom_line(
    data = df_theta_CV,
    aes(x = n, y = estimate),
    color = "steelblue4"
  ) +
  geom_hline(
    yintercept = I_true,
    color = "darkred"
  ) +
  scale_fill_manual(
    name = "Empirical range",
    values = c(
      "Ordinary MC" = "grey65",
      "Control variates" = "steelblue"
    )
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(
    x = "Number of simulations",
    y = "Estimated integral"
  )

# Empirical variance comparison
var_MC = var(theta_MC[, n])
var_CV = var(theta_CV[, n])

c(
  MC = var_MC,
  CV = var_CV,
  Variance_reduction = var_MC / var_CV
)
##########################################################