#################################################################
# Monte Carlo Integration                    
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Stochastic Simulation 
#################################################################

#################################################################
# Libraries 
library(ggplot2)      # Version 4.0.2
library(ggthemes)     # Version 5.2.0
library(dplyr)        # Version 1.2.0
#################################################################

#################################################################
# Integrate exp(-x^2) from 0 to infinity

# Number of samples considered and size of each sample
n = 1000
m = 100

# Matrix of uniform random variables
set.seed(314)
sims = matrix(runif(n*m),nrow=m,ncol=n)

# Matrix of estimates for each n and m
theta = matrix(nrow=m,ncol=n)
for(i in 1:m){
  hy = exp(-((1/sims[i,])-1)^2)/(sims[i,]^2)
  theta[i,]= cumsum(hy)/(1:n)
}

# Point estimate and empirical range
df_theta = data.frame(
  n = 1:n,
  estimate = theta[1, ],
  lower = apply(theta, 2, min),
  upper = apply(theta, 2, max)
)

ggplot(df_theta, aes(x = n)) +
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
    yintercept = sqrt(pi) / 2,
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )
#################################################################

#################################################################
# Integrate x^{alpha-1}exp(-x) from 0 to infinity
# using exponential iid random variables
alpha = 2.9

# Number of samples considered and size of each sample
n = 500
m = 100

# Matrix of uniform random variables
set.seed(31415)
sims = matrix(rexp(n*m),nrow=m,ncol=n)

# Matrix of estimates for each n and m
theta = matrix(nrow=m,ncol=n)
for(i in 1:m){
    theta[i,]= cumsum(sims[i,]^(alpha-1))/(1:n)
}

# Point estimate and empirical range
df_theta = data.frame(
  n = 1:n,
  estimate = theta[1, ],
  lower = apply(theta, 2, min),
  upper = apply(theta, 2, max)
)

ggplot(df_theta, aes(x = n)) +
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
    yintercept = gamma(alpha),
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )
#################################################################