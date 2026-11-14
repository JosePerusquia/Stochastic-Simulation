##########################################################
# Importance Sampling                    
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
# Objective: P(X >= a) for X ~ Cauchy(0, 1).
# Exact probability: 1/2 - atan(a)/pi.
set.seed(314)
a = 2
theta = .5 - atan(a)/pi

# We consider m independent Monte Carlo experiments
# each one of size n
m = 100
n = 1000
sims_MC = matrix(rcauchy(m*n),nrow=m,ncol=n,byrow=T)

# We count the number of Cauchy draws that are in the
# objective region and create the matrix of point 
# estimates
flag_MC = sims_MC>=2
theta_MC = t(apply(flag_MC,1,cummean))

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
    yintercept = theta,
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )

# We can analyse the number of draws that actually
# contribute to the Monte Carlo point estimate
in_region_MC = apply(flag_MC,1,sum)
df_in_region_MC = data.frame(
  x = in_region_MC
)

# We have a binomial distribution of parameters
# n and theta, which gives a mean of n*theta
ggplot(data=df_in_region_MC,aes(x))+
  geom_histogram(bins=10,col='black',fill='darkred')+
  labs(x='',y='')+
  theme_minimal()
##########################################################

##########################################################
# Monte Carlo estimate using importance sampling using
# as proposal Cauchy(2,1)
sims_IS = matrix(rcauchy(n*m,location = 2),nrow=m,
                 ncol=n,byrow=T)

# We count the number of Cauchy draws that are in the
# objective region and create the matrix of point 
# estimates
flag_IS = sims_IS>=a
theta_IS = matrix(nrow=m,ncol=n)
for(i in 1:m){
  weight = (1+(sims_IS[i,]-2)^2)/(1+sims_IS[i,]^2)
  theta_IS[i,]= cumsum(flag_IS[i,]*weight)/(1:n)
}

# Point estimate and empirical range
df_theta_IS = data.frame(
  n = 1:n,
  estimate = theta_IS[1, ],
  lower = apply(theta_IS, 2, min),
  upper = apply(theta_IS, 2, max)
)

ggplot(df_theta_IS, aes(x = n)) +
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
    yintercept = .5 - atan(2)/pi,
    col = "darkred"
  ) +
  theme_minimal() +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )

# Comparison between ordinary MC and IS
ggplot() +
  geom_ribbon(
    data = df_theta_MC,
    aes(x = n, ymin = lower, ymax = upper,
        fill = "Ordinary MC"),
    alpha = 0.55
  ) +
  geom_ribbon(
    data = df_theta_IS,
    aes(x = n, ymin = lower, ymax = upper,
        fill = "Importance sampling"),
    alpha = 0.40
  ) +
  geom_line(
    data = df_theta_IS,
    aes(x = n, y = estimate),
    color = "steelblue4",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = theta,
    color = "darkred"
  ) +
  scale_fill_manual(
    name = "Empirical range",
    values = c(
      "Ordinary MC" = "grey65",
      "Importance sampling" = "steelblue"
    ),
    breaks = c("Ordinary MC", "Importance sampling")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Number of simulations",
    y = "Monte Carlo estimate"
  )

# In region values of the proposal distribution
in_region_IS = apply(flag_IS,1,sum)
df_in_region_IS = data.frame(
  x = in_region_IS
)

# We have a binomial distribution of parameters
# n and 0.5, which gives a mean of n/2
ggplot(data=df_in_region_IS,aes(x))+
  geom_histogram(bins=10,col='black',fill='darkred')+
  labs(x='',y='')+
  theme_minimal()
##########################################################