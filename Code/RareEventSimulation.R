##########################################################
# Rare event simulation                    
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
# Objective: P(X >= a) for X ~ N(0, 1)
# Exact probability is approximately 0.00003167
set.seed(314)
a = 4
theta = pnorm(a,lower.tail = F)

# We consider m independent Monte Carlo experiments
# each one of size n
m = 100
n = 1000
sims_MC = matrix(rnorm(m*n),nrow=m,ncol=n,byrow=T)

# We count the number of Cauchy draws that are in the
# objective region and create the matrix of point 
# estimates
flag_MC = sims_MC>=a
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
# n and theta, which gives a mean of 0.03167124
ggplot(data=df_in_region_MC,aes(x))+
  geom_histogram(bins=10,col='black',fill='darkred')+
  labs(x='',y='')+
  theme_minimal()
##########################################################

##########################################################
# Monte Carlo estimate using importance sampling using
# as proposal Cauchy(2,1)
sims_IS = matrix(rnorm(n*m,mean = 4),nrow=m,
                 ncol=n,byrow=T)

# We count the number of Cauchy draws that are in the
# objective region and create the matrix of point 
# estimates
flag_IS = sims_IS>=a
theta_IS = matrix(nrow=m,ncol=n)
for(i in 1:m){
  weight = exp(-4*sims_IS[i,]+8)
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
    yintercept = theta,
    col = "darkred"
  ) +
  theme_minimal() +
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
# n and 0.5, which gives a mean of 500
ggplot(data=df_in_region_IS,aes(x))+
  geom_histogram(bins=10,col='black',fill='darkred')+
  labs(x='',y='')+
  theme_minimal()

# Plot of 20 weights
X=seq(0,8,by=.01)
W=exp(-4*X + 8)
df_w=data.frame(X,W)

X=rep(0,20)
Y=sort(rnorm(20,4))
W_y=exp(-4*Y + 8)
df_y=data.frame(X,Y,W_y)

ggplot(data=df_w,aes(x=X,y=W))+
  geom_line()+
  geom_point(data=df_y,aes(x=Y,y=W_y))+
  geom_segment(data=df_y,aes(x=Y,xend=Y,y=X,yend=W_y))+
  theme_minimal()+
  labs(x='',y='')
##########################################################