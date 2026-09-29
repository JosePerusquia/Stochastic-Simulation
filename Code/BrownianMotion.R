############################################################
# Brownian motion
# Author: Jose A. Perusquia Cortes
# Affil:  Facultad de Ciencias-UNAM
# Module: Stochastic Simulation 
############################################################

############################################################
# Libraries
library(ggplot2)         # Version 4.0.2
library(ggthemes)        # Version 5.2.0
library(dplyr)           # Version 1.2.0
############################################################

############################################################
# Function that simulates a path on [0,t] of the process
# Xt = mu*t + sigma*Bt, where Bt is a standard Brownian motion.
# Here mu is the drift and sigma is the diffusion coefficient.
# The interval [0,t] is divided into n equidistant subintervals.
BrownianMotion = function(mu,sigma,t,n){
  
  # Check the parameters
  if(t<=0){
    stop('t must be positive')
  }
  
  if(sigma<=0){
    stop('sigma must be positive')
  }
  
  if(n<=0||n%%1!=0){
    stop('n must be a positive integer')
  }
  
  # Times
  deltaT = t/n
  times = seq(0,t,length.out=n+1)
  
  # Normal random variables 
  Z = rnorm(n)
  
  # Sequentially obtain the stochastic process
  Xt = numeric(n+1)
  
  for(i in 2:(n+1)){
    Xt[i] = Xt[i-1] + mu*deltaT + sigma*sqrt(deltaT)*Z[i-1]
  }
  
  # Equivalent vectorised implementation using cumulative sums
  # increments = mu*deltaT + sigma*sqrt(deltaT)*Z
  # Xt = c(0,cumsum(increments))
  
  # Return a data frame with times and process values
  res = data.frame(
    t = times,
    X = Xt
  )
  
  return(res)
}
############################################################

############################################################
# Sample path of standard Brownian motion
set.seed(314159)
mu = 0
sigma = 1
t  = 1
n = 100
Bt = BrownianMotion(mu,sigma,t,n)

# Plot the path
ggplot(data=Bt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(B[t]))+
  theme_minimal()
############################################################

############################################################
# Sample path of Brownian motion with diffusion coefficient
set.seed(314159)
mu = 0
sigma = 2.5
t  = 1
n = 100
Bt = BrownianMotion(mu,sigma,t,n)

# Plot the path
ggplot(data=Bt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(B[t]))+
  theme_minimal()
############################################################

############################################################
# Sample path of Xt = mu*t + sigma*Bt
set.seed(314159)
mu = 5
sigma = 2.5
t  = 1
n = 100
Bt = BrownianMotion(mu,sigma,t,n)

# Plot the path
ggplot(data=Bt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################