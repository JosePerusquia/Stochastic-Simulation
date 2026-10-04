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

# Function that simulates the Bessel process
BesselProcess = function(m,t,n){
  
  # Check the parameters
  if(t<=0){
    stop('t must be positive')
  }
  
  if(m<=0||m%%1!=0){
    stop('m must be a positive integer')
  }
  
  if(n<=0||n%%1!=0){
    stop('n must be a positive integer')
  }
  
  # We sample the Brownian motions
  mu = rep(0,m)
  sigma = rep(1,m)
  Bt = matrix(nrow=m,ncol=n+1)
  
  for(i in 1:m){
    Bt[i,]=BrownianMotion(mu[i],sigma[i],t,n)$X
  }
  
  # Bessel process
  Rt = sqrt(colSums(Bt^2))
  Rt = data.frame(t=seq(0,t,length.out=n+1),
                  X=Rt)
  
  return(Rt)
  
}

# Brownian bridge
BrownianBridge = function(sigma,n){
  
  # Check the parameters
  if(sigma<=0){
    stop('sigma must be positive')
  }
  
  if(n<=0||n%%1!=0){
    stop('n must be a positive integer')
  }
  
  # Generate the Brownian motion in [0,1]
  Bt = BrownianMotion(0,sigma,1,n)
  t = Bt$t
  Xt = Bt$X-t*Bt$X
  
  # Data frame with the path
  df = data.frame(t=t,X=Xt)
  return(df)
}
############################################################

############################################################
# Sample path of standard Brownian motion
set.seed(314159)
mu = 0
sigma = 1
t  = 10
n = 1000
Bt = BrownianMotion(mu,sigma,t,n)

# Plot the path
ggplot(data=Bt,aes(x=t,y=X))+
  geom_line()+
  geom_hline(yintercept = 0,linetype='dashed',
             col='darkred')+
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
set.seed(31415)
mu = 5
sigma = 1
t  = 1
n = 1000
Bt = BrownianMotion(mu,sigma,t,n)

# Plot the path
ggplot(data=Bt,aes(x=t,y=X))+
  geom_line()+
  geom_abline(slope = mu,intercept = 0,col='darkred')+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################

############################################################
# Bivariate standard Brownian motion
set.seed(314159)
mu = c(0,0)
sigma = c(1,1)
t  = 10
n = 1000
Bt1 = BrownianMotion(mu[1],sigma[1],t,n)
Bt2 = BrownianMotion(mu[2],sigma[2],t,n)

# Plot the path
Bt=data.frame(x=Bt1$X,y=Bt2$X)

ggplot(data=Bt,aes(x=x,y=y))+
  geom_path()+
  geom_point(x=0,y=0,col='red')+
  geom_point(x=Bt$x[n+1],y=Bt$y[n+1],col='red')+
  labs(x=expression(B[1](t)),y=expression(B[2](t)))+
  theme_minimal()
############################################################

############################################################
# Bivariate Brownian motion with different sigmas
set.seed(314159)
mu = c(0,0)
sigma = c(.5,5)
t  = 10
n = 1000
Bt1 = BrownianMotion(mu[1],sigma[1],t,n)
Bt2 = BrownianMotion(mu[2],sigma[2],t,n)

# Plot the path
Bt=data.frame(x=Bt1$X,y=Bt2$X)

ggplot(data=Bt,aes(x=x,y=y))+
  geom_path()+
  geom_point(x=0,y=0,col='red')+
  geom_point(x=Bt$x[n+1],y=Bt$y[n+1],col='red')+
  labs(x=expression(B[1](t)),y=expression(B[2](t)))+
  theme_minimal()
############################################################

############################################################
# Two dimensional Bessel process
set.seed(314159)
m = 2
t  = 10
n = 1000

# Bessel process
Rt = BesselProcess(m,t,n)

ggplot(data=Rt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(R[t]))+
  theme_minimal()
############################################################

############################################################
# Three dimensional Bessel process
set.seed(31415)
m = 3
t  = 10
n = 1000

# Bessel process
Rt = BesselProcess(m,t,n)

ggplot(data=Rt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(R[t]))+
  theme_minimal()
############################################################

############################################################
# Ten dimensional Bessel process
set.seed(3141)
m = 10
t  = 10
n = 1000

# Bessel process
Rt = BesselProcess(m,t,n)

ggplot(data=Rt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(R[t]))+
  theme_minimal()
############################################################

############################################################
# Brownian bridge
set.seed(314)
sigma = 1
n = 1000

Xt = BrownianBridge(sigma,n)

ggplot(data=Xt,aes(x=t,y=X))+
  geom_line()+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################
