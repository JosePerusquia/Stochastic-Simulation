############################################################
# Poisson process
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
# Function that simulates a path of the Poisson process
PoisProcess = function(lambda,t){
  
  # Check the parameters
  if(lambda<=0){
    stop('lambda must be positive')
  }
  
  if(t<=0){
    stop('t must be positive')
  }
  
  # Simulate the number of events
  Xt = rpois(1,lambda*t)
  
  # Given the number of events we sample uniform times
  # and sort them
  times = sort(runif(Xt,0,t))
  
  # Return a data frame with times and events
  res = data.frame(
    t = c(0,times),
    N = c(0:Xt)
  )
  
  return(res)
}
############################################################

############################################################
# Sample path
set.seed(314159)
lambda = 1
t = 10
Nt = PoisProcess(lambda,t)


# Plot the path
segments = data.frame(
  x = Nt$t,
  xend = c(Nt$t[-1],t),
  y = Nt$N,
  yend = Nt$N
)

ggplot(data=Nt,aes(x=t,y=N))+
  geom_point()+
  geom_segment(data=segments,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(N[t]))+
  theme_minimal()

############################################################