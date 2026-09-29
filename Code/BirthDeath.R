############################################################
# Birth-death process
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
# Function that simulates a path of the birth-death process
# with constante birth and death rates lambda, mu
BirthDeath = function(mu,lambda,t){
  
  # Check the parameters
  if(mu<=0){
    stop('mu must be positive')
  }
  
  if(lambda<=0){
    stop('lambda must be positive')
  }
  
  if(t<=0){
    stop('t must be positive')
  }
  
  # Exponential parameters
  q_0 = lambda
  q_i = mu+lambda
  
  # Transition probabilities
  prob_i = c(mu,lambda)/(mu+lambda)
  
  # Start the chain
  Xt = c(0)
  n = 1
  
  time = c(0)
  cumtime = sum(time)
  
  while(cumtime<t){
    
    if(Xt[n]==0){
      time_aux = rexp(1,q_0)
      jump = 1
    }else{
      time_aux = rexp(1,q_i)
      jump = sample(x=c(Xt[n]-1,Xt[n]+1),size = 1,
                    prob = prob_i)
    }
    
    Xt = c(Xt,jump)
    n = n+1
    time = c(time,time_aux)
    cumtime = sum(time)
  }
  
  time = cumsum(time)
  time[n] = t
  Xt[n] = Xt[n-1]
  
  
  # Return a data frame with times and events
  res = data.frame(
    t = time,
    X = Xt
  )
  
  return(res)
}
############################################################

############################################################
# Sample path with mu=lambda=1
set.seed(3141592)
t = 10
mu = 1
lambda = 1

Xt = BirthDeath(mu,lambda,t)
n = nrow(Xt)

# Plot the path of the jump process
segmentsJP = data.frame(
  x = Xt$t[-n],
  xend = c(Xt$t[-1]),
  y = Xt$X[-n],
  yend = Xt$X[-n]
)

ggplot(data=Xt,aes(x=t,y=X))+
  geom_point(data=Xt[-n,])+
  geom_segment(data=segmentsJP,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################

############################################################
# Sample path with mu=lambda=0.25
set.seed(3141592)
t = 20
mu = 0.25
lambda = 0.25

Xt = BirthDeath(mu,lambda,t)
n = nrow(Xt)

# Plot the path of the jump process
segmentsJP = data.frame(
  x = Xt$t[-n],
  xend = c(Xt$t[-1]),
  y = Xt$X[-n],
  yend = Xt$X[-n]
)

ggplot(data=Xt,aes(x=t,y=X))+
  geom_point(data=Xt[-n,])+
  geom_segment(data=segmentsJP,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################

############################################################
# Sample path with 1=mu<lambda=0.25
set.seed(3141592)
t = 20
mu = 1
lambda = 0.25

Xt = BirthDeath(mu,lambda,t)
n = nrow(Xt)

# Plot the path of the jump process
segmentsJP = data.frame(
  x = Xt$t[-n],
  xend = c(Xt$t[-1]),
  y = Xt$X[-n],
  yend = Xt$X[-n]
)

ggplot(data=Xt,aes(x=t,y=X))+
  geom_point(data=Xt[-n,])+
  geom_segment(data=segmentsJP,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################

############################################################
# Sample path with 0.25=mu>lambda=1
set.seed(3141592)
t = 20
mu = 0.25
lambda = 1

Xt = BirthDeath(mu,lambda,t)
n = nrow(Xt)

# Plot the path of the jump process
segmentsJP = data.frame(
  x = Xt$t[-n],
  xend = c(Xt$t[-1]),
  y = Xt$X[-n],
  yend = Xt$X[-n]
)

ggplot(data=Xt,aes(x=t,y=X))+
  geom_point(data=Xt[-n,])+
  geom_segment(data=segmentsJP,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(X[t]))+
  theme_minimal()
############################################################