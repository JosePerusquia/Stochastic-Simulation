############################################################
# Inhomogeneous Poisson process
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

# Function that simulates a path of the Poisson process
# with intensity lamda(t) = 1 + t in the interval [0,T]
IPoisProcess = function(Tmax){
  
  # Check the parameter
  if(Tmax<=0){
    stop('T must be positive')
  }
  
  # Generate the Poisson process
  M = 1+Tmax
  PP = PoisProcess(M,Tmax)
  
  # Candidates for the times
  times = PP$t[-1]
  
  # Acceptance
  u = runif(length(times))
  accepted = u<((1+times)/M)
  
  # Results
  res = data.frame(
    times=c(0,times),
    accepted = c(TRUE,accepted),
    N = PP$N
  )
  
  return(res)
  
}
############################################################

############################################################
# Sample path
set.seed(3141)
Tmax = 10
Nt = IPoisProcess(Tmax)

# Proportion accepted (the first point is always accepted
# so we don't take it into account)
(sum(Nt$accepted)-1)/(nrow(Nt)-1)

# Plot the path of the homogeneous Poisson process
segmentsPP = data.frame(
  x = Nt$times,
  xend = c(Nt$times[-1],Tmax),
  y = Nt$N,
  yend = Nt$N
)

ggplot(data=Nt,aes(x=times,y=N))+
  geom_point()+
  geom_segment(data=segmentsPP,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(N[t]))+
  theme_minimal()

# Thinning
ggplot(data=Nt[-1,],aes(x=times,y=N))+
       geom_point(aes(shape=accepted,
                      col = accepted),size=2)+
       labs(x=expression(t),y=expression(N[t]))+
       theme_minimal()

# Filter the accepted values
NtAccepted = Nt%>%
  filter(accepted)
NtAccepted$N = 0:(nrow(NtAccepted)-1)

# Plot the path
segmentsIPP = data.frame(
  x = NtAccepted$times,
  xend = c(NtAccepted$times[-1],Tmax),
  y = NtAccepted$N,
  yend = NtAccepted$N
)

ggplot(data=NtAccepted,aes(x=times,y=N))+
  geom_point()+
  geom_segment(data=segmentsIPP,
               aes(x=x,y=y,xend=xend,yend=yend))+
  labs(x=expression(t),y=expression(N[t]))+
  theme_minimal()
############################################################