############################################################
# Gaussian mixture models
# Author: Jose A. Perusquia Cortes
# Affil:  Facultad de Ciencias - UNAM
# Module: Stochastic Simulation 
############################################################

############################################################
# Libraries
library(ggplot2)         # Version 4.0.2
library(ggthemes)        # Version 5.2.0
library(dplyr)           # Version 1.2.0
############################################################

############################################################
# Function that creates standard normal random variables 
# using Box-Muller transformation
gmmSimulation = function(n,w,mu,sigmaSq){
  
  # Number of components
  p = length(w)
  
  # Check parameters
  if(any(w<0) | any(w>1) | sum(w)!=1){
    stop('The weights are not well-defined.')
  }
  
  if(any(sigmaSq<=0)){
    stop('Variances need to be positive')
  }
  
  # Check dimensiones
  if(length(mu)!=p | length(sigmaSq)!=p){
    stop('mu and sigmaSq need to be of the same length as w')
  }
  
  # Generate the latent allocations variables
  Z = sample(size=n,x=c(1:p),replace=T)
  
  # Generate the mixture model
  sigma = sqrt(sigmaSq)
  X = rnorm(n,mean = mu[Z], sd = sigma[Z])
  
  # Results
  df = data.frame(Z=Z,X=X)
  
  return(df)
}
############################################################

############################################################
# Sample
n = 1000
mu = c(-6,0,3)
sigma = c(3,1,2)
w = c(1/2,1/4,1/4)

# Sample from the mixture model
set.seed(314159)
gmm = gmmSimulation(n,w,mu,sigma)

# Theoretical density
x  = seq(-10,6,by=0.01)
fx = (dnorm(x,mu[1],sigma[1])+
        dnorm(x,mu[2],sigma[2])+
        dnorm(x,mu[3],sigma[3]))/3
df = data.frame(x,fx)

# Histogram with density
ggplot(data=gmm,aes(x=X,y=after_stat(density)))+
  geom_histogram(bins=20,col='black',fill='darkblue',
                 alpha= 0.3)+
  geom_line(data=df, aes(x=x,y=fx))+
  labs(x='',y='')+
  theme_minimal()
############################################################