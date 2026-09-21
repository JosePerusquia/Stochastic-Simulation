####################################################################
# Beta-binomial predictive distribution
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Stochastic simulation 
####################################################################

####################################################################
# Libraries
library(ggplot2)    # Version 4.0.2
library(ggthemes)   # Version 5.2.0
library(dplyr)      # Version 1.2.0
####################################################################

####################################################################
# Function that generates a sample of size N from the 
# beta-binomial predictive distribution using the mixture
# representation.
betaBinomialPredictive = function(N,n,m,alpha,beta,x){
  
  # Check on the parameters
  if(alpha <= 0){
    stop("alpha must be positive")
  }
  
  if(beta <= 0){
    stop("beta must be positive")
  }
  
  if(N <= 0 || N%%1 != 0){
    stop("N must be a positive integer")
  }
  
  if(n < 0 || n%%1 != 0){
    stop("n must be a non-negative integer")
  }
  
  if(m < 0 || m%%1 != 0){
    stop("m must be a non-negative integer")
  }
  
  if(x < 0 || x > n || x%%1 != 0){
    stop("x must be an integer between 0 and n")
  }
  
  # Generate the sample
  theta = rbeta(N,alpha+x,beta+n-x)
  Z = rbinom(N,m,theta)
  
  res = data.frame(Z)
  return(res)
}
####################################################################

####################################################################
# Parameters
N = 100000
n = 10
m = 15
alpha = 1
beta = 3
x = 7

# beta-binomial predictive
Z = betaBinomialPredictive(N,n,m,alpha,beta,x)
ggplot(data=Z,aes(x=Z))+
  geom_bar(col='black',fill='skyblue3')+
  theme_minimal()+
  labs(x="",y="")

# Theoretical probabilities vs empirical
a = alpha+x
b = beta+n-x

theoretical = data.frame(
  z = 0:m,
  prob = choose(m,0:m)*beta(0:m+a,m-(0:m)+b)/beta(a,b),
  source = "Theoretical"
)

empirical = Z%>%
  count(Z)%>%
  mutate(prob = n/N,source = "Empirical")%>%
  rename(z=Z)

probabilities = bind_rows(theoretical,empirical)

ggplot(data=probabilities,
       aes(x=z,y=prob,col=source,shape=source))+
  geom_point(size=2)+
  scale_color_manual(values=c("Theoretical"="blue",
                              "Empirical"="darkred"))+
  scale_shape_manual(values=c("Theoretical"=3,
                              "Empirical"=2))+
  theme_minimal()+
  labs(x="",y="",col="",shape="")

# Empirical mean and variance
mean(Z$Z)
var(Z$Z)

# Theoretical mean and variance
mu = m*((alpha+x)/(alpha+beta+n))
mu

sigma2 = (m^2)*((alpha+x)*(beta+n-x))/(((alpha+beta+n)^2)*(alpha+beta+n+1))+
  m*(beta(alpha+x+1,beta+n-x+1)/beta(alpha+x,beta+n-x))
sigma2
####################################################################